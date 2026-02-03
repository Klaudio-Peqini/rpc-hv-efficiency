#!/usr/bin/env python3
"""efficiency_vs_hv.py

Stage 2b (HV scan analyzer) in a clean pipeline:

  Stage 1:  tracklets_reconstruction_efficiency_v2.py -> eff_map_HV*.root
  Stage 2a: make_eff_summary.py                      -> summary_*.json
  Stage 2b: efficiency_vs_hv.py                      -> plots + sigmoid fit

This script reads one or more JSON summaries (produced by make_eff_summary.py),
builds efficiency vs HV arrays (ROI and/or global), plots them, and optionally
fits a sigmoid.

It intentionally does not reprocess the original TTrees.

Examples
  python3 efficiency_vs_hv.py --glob 'summary_*.json' --use roi --fit

  python3 efficiency_vs_hv.py --files s1.json s2.json s3.json --use global --fit

Outputs
  - PNG plot(s)
  - Optionally a JSON with fit parameters

"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time
from typing import Any, Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


def sigmoid(x: np.ndarray, plateau: float, x0: float, k: float, baseline: float) -> np.ndarray:
    """4-parameter logistic (baseline + (plateau-baseline)/(1+exp(-(x-x0)/k)))."""
    # k sets the slope scale in kV
    return baseline + (plateau - baseline) / (1.0 + np.exp(-(x - x0) / k))


def fit_sigmoid(hv: np.ndarray, eff: np.ndarray, err: np.ndarray | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Fit sigmoid with scipy if available; fallback to a coarse grid search."""
    hv = np.asarray(hv, dtype=float)
    eff = np.asarray(eff, dtype=float)

    # Reasonable initial guesses
    plateau0 = float(np.nanmax(eff)) if np.isfinite(np.nanmax(eff)) else 1.0
    baseline0 = float(np.nanmin(eff)) if np.isfinite(np.nanmin(eff)) else 0.0
    x00 = float(np.nanmedian(hv))
    k0 = 0.1

    p0 = np.array([plateau0, x00, k0, baseline0], dtype=float)

    # Try scipy
    try:
        from scipy.optimize import curve_fit  # type: ignore

        sigma = None
        if err is not None:
            err = np.asarray(err, dtype=float)
            # avoid zero weights
            sigma = np.where(err > 0, err, np.nanmedian(err[err > 0]) if np.any(err > 0) else 1.0)

        popt, pcov = curve_fit(
            sigmoid,
            hv,
            eff,
            p0=p0,
            sigma=sigma,
            absolute_sigma=True if sigma is not None else False,
            maxfev=20000,
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except Exception:
        # Fallback: coarse grid over x0 and k, solve plateau/baseline by least squares
        # This is less precise but robust with no scipy.
        hv_min, hv_max = float(np.min(hv)), float(np.max(hv))
        x0_grid = np.linspace(hv_min, hv_max, 80)
        k_grid = np.linspace(0.02, 0.5, 80)

        best = None
        best_params = None

        w = None
        if err is not None:
            err = np.asarray(err, dtype=float)
            w = np.where(err > 0, 1.0 / (err * err), 1.0)

        for x0 in x0_grid:
            for k in k_grid:
                g = 1.0 / (1.0 + np.exp(-(hv - x0) / k))
                # Model: baseline + (plateau-baseline)*g = baseline*(1-g) + plateau*g
                # Linear in (plateau, baseline)
                A = np.vstack([g, 1.0 - g]).T
                if w is not None:
                    Aw = A * np.sqrt(w[:, None])
                    yw = eff * np.sqrt(w)
                else:
                    Aw = A
                    yw = eff
                try:
                    coeff, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
                    plateau, baseline = float(coeff[0]), float(coeff[1])
                except Exception:
                    continue

                pred = baseline + (plateau - baseline) * g
                if w is not None:
                    sse = float(np.sum(w * (eff - pred) ** 2))
                else:
                    sse = float(np.sum((eff - pred) ** 2))

                if best is None or sse < best:
                    best = sse
                    best_params = (plateau, x0, k, baseline)

        if best_params is None:
            return p0, np.full_like(p0, np.nan)

        # No covariance estimate in fallback
        return np.array(best_params, dtype=float), np.full(4, np.nan, dtype=float)


def load_summaries(files: List[str]) -> List[Dict[str, Any]]:
    out = []
    for fn in files:
        with open(fn, "r", encoding="utf-8") as f:
            out.append(json.load(f))
    return out


def main() -> None:
    t_start = time.perf_counter()

    ap = argparse.ArgumentParser(description="Plot and fit efficiency vs HV from summary JSON files.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--glob", dest="glob_pat", help="Glob pattern for summary JSON files")
    src.add_argument("--files", nargs="+", help="Explicit list of summary JSON files")

    ap.add_argument("--use", choices=["roi", "global"], default="roi", help="Which efficiency to plot")
    ap.add_argument("--roi-name", default=None, help="If summaries contain multiple ROIs (future), select by name")
    ap.add_argument("--out", default="eff_vs_hv.png", help="Output plot PNG")
    ap.add_argument("--fit", action="store_true", help="Fit a sigmoid and overlay")
    ap.add_argument("--fit-out", default="sigmoid_fit.json", help="Write fit params JSON (if --fit)")
    ap.add_argument("--title", default=None, help="Plot title")

    args = ap.parse_args()

    # Resolve file list
    t0 = time.perf_counter()
    if args.glob_pat:
        files = sorted(glob.glob(args.glob_pat))
    else:
        files = args.files
    if not files:
        raise SystemExit("No summary files found.")

    summaries = load_summaries(files)
    t1 = time.perf_counter()

    # Extract arrays
    hv_list: List[float] = []
    eff_list: List[float] = []
    err_list: List[float] = []
    tag_list: List[str] = []

    for s in summaries:
        hv = float(s["hv_kv"])
        tag = str(s.get("tag", ""))
        if args.use == "roi":
            eff = float(s["roi_efficiency"]["eff"])
            err = float(s["roi_efficiency"]["err_binom"])
        else:
            eff = float(s["global_efficiency"]["eff"])
            err = float(s["global_efficiency"]["err_binom"])

        hv_list.append(hv)
        eff_list.append(eff)
        err_list.append(err)
        tag_list.append(tag)

    # Sort by HV
    order = np.argsort(hv_list)
    hv = np.asarray([hv_list[i] for i in order], dtype=float)
    eff = np.asarray([eff_list[i] for i in order], dtype=float)
    err = np.asarray([err_list[i] for i in order], dtype=float)

    # Plot
    t2 = time.perf_counter()
    plt.figure()
    plt.errorbar(hv, eff, yerr=err, fmt="o", capsize=3)
    plt.xlabel("HV [kV]")
    plt.ylabel("Efficiency")

    if args.title:
        plt.title(args.title)
    else:
        plt.title(f"{args.use.upper()} efficiency vs HV")

    fit_payload = None
    if args.fit and len(hv) >= 4:
        popt, perr = fit_sigmoid(hv, eff, err)
        hv_grid = np.linspace(float(np.min(hv)), float(np.max(hv)), 400)
        plt.plot(hv_grid, sigmoid(hv_grid, *popt))

        fit_payload = {
            "model": "4-parameter logistic",
            "params": {
                "plateau": float(popt[0]),
                "x0": float(popt[1]),
                "k": float(popt[2]),
                "baseline": float(popt[3]),
            },
            "params_err": {
                "plateau": float(perr[0]) if np.isfinite(perr[0]) else None,
                "x0": float(perr[1]) if np.isfinite(perr[1]) else None,
                "k": float(perr[2]) if np.isfinite(perr[2]) else None,
                "baseline": float(perr[3]) if np.isfinite(perr[3]) else None,
            },
            "n_points": int(len(hv)),
            "source_files": [os.path.abspath(x) for x in files],
        }

    plt.ylim(-0.05, 1.05)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    t3 = time.perf_counter()

    if fit_payload is not None:
        with open(args.fit_out, "w", encoding="utf-8") as f:
            json.dump(fit_payload, f, indent=2)
        print(f"Wrote fit params: {args.fit_out}")

    print(f"Wrote plot: {args.out}")
    print(f"Loaded {len(files)} summary files")
    print(f"Timing: load_json={(t1 - t0):.3f}s, plot+fit={(t3 - t2):.3f}s, total={(time.perf_counter() - t_start):.3f}s")


if __name__ == "__main__":
    main()
