#!/usr/bin/env python3
"""
rpc_plot_eff_vs_hv_cms.py

Create a CMS/RPC-group style "Efficiency vs HV" plot from summary_HV*.json files
produced by rpc_make_summary.py.

Expected summary JSON schema (as produced by rpc_make_summary.py):
  - hv_kv
  - roi_efficiency:    { eff, err_binom }
  - global_efficiency: { eff, err_binom }

The script can plot ROI or global efficiency and reports:
  - plateau[%] (from sigmoid fit eps_max)
  - WP: HV@95%plateau + offset (default 100 V)
  - Eff(WP)[%]

Notes:
  - Uses non-interactive backend to avoid Wayland/Qt GUI warnings.
  - Uses scipy for robust fit if available; otherwise uses a reasonable fallback.
"""
import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HV_MAP = {
    "HV1":  6.0, "HV2":  6.5, "HV3":  6.7, "HV4":  6.8, "HV5":  6.9,
    "HV6":  7.0, "HV7":  7.1, "HV8":  7.2, "HV9":  7.3, "HV10": 7.4, "HV11": 7.5,
}

def sigmoid(hv, eps0, eps_max, hv50, k):
    return eps0 + (eps_max - eps0) / (1.0 + np.exp(-(hv - hv50) / k))

def fit_sigmoid(hv, eff):
    hv = np.asarray(hv, float)
    eff = np.asarray(eff, float)

    eps0_guess = max(0.0, float(np.min(eff)))
    eps_max_guess = min(1.0, float(np.max(eff)))
    hv50_guess = float(hv[np.argmin(np.abs(eff - 0.5 * (eps0_guess + eps_max_guess)))])
    k_guess = 0.05  # ~50 V

    p0 = [eps0_guess, eps_max_guess, hv50_guess, k_guess]

    try:
        from scipy.optimize import curve_fit
        bounds = ([0.0, 0.0, hv.min() - 1.0, 1e-3],
                  [1.0, 1.0, hv.max() + 1.0, 1.0])
        popt, _ = curve_fit(sigmoid, hv, eff, p0=p0, bounds=bounds, maxfev=20000)
        return tuple(map(float, popt))
    except Exception:
        return tuple(map(float, p0))

def hv_at_fraction_of_plateau(eps0, eps_max, hv50, k, frac):
    frac = min(max(float(frac), 1e-6), 1 - 1e-6)
    return hv50 + k * np.log(frac / (1.0 - frac))

def infer_hv_from_filename(p: Path) -> float:
    m = re.search(r"(HV\d+)", p.name)
    if m and m.group(1) in HV_MAP:
        return HV_MAP[m.group(1)]
    m = re.search(r"(HV\d+)", str(p.parent))
    if m and m.group(1) in HV_MAP:
        return HV_MAP[m.group(1)]
    raise RuntimeError(f"Cannot infer HV from filename: {p}")

def load_series(glob_pattern: str, use_roi: bool):
    paths = sorted(Path().glob(glob_pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {glob_pattern}")

    hv, eff, err = [], [], []
    for p in paths:
        d = json.loads(p.read_text())
        hv_kv = float(d.get("hv_kv", infer_hv_from_filename(p)))

        if use_roi:
            e = float(d["roi_efficiency"]["eff"])
            de = float(d["roi_efficiency"]["err_binom"])
        else:
            e = float(d["global_efficiency"]["eff"])
            de = float(d["global_efficiency"]["err_binom"])

        hv.append(hv_kv); eff.append(e); err.append(de)

    hv = np.array(hv, float); eff = np.array(eff, float); err = np.array(err, float)
    order = np.argsort(hv)
    return hv[order], eff[order], err[order]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="+", required=True,
                    help="label::globpattern, e.g. 'RxLW = 12 BX::out_hv/summaries/summary_HV*.json'")
    ap.add_argument("--out", default="efficiency_cms_style.png")
    ap.add_argument("--use-roi", action="store_true", help="Plot ROI efficiency (default: global).")

    ap.add_argument("--cms-left", default="CMS Preliminary")
    ap.add_argument("--cms-right", default="CERN 904 Lab")
    ap.add_argument("--textbox", default="chamber 236\n1.4mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40fC (dac10)")

    ap.add_argument("--ymax", type=float, default=110.0)
    ap.add_argument("--wp-frac", type=float, default=0.95)
    ap.add_argument("--wp-offset-v", type=float, default=100.0)

    args = ap.parse_args()

    fig = plt.figure(figsize=(7.2, 7.2))
    ax = plt.gca()

    ax.set_xlabel("HV [kV]", fontsize=18)
    ax.set_ylabel("Efficiency [%]", fontsize=22)
    ax.set_xlim(6.0, 7.5)
    ax.set_ylim(0.0, args.ymax)
    ax.grid(True, which="both", alpha=0.25)
    ax.minorticks_on()
    ax.axhline(100.0, linestyle="--", linewidth=1.2, alpha=0.6)

    ax.text(0.02, 1.03, args.cms_left, transform=ax.transAxes,
            fontsize=22, fontweight="bold", va="bottom")
    ax.text(0.98, 1.03, args.cms_right, transform=ax.transAxes,
            fontsize=18, ha="right", va="bottom")

    ax.text(0.10, 0.55, args.textbox, transform=ax.transAxes,
            fontsize=14, va="top", ha="left")

    for spec in args.series:
        if "::" not in spec:
            raise SystemExit(f"Bad --series '{spec}'. Use label::globpattern")
        label, globpat = spec.split("::", 1)

        hv, eff, err = load_series(globpat, use_roi=args.use_roi)
        eps0, eps_max, hv50, k = fit_sigmoid(hv, eff)

        plateau_pct = 100.0 * eps_max
        hv95 = hv_at_fraction_of_plateau(eps0, eps_max, hv50, k, args.wp_frac)
        wp_kv = hv95 + args.wp_offset_v / 1000.0
        wp_v = wp_kv * 1000.0
        eff_wp_pct = 100.0 * sigmoid(wp_kv, eps0, eps_max, hv50, k)

        ax.errorbar(
            hv, 100.0 * eff, yerr=100.0 * err, fmt="x", capsize=3,
            label=f"{label}, plateau[%] = {plateau_pct:.2f}, WP = {wp_v:.2f} V, Eff(WP)[%] = {eff_wp_pct:.2f}"
        )

        hv_fine = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 400)
        ax.plot(hv_fine, 100.0 * sigmoid(hv_fine, eps0, eps_max, hv50, k), linewidth=1.5, alpha=0.85)

    ax.legend(loc="upper left", fontsize=10, frameon=False)
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Saved: {args.out}")

if __name__ == "__main__":
    main()
