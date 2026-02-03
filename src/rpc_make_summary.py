#!/usr/bin/env python3
"""make_eff_summary.py

Stage 2a (summarizer) in a clean pipeline:

  Stage 1: tracklets_reconstruction_efficiency_v2.py -> eff_map_HV*.root
  Stage 2a: make_eff_summary.py                     -> summary_HV*.json
  Stage 2b: efficiency_vs_hv.py                     -> plots + sigmoid fit

This script reads a ROOT file produced by Stage 1, integrates the numerator
and denominator TH2 histograms in a chosen ROI, and writes a small JSON summary.

Why JSON summaries?
- HV scan plotting & fitting becomes instant (no event loops).
- You can change ROI / fits / plots without rerunning reconstruction.

The efficiency uncertainty is computed with a simple binomial approximation:
  err = sqrt(e*(1-e)/N) for N>0.

Usage examples
  python3 make_eff_summary.py --root eff_map_HV6p8.root --hv 6.8 \
      --roi-name scint_overlap --roi-rect 0 22 0 130

  python3 make_eff_summary.py --root eff_map_HV7p0.root --hv 7.0 \
      --roi-name full_active --roi-rect 0 60 0 130 \
      --out summary_HV7p0.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Any, Dict, Tuple

import ROOT


def _binomial_err(matched: float, eligible: float) -> float:
    if eligible <= 0:
        return 0.0
    e = matched / eligible
    # Guard numerical issues
    if e < 0:
        e = 0.0
    if e > 1:
        e = 1.0
    return math.sqrt(e * (1.0 - e) / eligible)


def _get_th2(file: ROOT.TFile, name: str) -> ROOT.TH2:
    h = file.Get(name)
    if not h:
        raise RuntimeError(f"Histogram '{name}' not found in ROOT file")
    # Detach from file
    h.SetDirectory(0)
    return h


def _roi_bin_range(h: ROOT.TH2, x_min: float, x_max: float, y_min: float, y_max: float) -> Tuple[int, int, int, int]:
    ax = h.GetXaxis()
    ay = h.GetYaxis()

    # FindBin uses inclusive behavior; clamp to axis range.
    bx1 = ax.FindBin(x_min)
    bx2 = ax.FindBin(x_max)
    by1 = ay.FindBin(y_min)
    by2 = ay.FindBin(y_max)

    bx1 = max(1, min(bx1, ax.GetNbins()))
    bx2 = max(1, min(bx2, ax.GetNbins()))
    by1 = max(1, min(by1, ay.GetNbins()))
    by2 = max(1, min(by2, ay.GetNbins()))

    if bx2 < bx1:
        bx1, bx2 = bx2, bx1
    if by2 < by1:
        by1, by2 = by2, by1

    return bx1, bx2, by1, by2


def _integral_2d(h: ROOT.TH2, bx1: int, bx2: int, by1: int, by2: int) -> float:
    return float(h.Integral(bx1, bx2, by1, by2))


def main() -> None:
    t_start = time.perf_counter()

    parser = argparse.ArgumentParser(description="Summarize per-HV efficiency maps into a compact JSON.")
    parser.add_argument("--root", required=True, help="Input ROOT file from v2 (contains TH2 histograms)")
    parser.add_argument("--hv", type=float, required=True, help="Working point high voltage (kV)")
    parser.add_argument("--tag", default="", help="Optional run tag (string)")
    parser.add_argument("--h-total", default="h_total", help="Denominator TH2 name (default: h_total)")
    parser.add_argument("--h-tracklet", default="h_tracklet", help="Numerator TH2 name (default: h_tracklet)")

    parser.add_argument("--roi-name", default="roi", help="ROI name stored in JSON")
    parser.add_argument(
        "--roi-rect",
        nargs=4,
        type=float,
        metavar=("X_MIN", "X_MAX", "Y_MIN", "Y_MAX"),
        default=None,
        help="Rectangle ROI in detector coordinates. If omitted, ROI=full histogram range.",
    )

    parser.add_argument("--out", default=None, help="Output JSON filename (default derived from input)")

    args = parser.parse_args()

    # I/O timing
    t0 = time.perf_counter()
    f = ROOT.TFile.Open(args.root)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {args.root}")

    h_total = _get_th2(f, args.h_total)
    h_num = _get_th2(f, args.h_tracklet)
    f.Close()
    t1 = time.perf_counter()

    # Determine ROI
    if args.roi_rect is None:
        ax = h_total.GetXaxis()
        ay = h_total.GetYaxis()
        x_min, x_max = ax.GetXmin(), ax.GetXmax()
        y_min, y_max = ay.GetXmin(), ay.GetXmax()
    else:
        x_min, x_max, y_min, y_max = args.roi_rect

    bx1, bx2, by1, by2 = _roi_bin_range(h_total, x_min, x_max, y_min, y_max)

    # Compute integrals
    t2 = time.perf_counter()
    eligible_roi = _integral_2d(h_total, bx1, bx2, by1, by2)
    matched_roi = _integral_2d(h_num, bx1, bx2, by1, by2)

    eligible_all = float(h_total.Integral())
    matched_all = float(h_num.Integral())

    eff_roi = (matched_roi / eligible_roi) if eligible_roi > 0 else 0.0
    err_roi = _binomial_err(matched_roi, eligible_roi)

    eff_all = (matched_all / eligible_all) if eligible_all > 0 else 0.0
    err_all = _binomial_err(matched_all, eligible_all)
    t3 = time.perf_counter()

    payload: Dict[str, Any] = {
        "hv_kv": float(args.hv),
        "tag": args.tag,
        "source_root": os.path.abspath(args.root),
        "h_total": args.h_total,
        "h_tracklet": args.h_tracklet,
        "roi": {
            "name": args.roi_name,
            "type": "rect",
            "x_min": float(x_min),
            "x_max": float(x_max),
            "y_min": float(y_min),
            "y_max": float(y_max),
            "bin_x1": int(bx1),
            "bin_x2": int(bx2),
            "bin_y1": int(by1),
            "bin_y2": int(by2),
        },
        "roi_counts": {
            "eligible": float(eligible_roi),
            "matched": float(matched_roi),
        },
        "roi_efficiency": {
            "eff": float(eff_roi),
            "err_binom": float(err_roi),
        },
        "global_counts": {
            "eligible": float(eligible_all),
            "matched": float(matched_all),
        },
        "global_efficiency": {
            "eff": float(eff_all),
            "err_binom": float(err_all),
        },
        "timing_sec": {
            "read_root": float(t1 - t0),
            "compute": float(t3 - t2),
            "total": float(time.perf_counter() - t_start),
        },
    }

    out_path = args.out
    if out_path is None:
        base = os.path.splitext(os.path.basename(args.root))[0]
        # Make filename stable across locales
        hv_str = str(args.hv).replace(".", "p")
        out_path = f"summary_{base}_HV{hv_str}.json"

    with open(out_path, "w", encoding="utf-8") as out:
        json.dump(payload, out, indent=2)

    print(f"Wrote JSON summary: {out_path}")
    print(
        f"ROI eff = {eff_roi:.4f} ± {err_roi:.4f}  (matched={matched_roi:.0f}, eligible={eligible_roi:.0f})"
    )
    print(
        f"Global eff = {eff_all:.4f} ± {err_all:.4f} (matched={matched_all:.0f}, eligible={eligible_all:.0f})"
    )
    print(
        f"Timing: read_root={payload['timing_sec']['read_root']:.3f}s, compute={payload['timing_sec']['compute']:.3f}s, total={payload['timing_sec']['total']:.3f}s"
    )


if __name__ == "__main__":
    main()
