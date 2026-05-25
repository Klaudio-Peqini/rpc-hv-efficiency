#!/usr/bin/env python3
"""
rpc_tracklet_efficiency.py

Option C-final N-chamber producer with diagnostic histograms.

Core outputs, kept compatible with the existing downstream workflow:
  - h_total     : denominator / eligible events in anchor coordinates
  - h_tracklet  : numerator / matched reconstructed tracklets in anchor coordinates
  - h_eff       : local efficiency map, h_tracklet / h_total

Additional diagnostic outputs written to the same ROOT file:
  Candidate-level residuals, before the matching decision:
    - h_dx_raw_candidate, h_dy_raw_candidate, h_dx_raw_dy_raw_candidate
    - h_dx_candidate,     h_dy_candidate,     h_dx_dy_candidate
    - h_dt_candidate

  Matched-tracklet diagnostics, after the matching decision:
    - h_dx, h_dy, h_dt, h_dx_dy
    - h_theta_x, h_theta_y, h_theta, h_theta_x_theta_y

Angle convention:
  - For >=3 chamber runs, theta_x/theta_y are computed from the fitted
    reference-track slopes x(z)=a_x z+b_x and y(z)=a_y z+b_y.
  - For 2 chamber runs, theta_x/theta_y are computed from the corrected
    residuals divided by dz(anchor,target). These are useful diagnostic
    pair angles, not independent reference-track angles.

The optional CSV/NDJSON outputs now also contain theta_x_deg, theta_y_deg,
and theta_deg, which makes correlation-matrix plotting straightforward.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import ROOT  # type: ignore


def stamp(t0: float, msg: str) -> None:
    print(f"[{time.time() - t0:8.2f} s] {msg}")


def pick_hit(ev, hit_mode: str = "earliest") -> Optional[Tuple[float, float, float]]:
    """Pick one representative hit from vector branches."""
    xs = np.asarray(ev.Strip_X, dtype=float)
    ys = np.asarray(ev.Strip_Y, dtype=float)
    ts = np.asarray(ev.Strip_Tdc_diff, dtype=float)

    if xs.size == 0 or ys.size == 0 or ts.size == 0:
        return None

    n = min(xs.size, ys.size, ts.size)
    xs, ys, ts = xs[:n], ys[:n], ts[:n]

    if hit_mode == "earliest":
        i = int(np.argmin(ts))
        return float(xs[i]), float(ys[i]), float(ts[i])

    if hit_mode == "centroid":
        return float(xs.mean()), float(ys.mean()), float(ts.mean())

    raise ValueError(f"Unknown hit_mode={hit_mode}")


def fit_line(zs: Sequence[float], vs: Sequence[float]) -> Tuple[float, float]:
    """Least-squares fit v(z)=a*z+b for >=2 points."""
    z = np.asarray(zs, dtype=float)
    v = np.asarray(vs, dtype=float)
    A = np.vstack([z, np.ones_like(z)]).T
    (a, b), *_ = np.linalg.lstsq(A, v, rcond=None)
    return float(a), float(b)


def compute_eff_hist(h_num, h_den, name: str = "h_eff"):
    h_eff = h_num.Clone(name)
    h_eff.Reset()
    h_eff.Divide(h_num, h_den, 1.0, 1.0, "B")
    return h_eff


def ttest_mean_zero(samples: Sequence[float]) -> Tuple[float, float, int]:
    """
    One-sample t-test for H0: mean=0, two-sided.

    Uses SciPy if available; otherwise falls back to a normal approximation,
    which is adequate for the large residual samples normally used here.
    Returns (mean, p_value, n).
    """
    x = np.asarray(samples, dtype=float)
    n = int(x.size)
    if n < 2:
        return (float(x.mean()) if n else 0.0, 1.0, n)

    mean = float(x.mean())
    s = float(x.std(ddof=1))
    if s == 0.0:
        return mean, (0.0 if mean != 0.0 else 1.0), n

    t = mean / (s / math.sqrt(n))
    df = n - 1

    try:
        from scipy.stats import t as student_t  # type: ignore

        p = 2.0 * float(student_t.sf(abs(t), df))
        return mean, p, n
    except Exception:
        p = 2.0 * 0.5 * math.erfc(abs(t) / math.sqrt(2.0))
        return mean, p, n


def make_diagnostic_histograms() -> Dict[str, object]:
    """Create residual and angular diagnostic histograms."""
    hists = {
        # Raw candidate residuals: before dx0/dy0 alignment correction.
        "h_dx_raw_candidate": ROOT.TH1D(
            "h_dx_raw_candidate",
            "Candidate raw residual #Delta x;#Delta x_{raw} [cm];Candidates",
            160,
            -40.0,
            40.0,
        ),
        "h_dy_raw_candidate": ROOT.TH1D(
            "h_dy_raw_candidate",
            "Candidate raw residual #Delta y;#Delta y_{raw} [cm];Candidates",
            160,
            -80.0,
            80.0,
        ),
        "h_dx_raw_dy_raw_candidate": ROOT.TH2D(
            "h_dx_raw_dy_raw_candidate",
            "Candidate raw #Delta x vs raw #Delta y;#Delta x_{raw} [cm];#Delta y_{raw} [cm]",
            120,
            -40.0,
            40.0,
            120,
            -80.0,
            80.0,
        ),
        # Corrected candidate residuals: after dx0/dy0 alignment correction.
        "h_dx_candidate": ROOT.TH1D(
            "h_dx_candidate",
            "Candidate corrected residual #Delta x;#Delta x [cm];Candidates",
            160,
            -40.0,
            40.0,
        ),
        "h_dy_candidate": ROOT.TH1D(
            "h_dy_candidate",
            "Candidate corrected residual #Delta y;#Delta y [cm];Candidates",
            160,
            -80.0,
            80.0,
        ),
        "h_dt_candidate": ROOT.TH1D(
            "h_dt_candidate",
            "Candidate timing residual #Delta t;#Delta t [TDC units];Candidates",
            160,
            -20.0,
            20.0,
        ),
        "h_dx_dy_candidate": ROOT.TH2D(
            "h_dx_dy_candidate",
            "Candidate corrected #Delta x vs #Delta y;#Delta x [cm];#Delta y [cm]",
            120,
            -40.0,
            40.0,
            120,
            -80.0,
            80.0,
        ),
        # Matched residuals.
        "h_dx": ROOT.TH1D(
            "h_dx",
            "Matched residual #Delta x;#Delta x [cm];Matched tracklets",
            160,
            -40.0,
            40.0,
        ),
        "h_dy": ROOT.TH1D(
            "h_dy",
            "Matched residual #Delta y;#Delta y [cm];Matched tracklets",
            160,
            -80.0,
            80.0,
        ),
        "h_dt": ROOT.TH1D(
            "h_dt",
            "Matched timing residual #Delta t;#Delta t [TDC units];Matched tracklets",
            160,
            -20.0,
            20.0,
        ),
        "h_dx_dy": ROOT.TH2D(
            "h_dx_dy",
            "Matched #Delta x vs #Delta y;#Delta x [cm];#Delta y [cm]",
            120,
            -40.0,
            40.0,
            120,
            -80.0,
            80.0,
        ),
        # Matched angular diagnostics.
        "h_theta_x": ROOT.TH1D(
            "h_theta_x",
            "Tracklet angle #theta_{x};#theta_{x} [deg];Matched tracklets",
            160,
            -45.0,
            45.0,
        ),
        "h_theta_y": ROOT.TH1D(
            "h_theta_y",
            "Tracklet angle #theta_{y};#theta_{y} [deg];Matched tracklets",
            160,
            -70.0,
            70.0,
        ),
        "h_theta": ROOT.TH1D(
            "h_theta",
            "Total tracklet inclination #theta;#theta [deg];Matched tracklets",
            160,
            0.0,
            80.0,
        ),
        "h_theta_x_theta_y": ROOT.TH2D(
            "h_theta_x_theta_y",
            "#theta_{x} vs #theta_{y};#theta_{x} [deg];#theta_{y} [deg]",
            120,
            -45.0,
            45.0,
            120,
            -70.0,
            70.0,
        ),
    }

    for hist in hists.values():
        hist.Sumw2()

    return hists


def finite_angle(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(float(value))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="RPC N-ch tracklet efficiency producer with dynamic-cut support"
    )
    ap.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="ROOT files, one per chamber, same event order for one HV point.",
    )
    ap.add_argument(
        "--z",
        nargs="+",
        required=True,
        type=float,
        help="z positions in the same order as --files.",
    )
    ap.add_argument("--tree", default="events")
    ap.add_argument("--tag", required=True, help="Output tag, e.g. HV6")
    ap.add_argument(
        "--target-index",
        type=int,
        default=-1,
        help="Target chamber index. Default: last chamber.",
    )
    ap.add_argument(
        "--anchor-index",
        type=int,
        default=0,
        help="Anchor chamber defining denominator map and reference time.",
    )

    # v2-like selection knobs / Option C semantics.
    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])
    ap.add_argument(
        "--regime",
        default="normal",
        choices=["normal", "beam"],
        help="Running regime label propagated by the HV pipeline.",
    )
    ap.add_argument(
        "--v2-eligibility",
        action="store_true",
        help="Kept for compatibility. Option C already follows v2-like eligibility semantics.",
    )
    ap.add_argument(
        "--require-ncluster-eq1",
        action="store_true",
        help="Require nCluster==1 for all reference chambers, excluding target.",
    )
    ap.add_argument(
        "--reject-target-ncluster-gt1",
        action="store_true",
        help="Kept for compatibility. Target nCluster>1 is rejected before eligible++.",
    )
    ap.add_argument(
        "--cluster-size-max",
        type=int,
        default=6,
        help="Cluster-size veto threshold using len(Strip_Tdc_diff).",
    )
    ap.add_argument(
        "--cluster-size-on-all-ref",
        action="store_true",
        help="Apply cluster-size veto to all reference chambers. Default: anchor only.",
    )
    ap.add_argument(
        "--theta-max-deg",
        type=float,
        default=None,
        help="Optional maximum reference-track inclination angle in degrees.",
    )

    # Matching tolerances.
    ap.add_argument("--tol-x", type=float, default=7.0)
    ap.add_argument("--tol-y", type=float, default=20.0)
    ap.add_argument("--tol-t", type=float, default=3.0)

    # Shift / hypothesis testing.
    ap.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="Significance level for rejecting H0: mean(Delta)=0.",
    )
    ap.add_argument("--disable-shift", action="store_true", help="Disable automatic dx/dy shifting.")
    ap.add_argument(
        "--min-residuals",
        type=int,
        default=50,
        help="Minimum number of residual samples needed for the hypothesis test.",
    )

    # Histogram configuration.
    ap.add_argument("--nbins-x", type=int, default=64)
    ap.add_argument("--nbins-y", type=int, default=128)
    ap.add_argument("--x-range", nargs=2, type=float, default=[0.0, 60.0])
    ap.add_argument("--y-range", nargs=2, type=float, default=[0.0, 140.0])

    # Outputs.
    ap.add_argument("--write-csv", action="store_true")
    ap.add_argument("--write-ndjson", action="store_true")
    ap.add_argument(
        "--no-diagnostic-hists",
        action="store_true",
        help="Do not write residual/angular diagnostic histograms into the ROOT output.",
    )

    args = ap.parse_args()

    t0 = time.time()
    stamp(t0, f"Starting Option C-final N-ch efficiency producer (regime={args.regime})")

    files = [Path(p) for p in args.files]
    zs = list(args.z)

    if len(files) != len(zs):
        raise SystemExit("ERROR: --files and --z must have the same length.")

    n_ch = len(files)
    if n_ch < 2:
        raise SystemExit("ERROR: Need at least 2 chambers.")

    tgt = args.target_index if args.target_index >= 0 else n_ch + args.target_index
    anc = args.anchor_index if args.anchor_index >= 0 else n_ch + args.anchor_index

    if not (0 <= tgt < n_ch) or not (0 <= anc < n_ch):
        raise SystemExit("ERROR: invalid target-index or anchor-index.")
    if tgt == anc:
        raise SystemExit("ERROR: target-index cannot equal anchor-index.")

    ref_idx = [i for i in range(n_ch) if i != tgt]
    if anc not in ref_idx:
        raise SystemExit("ERROR: internal: anchor must be in reference set.")

    # Open files / trees.
    roots = []
    trees = []
    for f in files:
        rf = ROOT.TFile.Open(str(f))
        if not rf or rf.IsZombie():
            raise SystemExit(f"ERROR: cannot open {f}")
        tr = rf.Get(args.tree)
        if not tr:
            available = [k.GetName() for k in rf.GetListOfKeys()]
            raise SystemExit(
                f"ERROR: tree '{args.tree}' not found in {f}. Available objects: {available}"
            )
        roots.append(rf)
        trees.append(tr)

    n_entries = min(int(t.GetEntries()) for t in trees)
    stamp(t0, f"Loaded {n_ch} trees, using first {n_entries} entries from tree '{args.tree}'")

    # Core histograms in anchor coordinates.
    h_den = ROOT.TH2D(
        "h_total",
        "Denominator;x_{anchor} [cm];y_{anchor} [cm]",
        args.nbins_x,
        args.x_range[0],
        args.x_range[1],
        args.nbins_y,
        args.y_range[0],
        args.y_range[1],
    )
    h_num = ROOT.TH2D(
        "h_tracklet",
        "Numerator;x_{anchor} [cm];y_{anchor} [cm]",
        args.nbins_x,
        args.x_range[0],
        args.x_range[1],
        args.nbins_y,
        args.y_range[0],
        args.y_range[1],
    )
    h_den.Sumw2()
    h_num.Sumw2()

    diagnostic_hists = {} if args.no_diagnostic_hists else make_diagnostic_histograms()

    # Record outputs only on the final pass.
    csv_rows: List[Tuple[object, ...]] = []
    ndjson_path = Path("tracklets.ndjson")
    if args.write_ndjson:
        ndjson_path.write_text("", encoding="utf-8")

    def compute_angles(
        *,
        ax_slope: Optional[float],
        ay_slope: Optional[float],
        dx: float,
        dy: float,
    ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """Return theta_x_deg, theta_y_deg, theta_total_deg."""
        if ax_slope is not None and ay_slope is not None:
            theta_x_deg = math.degrees(math.atan(ax_slope))
            theta_y_deg = math.degrees(math.atan(ay_slope))
            theta_total_deg = math.degrees(math.atan2(math.sqrt(ax_slope * ax_slope + ay_slope * ay_slope), 1.0))
            return theta_x_deg, theta_y_deg, theta_total_deg

        dz_anchor_target = float(zs[tgt] - zs[anc])
        if abs(dz_anchor_target) <= 0.0:
            return None, None, None

        theta_x_deg = math.degrees(math.atan2(dx, dz_anchor_target))
        theta_y_deg = math.degrees(math.atan2(dy, dz_anchor_target))
        theta_total_deg = math.degrees(
            math.atan2(math.sqrt(dx * dx + dy * dy), abs(dz_anchor_target))
        )
        return theta_x_deg, theta_y_deg, theta_total_deg

    def loop(
        *,
        apply_shift: bool,
        dx0: float,
        dy0: float,
        fill_hists: bool,
        collect_residuals: bool,
    ) -> Tuple[int, int, int, List[float], List[float]]:
        eligible = 0
        matched = 0
        ghosts = 0
        dx_list: List[float] = []
        dy_list: List[float] = []

        for i in range(n_entries):
            for t in trees:
                t.GetEntry(i)

            evA = trees[anc]
            evT = trees[tgt]

            # v2-like eligibility preselection.
            # Anchor must be clean. Target nCluster>1 is rejected before eligible++.
            if int(evA.nCluster) != 1:
                continue

            if int(evT.nCluster) > 1:
                continue

            # Anchor cluster-size veto.
            if len(evA.Strip_Tdc_diff) > args.cluster_size_max:
                continue

            # Additional reference chambers used for prediction.
            hits_ref: Dict[int, Tuple[float, float, float]] = {}
            ok = True
            for j in ref_idx:
                ev = trees[j]

                if j != anc:
                    if args.require_ncluster_eq1 and int(ev.nCluster) != 1:
                        ok = False
                        break
                    if args.cluster_size_on_all_ref and len(ev.Strip_Tdc_diff) > args.cluster_size_max:
                        ok = False
                        break

                hit = pick_hit(ev, args.hit_mode)
                if hit is None:
                    ok = False
                    break
                hits_ref[j] = hit

            if not ok:
                continue

            # Eligible event after clean reference preselection.
            eligible += 1

            hitA = hits_ref.get(anc)
            if hitA is None:
                continue
            xA, yA, tA = hitA

            if fill_hists:
                h_den.Fill(xA, yA)

            # Ghost check after denominator filling.
            if int(evT.nCluster) == 0:
                ghosts += 1
                continue

            # Predict target position from references.
            z_ref: List[float] = []
            x_ref: List[float] = []
            y_ref: List[float] = []
            for j in ref_idx:
                z_ref.append(float(zs[j]))
                x_ref.append(hits_ref[j][0])
                y_ref.append(hits_ref[j][1])

            z_t = float(zs[tgt])
            ax_slope: Optional[float] = None
            ay_slope: Optional[float] = None

            if len(z_ref) >= 2:
                ax_slope, bx = fit_line(z_ref, x_ref)
                ay_slope, by = fit_line(z_ref, y_ref)
                x_pred = ax_slope * z_t + bx
                y_pred = ay_slope * z_t + by
                theta_ref_deg = math.degrees(math.atan2(math.sqrt(ax_slope * ax_slope + ay_slope * ay_slope), 1.0))
                if args.theta_max_deg is not None and theta_ref_deg > args.theta_max_deg:
                    continue
            else:
                x_pred, y_pred = xA, yA

            t_pred = tA

            hitT = pick_hit(evT, args.hit_mode)
            if hitT is None:
                ghosts += 1
                continue
            xT, yT, tT = hitT

            dx_raw = xT - x_pred
            dy_raw = yT - y_pred
            dt = tT - t_pred

            dx = dx_raw - dx0 if apply_shift else dx_raw
            dy = dy_raw - dy0 if apply_shift else dy_raw

            theta_x_deg, theta_y_deg, theta_deg = compute_angles(
                ax_slope=ax_slope,
                ay_slope=ay_slope,
                dx=dx,
                dy=dy,
            )

            # Residual samples for alignment hypothesis test are collected from
            # candidate pairs before applying the final matching window.
            if collect_residuals:
                dx_list.append(dx_raw)
                dy_list.append(dy_raw)

            if fill_hists and diagnostic_hists:
                diagnostic_hists["h_dx_raw_candidate"].Fill(dx_raw)
                diagnostic_hists["h_dy_raw_candidate"].Fill(dy_raw)
                diagnostic_hists["h_dx_raw_dy_raw_candidate"].Fill(dx_raw, dy_raw)
                diagnostic_hists["h_dx_candidate"].Fill(dx)
                diagnostic_hists["h_dy_candidate"].Fill(dy)
                diagnostic_hists["h_dt_candidate"].Fill(dt)
                diagnostic_hists["h_dx_dy_candidate"].Fill(dx, dy)

            is_match = (abs(dt) < args.tol_t) and (abs(dx) < args.tol_x) and (abs(dy) < args.tol_y)
            if not is_match:
                continue

            matched += 1

            if fill_hists:
                h_num.Fill(xA, yA)

                if diagnostic_hists:
                    diagnostic_hists["h_dx"].Fill(dx)
                    diagnostic_hists["h_dy"].Fill(dy)
                    diagnostic_hists["h_dt"].Fill(dt)
                    diagnostic_hists["h_dx_dy"].Fill(dx, dy)

                    if finite_angle(theta_x_deg):
                        diagnostic_hists["h_theta_x"].Fill(float(theta_x_deg))
                    if finite_angle(theta_y_deg):
                        diagnostic_hists["h_theta_y"].Fill(float(theta_y_deg))
                    if finite_angle(theta_deg):
                        diagnostic_hists["h_theta"].Fill(float(theta_deg))
                    if finite_angle(theta_x_deg) and finite_angle(theta_y_deg):
                        diagnostic_hists["h_theta_x_theta_y"].Fill(float(theta_x_deg), float(theta_y_deg))

                if args.write_csv:
                    csv_rows.append(
                        (
                            i,
                            xA,
                            yA,
                            x_pred,
                            y_pred,
                            xT,
                            yT,
                            dx,
                            dy,
                            dt,
                            dx_raw,
                            dy_raw,
                            dx0,
                            dy0,
                            theta_x_deg,
                            theta_y_deg,
                            theta_deg,
                        )
                    )

                if args.write_ndjson:
                    rec = {
                        "entry": i,
                        "anchor_index": anc,
                        "target_index": tgt,
                        "n_chambers": n_ch,
                        "regime": args.regime,
                        "xA": xA,
                        "yA": yA,
                        "tA": tA,
                        "x_pred": x_pred,
                        "y_pred": y_pred,
                        "t_pred": t_pred,
                        "xT": xT,
                        "yT": yT,
                        "tT": tT,
                        "dx": dx,
                        "dy": dy,
                        "dt": dt,
                        "dx_raw": dx_raw,
                        "dy_raw": dy_raw,
                        "shift_applied": bool(apply_shift),
                        "dx0": dx0,
                        "dy0": dy0,
                        "theta_x_deg": theta_x_deg,
                        "theta_y_deg": theta_y_deg,
                        "theta_deg": theta_deg,
                    }
                    with ndjson_path.open("a", encoding="utf-8") as f:
                        f.write(json.dumps(rec) + "\n")

        return eligible, matched, ghosts, dx_list, dy_list

    # PASS 1: collect raw residuals for the misalignment hypothesis test.
    dx0 = 0.0
    dy0 = 0.0
    shift_applied_final = False

    if args.disable_shift:
        print("Null hypothesis not rejected and shifting not carried over")
    else:
        _, _, _, dx_list, dy_list = loop(
            apply_shift=False,
            dx0=0.0,
            dy0=0.0,
            fill_hists=False,
            collect_residuals=True,
        )

        if (len(dx_list) < args.min_residuals) or (len(dy_list) < args.min_residuals):
            print("Null hypothesis not rejected and shifting not carried over")
        else:
            mean_dx, p_dx, n_dx = ttest_mean_zero(dx_list)
            mean_dy, p_dy, n_dy = ttest_mean_zero(dy_list)
            reject = (p_dx < args.alpha) or (p_dy < args.alpha)

            print("\n--- Hypothesis tests for misalignment (H0: mean = 0) ---")
            print(f"Delta x: n = {n_dx}, mean = {mean_dx:.4f}, p = {p_dx:.3e}")
            print(f"Delta y: n = {n_dy}, mean = {mean_dy:.4f}, p = {p_dy:.3e}")

            if reject:
                dx0 = mean_dx
                dy0 = mean_dy
                shift_applied_final = True
                print("Null hypothesis rejected and shifting carried over")
                print(f"Applied position offsets: Delta X0 = {dx0:.4f} cm, Delta Y0 = {dy0:.4f} cm")
            else:
                print("Null hypothesis not rejected and shifting not carried over")

    # PASS 2: final filling.
    eligible, matched, ghosts, _, _ = loop(
        apply_shift=shift_applied_final,
        dx0=dx0,
        dy0=dy0,
        fill_hists=True,
        collect_residuals=False,
    )

    # Write ROOT output.
    out_root = Path(f"tracklets_efficiency_{args.tag}.root")
    rf_out = ROOT.TFile(str(out_root), "RECREATE")
    h_den.Write()
    h_num.Write()
    h_eff = compute_eff_hist(h_num, h_den, "h_eff")
    h_eff.Write()

    for hist in diagnostic_hists.values():
        hist.Write()

    rf_out.Close()

    # Optional CSV.
    if args.write_csv:
        csv_path = Path("tracklets_full.csv")
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "entry",
                    "xA",
                    "yA",
                    "x_pred",
                    "y_pred",
                    "xT",
                    "yT",
                    "dx",
                    "dy",
                    "dt",
                    "dx_raw",
                    "dy_raw",
                    "dx0",
                    "dy0",
                    "theta_x_deg",
                    "theta_y_deg",
                    "theta_deg",
                ]
            )
            w.writerows(csv_rows)

    # Summary printout.
    eff = (matched / eligible) if eligible else 0.0
    err = math.sqrt(eff * (1.0 - eff) / eligible) if eligible else 0.0

    stamp(t0, f"ROOT output: {out_root}")
    if args.write_csv:
        stamp(t0, "CSV output: tracklets_full.csv")
    if args.write_ndjson:
        stamp(t0, "NDJSON output: tracklets.ndjson")

    print(f"\nEligible events: {eligible}")
    print(f"Tracklets found: {matched}")
    print(f"Ghosts found: {ghosts}")
    print(f"Global efficiency = {eff:.4f} +/- {err:.4f}")
    print(f"Regime = {args.regime}")
    print(f"Tree = {args.tree}")
    print(f"Diagnostic histograms = {'no' if args.no_diagnostic_hists else 'yes'}")
    if args.theta_max_deg is not None:
        print(f"Theta max = {args.theta_max_deg:.4f} deg")
    if eligible:
        print(f"Ghost rate = {ghosts / eligible:.4f}")

    for rf in roots:
        rf.Close()

    stamp(t0, "Done")


if __name__ == "__main__":
    main()
