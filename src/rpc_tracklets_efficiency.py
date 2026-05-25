#!/usr/bin/env python3
"""
tracklets_reconstruction_efficiency_nch_optionC.py

Option C-final N-chamber producer:
- Keeps *v2* eligibility semantics (denominator definition & preselection style)
- Uses multi-chamber prediction (>=2 reference planes -> line fit; 2-ch fallback uses anchor hit)
- Includes automatic Δx/Δy shifting by hypothesis test (two-pass), applied to matching windows
- Produces the same core ROOT outputs used downstream:
    - h_total (denominator) filled in anchor coordinates (xA,yA)
    - h_tracklet (numerator) filled in anchor coordinates (xA,yA)
    - h_eff = h_tracklet / h_total

Notes:
- "Eligible" counter matches v2 logic: incremented after clean preselection,
  even if later hit extraction fails (rare), exactly like v2.
- Denominator histogram (h_total) is filled only after anchor hit is available,
  as in v2.
"""

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
import ROOT


def stamp(t0, msg):
    print(f"[{time.time() - t0:8.2f} s] {msg}")


def pick_hit(ev, hit_mode="earliest"):
    """
    Pick one hit from vector branches:
      Strip_X, Strip_Y, Strip_Tdc_diff
    """
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


def fit_line(zs, vs):
    """Least squares fit v(z)=a*z+b for >=2 points."""
    z = np.asarray(zs, dtype=float)
    v = np.asarray(vs, dtype=float)
    A = np.vstack([z, np.ones_like(z)]).T
    (a, b), *_ = np.linalg.lstsq(A, v, rcond=None)
    return float(a), float(b)


def compute_eff_hist(h_num, h_den, name="h_eff"):
    h_eff = h_num.Clone(name)
    h_eff.Reset()
    h_eff.Divide(h_num, h_den, 1.0, 1.0, "B")
    return h_eff


def ttest_mean_zero(samples):
    """
    One-sample t-test for H0: mean=0 (two-sided).
    Uses SciPy if available; otherwise falls back to normal approximation (good for large N).
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


def main():
    ap = argparse.ArgumentParser(
        description="Option C-final: v2-eligible N-ch efficiency + auto Δx/Δy shift"
    )
    ap.add_argument("--files", nargs="+", required=True,
                    help="ROOT files, one per chamber, same event order (HV point).")
    ap.add_argument("--z", nargs="+", required=True, type=float,
                    help="z positions (same count/order as --files).")
    ap.add_argument("--tree", default="events")
    ap.add_argument("--tag", required=True, help="e.g. HV6")

    ap.add_argument("--target-index", type=int, default=-1,
                    help="Target chamber index (default: last).")
    ap.add_argument("--anchor-index", type=int, default=0,
                    help="Anchor chamber defining denominator map & reference time (default: 0).")

    # v2-like selection knobs (Option C)
    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])

    ap.add_argument("--v2-eligibility", action="store_true",
                    help="Enforce v2 eligibility semantics (recommended / Option C).")
    ap.add_argument("--require-ncluster-eq1", action="store_true",
                    help="Require nCluster==1 for ALL reference chambers (excluding target).")
    ap.add_argument("--reject-target-ncluster-gt1", action="store_true",
                    help="Reject events with target nCluster>1 BEFORE eligible++ (v2 semantics).")

    ap.add_argument("--cluster-size-max", type=int, default=999999,
                    help="Cluster size veto threshold (len(Strip_Tdc_diff)).")
    ap.add_argument("--cluster-size-on-all-ref", action="store_true",
                    help="Apply cluster-size veto to ALL reference chambers (default: anchor only, like v2).")

    # Matching tolerances
    ap.add_argument("--tol-x", type=float, default=2.0)
    ap.add_argument("--tol-y", type=float, default=2.0)
    ap.add_argument("--tol-t", type=float, default=4.0)

    # Shift / hypothesis-testing
    ap.add_argument("--alpha", type=float, default=0.01,
                    help="Significance level for rejecting H0: mean(Δ)=0. Default: 0.01.")
    ap.add_argument("--disable-shift", action="store_true",
                    help="Disable automatic Δx/Δy shifting (debug).")
    ap.add_argument("--min-residuals", type=int, default=50,
                    help="Minimum number of residual samples needed to run hypothesis test.")

    # Histogram config
    ap.add_argument("--nbins-x", type=int, default=64)
    ap.add_argument("--nbins-y", type=int, default=128)
    ap.add_argument("--x-range", nargs=2, type=float, default=[0.0, 60.0])
    ap.add_argument("--y-range", nargs=2, type=float, default=[0.0, 140.0])

    # Outputs
    ap.add_argument("--write-csv", action="store_true")
    ap.add_argument("--write-ndjson", action="store_true")

    args = ap.parse_args()
    t0 = time.time()
    stamp(t0, "Starting Option C-final N-ch efficiency producer")

    files = [Path(p) for p in args.files]
    zs = list(args.z)
    if len(files) != len(zs):
        raise SystemExit("ERROR: --files and --z must have same length.")
    n_ch = len(files)
    if n_ch < 2:
        raise SystemExit("ERROR: Need at least 2 chambers.")

    tgt = args.target_index
    if tgt < 0:
        tgt = n_ch + tgt
    anc = args.anchor_index
    if anc < 0:
        anc = n_ch + anc
    if not (0 <= tgt < n_ch) or not (0 <= anc < n_ch):
        raise SystemExit("ERROR: invalid target-index or anchor-index.")
    if tgt == anc:
        raise SystemExit("ERROR: target-index cannot equal anchor-index.")

    ref_idx = [i for i in range(n_ch) if i != tgt]  # reference set includes anchor
    if anc not in ref_idx:
        raise SystemExit("ERROR: internal: anchor must be in reference set.")

    # Open files / trees
    roots, trees = [], []
    for f in files:
        rf = ROOT.TFile.Open(str(f))
        if not rf or rf.IsZombie():
            raise SystemExit(f"ERROR: cannot open {f}")
        tr = rf.Get(args.tree)
        if not tr:
            raise SystemExit(f"ERROR: tree '{args.tree}' not found in {f}")
        roots.append(rf)
        trees.append(tr)

    n_entries = min(int(t.GetEntries()) for t in trees)
    stamp(t0, f"Loaded {n_ch} trees, using first {n_entries} entries")

    # Histograms (anchor coords, like v2)
    h_den = ROOT.TH2D(
        "h_total", "Denominator; x_anchor; y_anchor",
        args.nbins_x, args.x_range[0], args.x_range[1],
        args.nbins_y, args.y_range[0], args.y_range[1]
    )
    h_num = ROOT.TH2D(
        "h_tracklet", "Numerator; x_anchor; y_anchor",
        args.nbins_x, args.x_range[0], args.x_range[1],
        args.nbins_y, args.y_range[0], args.y_range[1]
    )

    # Record outputs only on the final pass (fill_hists=True)
    csv_rows = []
    ndjson_path = Path("tracklets.ndjson")
    if args.write_ndjson:
        ndjson_path.write_text("")

    # Option C defaults: enforce v2 eligibility unless user explicitly chooses otherwise

    def loop(apply_shift: bool, dx0: float, dy0: float, fill_hists: bool, collect_residuals: bool):
        eligible = matched = ghosts = 0
        dx_list, dy_list = [], []

        for i in range(n_entries):
            for t in trees:
                t.GetEntry(i)

            evA = trees[anc]
            evT = trees[tgt]

            # --- v2 eligibility preselection (Option C) ---
            # v2: if evA.nCluster != 1 or evB.nCluster > 1: continue
            # We enforce anchor clean always. Target>1 is rejected *before* eligible++.
            if int(evA.nCluster) != 1:
                continue

            if int(evT.nCluster) > 1:
                # v2 semantics: reject before eligible++
                # (also matches the user's observed discrepancy source)
                continue

            # v2: cluster size veto is applied on anchor
            if len(evA.Strip_Tdc_diff) > args.cluster_size_max:
                continue

            # Additional reference chambers (for prediction): require hits (and optionally nCluster==1)
            hits_ref = {}  # idx -> (x,y,t) for reference chambers only (exclude target)
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

            # v2: eligible++ happens after the clean preselection
            eligible += 1

            # Anchor hit (may be None in pathological cases; v2 keeps eligible but doesn't fill denom)
            hitA = hits_ref.get(anc)
            if hitA is None:
                continue
            xA, yA, tA = hitA

            if fill_hists:
                h_den.Fill(xA, yA)  # denom map at anchor coords (like v2)

            # v2: ghost check happens AFTER eligible++ AND AFTER denominator fill (DENOM must not depend on target)
            if int(evT.nCluster) == 0:
                ghosts += 1
                continue


            # --- Predict target from references (multi-ch) ---
            # Use all available ref planes (excluding target). With >=2 planes we fit; with 1 we fallback.
            z_ref, x_ref, y_ref = [], [], []
            for j in ref_idx:
                z_ref.append(zs[j])
                x_ref.append(hits_ref[j][0])
                y_ref.append(hits_ref[j][1])

            z_t = zs[tgt]
            if len(z_ref) >= 2:
                ax, bx = fit_line(z_ref, x_ref)
                ay, by = fit_line(z_ref, y_ref)
                x_pred = ax * z_t + bx
                y_pred = ay * z_t + by
            else:
                x_pred, y_pred = xA, yA

            # v2-like timing anchor
            t_pred = tA

            # Target hit
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

            # Matching window
            if (abs(dt) < args.tol_t) and (abs(dx) < args.tol_x) and (abs(dy) < args.tol_y):
                matched += 1

                if fill_hists:
                    h_num.Fill(xA, yA)

                    if args.write_csv:
                        csv_rows.append(
                            (i, xA, yA, x_pred, y_pred, xT, yT, dx, dy, dt, dx_raw, dy_raw, dx0, dy0)
                        )

                    if args.write_ndjson:
                        rec = {
                            "entry": i,
                            "anchor_index": anc,
                            "target_index": tgt,
                            "n_chambers": n_ch,
                            "xA": xA, "yA": yA, "tA": tA,
                            "x_pred": x_pred, "y_pred": y_pred, "t_pred": t_pred,
                            "xT": xT, "yT": yT, "tT": tT,
                            "dx": dx, "dy": dy, "dt": dt,
                            "dx_raw": dx_raw, "dy_raw": dy_raw,
                            "shift_applied": bool(apply_shift),
                            "dx0": dx0, "dy0": dy0,
                        }
                        with ndjson_path.open("a") as f:
                            f.write(json.dumps(rec) + "\n")

                if collect_residuals:
                    dx_list.append(dx_raw)
                    dy_list.append(dy_raw)

        return eligible, matched, ghosts, dx_list, dy_list

    # --- PASS 1: collect residuals for hypothesis test ---
    dx0 = dy0 = 0.0
    if args.disable_shift:
        print("Null hypothesis not rejected and shifting not carried over")
    else:
        e1, m1, g1, dx_list, dy_list = loop(
            apply_shift=False, dx0=0.0, dy0=0.0, fill_hists=False, collect_residuals=True
        )

        if (len(dx_list) < args.min_residuals) or (len(dy_list) < args.min_residuals):
            # Not enough stats to make a decision
            print("Null hypothesis not rejected and shifting not carried over")
        else:
            mean_dx, p_dx, _ = ttest_mean_zero(dx_list)
            mean_dy, p_dy, _ = ttest_mean_zero(dy_list)

            reject = (p_dx < args.alpha) or (p_dy < args.alpha)
            if reject:
                dx0, dy0 = mean_dx, mean_dy
                print("Null hypothesis rejected and shifting carried over")
                print(f"Applied position offsets: ΔX₀ = {dx0:.4f} cm, ΔY₀ = {dy0:.4f} cm")
            else:
                print("Null hypothesis not rejected and shifting not carried over")

    # --- PASS 2: final filling (shift applied if decided) ---
    eligible, matched, ghosts, _, _ = loop(
        apply_shift=(not args.disable_shift and (dx0 != 0.0 or dy0 != 0.0)),
        dx0=dx0, dy0=dy0,
        fill_hists=True,
        collect_residuals=False
    )

    # Write ROOT output
    out_root = Path(f"tracklets_efficiency_{args.tag}.root")
    rf_out = ROOT.TFile(str(out_root), "RECREATE")
    h_den.Write()
    h_num.Write()
    h_eff = compute_eff_hist(h_num, h_den, "h_eff")
    h_eff.Write()
    rf_out.Close()

    # Optional CSV
    if args.write_csv:
        csv_path = Path("tracklets_full.csv")
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "entry", "xA", "yA", "x_pred", "y_pred", "xT", "yT",
                "dx", "dy", "dt", "dx_raw", "dy_raw", "dx0", "dy0"
            ])
            w.writerows(csv_rows)

    # Summary print (match v2 style)
    eff = (matched / eligible) if eligible else 0.0
    err = math.sqrt(eff * (1 - eff) / eligible) if eligible else 0.0

    stamp(t0, f"ROOT output:   {out_root}")
    if args.write_csv:
        stamp(t0, "CSV output:    tracklets_full.csv")
    if args.write_ndjson:
        stamp(t0, "NDJSON output: tracklets.ndjson")

    print(f"\nEligible events: {eligible}")
    print(f"Tracklets found: {matched}")
    print(f"Ghosts found: {ghosts}")
    print(f"Global efficiency = {eff:.4f} ± {err:.4f}")
    if eligible:
        print(f"Ghost rate = {ghosts/eligible:.4f}")

    stamp(t0, "Done")


if __name__ == "__main__":
    main()
