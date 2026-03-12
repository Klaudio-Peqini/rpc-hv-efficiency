#!/usr/bin/env python3
"""
rpc_estimate_tracklet_cuts.py

Estimate per-HV matching cuts for RPC tracklet reconstruction.

This script performs a loose preselection, computes target residuals relative to
reference-track prediction, and derives robust cut values per HV using the
median/MAD estimator. It is designed to plug into hv_pipeline_nch_extended_regime.py.
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import ROOT

HV_ORDER = {f"HV{i}": i for i in range(1, 100)}


def stamp(t0, msg):
    print(f"[{time.time() - t0:8.2f} s] {msg}")


def pick_hit(ev, hit_mode="earliest"):
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
    z = np.asarray(zs, dtype=float)
    v = np.asarray(vs, dtype=float)
    A = np.vstack([z, np.ones_like(z)]).T
    (a, b), *_ = np.linalg.lstsq(A, v, rcond=None)
    return float(a), float(b)


def robust_center_sigma(values):
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return 0.0, 0.0
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    sigma = 1.4826 * mad
    if sigma == 0.0 and x.size >= 2:
        sigma = float(np.std(x, ddof=1))
    return med, sigma


def capped_abs_cut(values, nsigma, min_cut, max_cut, percentile=0.95):
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return float(min_cut), 0.0, 0.0, 0.0
    center, sigma = robust_center_sigma(x)
    cut_robust = nsigma * sigma
    cut_quant = float(np.quantile(np.abs(x - center), percentile)) if x.size >= 5 else cut_robust
    raw_cut = max(cut_robust, cut_quant)
    cut = min(max(raw_cut, min_cut), max_cut)
    return float(cut), float(center), float(sigma), float(raw_cut)


def angle_cut(dx_vals, dy_vals, dz_span, min_deg=0.0, max_deg=25.0, percentile=0.98):
    if dz_span <= 0 or len(dx_vals) == 0:
        return float(max_deg)
    dx = np.asarray(dx_vals, dtype=float)
    dy = np.asarray(dy_vals, dtype=float)
    theta = np.degrees(np.arctan2(np.sqrt(dx * dx + dy * dy), dz_span))
    q = float(np.quantile(theta, percentile)) if theta.size >= 5 else float(np.max(theta))
    return min(max(q, min_deg), max_deg)


def iter_hv_tags(base: Path):
    return sorted([d.name for d in base.glob("HV*") if d.is_dir()], key=lambda s: (HV_ORDER.get(s, 9999), s))


def load_tree(path: Path, tree_name: str):
    rf = ROOT.TFile.Open(str(path))
    if not rf or rf.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    tr = rf.Get(tree_name)
    if not tr:
        raise RuntimeError(f"Tree '{tree_name}' not found in {path}")
    return rf, tr


def collect_for_hv(files, zs, args):
    roots, trees = [], []
    for f in files:
        rf, tr = load_tree(Path(f), args.tree)
        roots.append(rf)
        trees.append(tr)

    n_entries = min(int(t.GetEntries()) for t in trees)
    n_ch = len(files)
    tgt = args.target_index if args.target_index >= 0 else n_ch + args.target_index
    anc = args.anchor_index if args.anchor_index >= 0 else n_ch + args.anchor_index
    ref_idx = [i for i in range(n_ch) if i != tgt]

    dx_list, dy_list, dt_list = [], [], []
    xpred_list, ypred_list = [], []
    xA_list, yA_list = [], []
    mx_list, my_list = [], []

    z_ref_all = [zs[j] for j in ref_idx]
    dz_span = float(max(z_ref_all) - min(z_ref_all)) if len(z_ref_all) >= 2 else float(abs(zs[tgt] - zs[anc]))

    for i in range(n_entries):
        for t in trees:
            t.GetEntry(i)

        evA = trees[anc]
        evT = trees[tgt]

        if int(evA.nCluster) != 1:
            continue
        if int(evT.nCluster) > 1:
            continue
        if len(evA.Strip_Tdc_diff) > args.cluster_size_max:
            continue

        hits_ref = {}
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

        hitA = hits_ref.get(anc)
        hitT = pick_hit(evT, args.hit_mode)
        if hitA is None or hitT is None:
            continue

        xA, yA, tA = hitA
        xT, yT, tT = hitT

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
            mx, my = ax, ay
        else:
            x_pred, y_pred = xA, yA
            dz = zs[tgt] - zs[anc]
            mx = (xT - xA) / dz if dz != 0 else 0.0
            my = (yT - yA) / dz if dz != 0 else 0.0

        if args.regime == "beam" and args.roi_rect is not None:
            xmin, xmax, ymin, ymax = args.roi_rect
            if not (xmin <= x_pred <= xmax and ymin <= y_pred <= ymax):
                continue

        dx_list.append(float(xT - x_pred))
        dy_list.append(float(yT - y_pred))
        dt_list.append(float(tT - tA))
        xpred_list.append(float(x_pred))
        ypred_list.append(float(y_pred))
        xA_list.append(float(xA))
        yA_list.append(float(yA))
        mx_list.append(float(mx))
        my_list.append(float(my))

    for rf in roots:
        rf.Close()

    return {
        "dx": dx_list,
        "dy": dy_list,
        "dt": dt_list,
        "x_pred": xpred_list,
        "y_pred": ypred_list,
        "x_anchor": xA_list,
        "y_anchor": yA_list,
        "mx": mx_list,
        "my": my_list,
        "dz_span": dz_span,
        "n_samples": len(dx_list),
    }


def build_hv_file_list(chamber_bases, hvtag):
    files = []
    for base in chamber_bases:
        f = Path(base) / hvtag / "data.root"
        if not f.exists():
            return None
        files.append(str(f))
    return files


def main():
    t0 = time.time()
    ap = argparse.ArgumentParser(description="Estimate per-HV tracklet matching cuts.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--regime", choices=["normal", "beam"], default="normal")
    ap.add_argument("--tree", default="events")
    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])
    ap.add_argument("--mode", default="optionc", choices=["v2", "optionc"])
    ap.add_argument("--ref-base")
    ap.add_argument("--test-base")
    ap.add_argument("--chambers", nargs="+")
    ap.add_argument("--z", nargs="+", type=float)
    ap.add_argument("--anchor-index", type=int, default=0)
    ap.add_argument("--target-index", type=int, default=-1)
    ap.add_argument("--require-ncluster-eq1", action="store_true")
    ap.add_argument("--cluster-size-max", type=int, default=6)
    ap.add_argument("--cluster-size-on-all-ref", action="store_true")
    ap.add_argument("--roi-name", default="scint_overlap")
    ap.add_argument("--roi-rect", nargs=4, type=float)
    ap.add_argument("--dx-nsigma", type=float, default=3.0)
    ap.add_argument("--dy-nsigma", type=float, default=3.0)
    ap.add_argument("--dt-nsigma", type=float, default=3.0)
    ap.add_argument("--min-samples", type=int, default=50)
    ap.add_argument("--min-tol-x", type=float, default=3.0)
    ap.add_argument("--max-tol-x", type=float, default=15.0)
    ap.add_argument("--min-tol-y", type=float, default=6.0)
    ap.add_argument("--max-tol-y", type=float, default=35.0)
    ap.add_argument("--min-tol-t", type=float, default=1.0)
    ap.add_argument("--max-tol-t", type=float, default=8.0)
    ap.add_argument("--angle-max-deg", type=float, default=25.0)
    ap.add_argument("--percentile", type=float, default=0.95)
    args = ap.parse_args()

    if args.chambers:
        chamber_bases = [Path(p) for p in args.chambers]
    else:
        if not args.ref_base or not args.test_base:
            raise SystemExit("Provide either --chambers or both --ref-base and --test-base")
        chamber_bases = [Path(args.ref_base), Path(args.test_base)]

    if args.z is None:
        if len(chamber_bases) == 2:
            zs = [0.0, 100.0]
        else:
            raise SystemExit("For N-ch mode, provide --z")
    else:
        zs = list(args.z)

    if len(zs) != len(chamber_bases):
        raise SystemExit("--z must have one value per chamber")

    hvtags = iter_hv_tags(chamber_bases[0])
    out = {}
    diag = {}

    stamp(t0, f"Estimating cuts for regime={args.regime} across {len(hvtags)} HV points")
    for hvtag in hvtags:
        files = build_hv_file_list(chamber_bases, hvtag)
        if files is None:
            stamp(t0, f"Skipping {hvtag}: missing one or more data.root files")
            continue

        stats = collect_for_hv(files, zs, args)
        n = stats["n_samples"]
        dx = np.asarray(stats["dx"], dtype=float)
        dy = np.asarray(stats["dy"], dtype=float)
        dt = np.asarray(stats["dt"], dtype=float)

        if n < args.min_samples:
            tol_x = max(args.min_tol_x, min(args.max_tol_x, 7.0))
            tol_y = max(args.min_tol_y, min(args.max_tol_y, 20.0))
            tol_t = max(args.min_tol_t, min(args.max_tol_t, 3.0))
            cx, sx = robust_center_sigma(dx)
            cy, sy = robust_center_sigma(dy)
            ct, st = robust_center_sigma(dt)
            theta_max = angle_cut(dx, dy, stats["dz_span"], max_deg=args.angle_max_deg)
        else:
            tol_x, cx, sx, raw_x = capped_abs_cut(dx, args.dx_nsigma, args.min_tol_x, args.max_tol_x, args.percentile)
            tol_y, cy, sy, raw_y = capped_abs_cut(dy, args.dy_nsigma, args.min_tol_y, args.max_tol_y, args.percentile)
            tol_t, ct, st, raw_t = capped_abs_cut(dt, args.dt_nsigma, args.min_tol_t, args.max_tol_t, args.percentile)
            theta_max = angle_cut(dx, dy, stats["dz_span"], max_deg=args.angle_max_deg)

        out[hvtag] = {
            "tol_x": round(float(tol_x), 4),
            "tol_y": round(float(tol_y), 4),
            "tol_t": round(float(tol_t), 4),
            "dx_center": round(float(cx), 4),
            "dy_center": round(float(cy), 4),
            "dt_center": round(float(ct), 4),
            "theta_max_deg": round(float(theta_max), 4),
            "n_samples": int(n),
            "regime": args.regime,
            "method": "robust_mad_plus_quantile",
        }
        diag[hvtag] = {
            "sigma_dx": round(float(sx), 4),
            "sigma_dy": round(float(sy), 4),
            "sigma_dt": round(float(st), 4),
            "dz_span": round(float(stats["dz_span"]), 4),
            "n_samples": int(n),
        }
        stamp(t0, f"{hvtag}: n={n} -> tol_x={tol_x:.2f}, tol_y={tol_y:.2f}, tol_t={tol_t:.2f}, theta_max={theta_max:.1f} deg")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    diag_path = out_path.with_name(out_path.stem + "_diagnostics.json")
    with diag_path.open("w", encoding="utf-8") as f:
        json.dump(diag, f, indent=2, sort_keys=True)

    stamp(t0, f"Wrote cuts JSON: {out_path}")
    stamp(t0, f"Wrote diagnostics JSON: {diag_path}")


if __name__ == "__main__":
    main()
