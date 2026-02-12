#!/usr/bin/env python3
"""
rpc_plot_tracklets_3d.py

3D visualization of matched tracklets as straight lines through multiple RPC chambers.

What it does
------------
- Opens N ROOT files (one per chamber), each containing the same TTree (same event order).
- For each event, picks one hit per chamber (earliest or centroid) from branches:
    Strip_X, Strip_Y, Strip_Tdc_diff
- Builds a multi-chamber "tracklet" by predicting the target chamber hit from the
  reference chambers (all chambers except target). With >=2 reference planes it fits
  x(z)=a_x*z+b_x and y(z)=a_y*z+b_y.
- Applies matching cuts on the target residuals (dx, dy) and timing difference (dt).
- For each matched event, plots:
    * A thin black polyline connecting all chamber points (ordered by z)
    * Anchor and Target points emphasized and colored identically per track
    * Intermediate chamber points in light gray (optional)

Outputs
-------
- Saves a PNG (and optionally PDF) with a 3D plot.

Examples
--------
2 chambers:
  python3 rpc_plot_tracklets_3d.py --files A.root B.root --z 0 100 --tag HV6

3 chambers:
  python3 rpc_plot_tracklets_3d.py --files A.root B.root C.root --z 0 50 100 --tag HV6

Note
----
If you produced NDJSON/CSV already, you can still use this tool directly from ROOT files
so it always has per-chamber points.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import ROOT


def pick_hit(ev, hit_mode: str = "earliest"):
    """Pick one hit from vector branches: Strip_X, Strip_Y, Strip_Tdc_diff."""
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


def main():
    ap = argparse.ArgumentParser(description="3D plot of RPC tracklets through N chambers.")
    ap.add_argument("--files", nargs="+", required=True, help="ROOT files, one per chamber (same event order).")
    ap.add_argument("--z", nargs="+", type=float, required=True, help="z positions (same count/order as --files).")
    ap.add_argument("--tree", default="events", help="TTree name (default: events)")
    ap.add_argument("--tag", default="tracklets", help="Tag used in output filenames.")

    ap.add_argument("--anchor-index", type=int, default=0, help="Anchor chamber index (default: 0).")
    ap.add_argument("--target-index", type=int, default=-1, help="Target chamber index (default: -1=last).")

    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])
    ap.add_argument("--require-ncluster-eq1", action="store_true",
                    help="Require nCluster==1 for ALL chambers (strict).")

    ap.add_argument("--tol-x", type=float, default=7.0, help="|dx| tolerance (same units as Strip_X)")
    ap.add_argument("--tol-y", type=float, default=16.0, help="|dy| tolerance (same units as Strip_Y)")
    ap.add_argument("--tol-t", type=float, default=3.0, help="|dt| tolerance (same units as Strip_Tdc_diff)")

    ap.add_argument("--max-tracks", type=int, default=200, help="Maximum number of tracks to draw.")
    ap.add_argument("--max-events", type=int, default=None, help="Maximum events to scan (default: all).")

    ap.add_argument("--out", default=None, help="Output PNG filename (default: tracklets_3d_<tag>.png)")
    ap.add_argument("--pdf", action="store_true", help="Also save a PDF.")
    ap.add_argument("--show-intermediate", action="store_true",
                    help="Show intermediate points (light gray). Default: only endpoints colored.")

    ap.add_argument("--title", default=None, help="Plot title.")
    ap.add_argument("--elev", type=float, default=18.0, help="3D view elevation angle.")
    ap.add_argument("--azim", type=float, default=-60.0, help="3D view azimuth angle.")
    ap.add_argument("--chamber-rect", nargs=4, type=float, default=[0.0, 54.0, 0.0, 120.0],
                    metavar=("XMIN","XMAX","YMIN","YMAX"),
                    help="Chamber contour as rectangle in X,Y (default: 0 54 0 120).")
    ap.add_argument("--chamber-style", default="k",
                    help="Matplotlib style/color for chamber contour lines (default: k).")
    ap.add_argument("--chamber-lw", type=float, default=0.8,
                    help="Line width for chamber contours (default: 0.8).")
    args = ap.parse_args()

    files = [Path(p) for p in args.files]
    zs = list(args.z)
    if len(files) != len(zs):
        raise SystemExit("ERROR: --files and --z must have the same number of entries.")
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

    # Open files and trees
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
    if args.max_events is not None:
        n_entries = min(n_entries, int(args.max_events))

    tracks = []
    scanned = 0

    ref_idx = [i for i in range(n_ch) if i != tgt]

    for i in range(n_entries):
        scanned += 1
        for t in trees:
            t.GetEntry(i)

        if args.require_ncluster_eq1:
            if any(int(t.nCluster) != 1 for t in trees):
                continue

        hits = {}
        ok = True
        for j in range(n_ch):
            ev = trees[j]
            hit = pick_hit(ev, args.hit_mode)
            if hit is None:
                ok = False
                break
            hits[j] = hit  # (x,y,t)
        if not ok:
            continue

        xA, yA, tA = hits[anc]
        xT, yT, tT = hits[tgt]

        # Predict target from reference chambers (excluding target)
        z_ref = [zs[j] for j in ref_idx]
        x_ref = [hits[j][0] for j in ref_idx]
        y_ref = [hits[j][1] for j in ref_idx]

        z_t = zs[tgt]
        if len(z_ref) >= 2:
            ax, bx = fit_line(z_ref, x_ref)
            ay, by = fit_line(z_ref, y_ref)
            x_pred = ax * z_t + bx
            y_pred = ay * z_t + by
        else:
            x_pred, y_pred = xA, yA

        dx = xT - x_pred
        dy = yT - y_pred
        dt = tT - tA

        if (abs(dx) > args.tol_x) or (abs(dy) > args.tol_y) or (abs(dt) > args.tol_t):
            continue

        order = np.argsort(zs)
        pts = [(hits[j][0], hits[j][1], zs[j]) for j in order]

        tracks.append({
            "entry": i,
            "anchor": (xA, yA, zs[anc]),
            "target": (xT, yT, zs[tgt]),
            "pts": pts
        })
        if len(tracks) >= args.max_tracks:
            break

    for rf in roots:
        rf.Close()

    if not tracks:
        raise SystemExit("No matched tracklets found. Try loosening tolerances or increasing max-events.")

    out_png = args.out or f"tracklets_3d_{args.tag}.png"
    out_pdf = out_png.replace(".png", ".pdf")

    fig = plt.figure(figsize=(9.0, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    # --- Draw chamber contours (one per z plane)
    x_min, x_max, y_min, y_max = args.chamber_rect
    rect_x = [x_min, x_max, x_max, x_min, x_min]
    rect_y = [y_min, y_min, y_max, y_max, y_min]
    for zc in zs:
        rect_z = [zc] * len(rect_x)
        ax.plot(rect_x, rect_y, rect_z, args.chamber_style, linewidth=args.chamber_lw, alpha=0.9)

    cmap = plt.get_cmap("tab20")
    ncol = cmap.N

    for k, tr in enumerate(tracks):
        pts = tr["pts"]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        zz = [p[2] for p in pts]

        ax.plot(xs, ys, zz, linewidth=0.7, alpha=0.8, color="black")

        c = cmap(k % ncol)
        xa, ya, za = tr["anchor"]
        xt, yt, zt = tr["target"]
        ax.scatter([xa], [ya], [za], s=35, color=c, depthshade=False)
        ax.scatter([xt], [yt], [zt], s=35, color=c, depthshade=False)

        if args.show_intermediate and len(pts) > 2:
            xi = xs[1:-1]
            yi = ys[1:-1]
            zi = zz[1:-1]
            ax.scatter(xi, yi, zi, s=10, color="0.65", depthshade=False)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_title(args.title or f"RPC tracklets (N={len(tracks)}) — {args.tag}")
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.grid(True)

    fig.tight_layout()
    fig.savefig(out_png, dpi=220)
    if args.pdf:
        fig.savefig(out_pdf)
    plt.close(fig)

    print(f"Matched tracklets plotted: {len(tracks)} (scanned {scanned} events)")
    print(f"Saved: {out_png}")
    if args.pdf:
        print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
