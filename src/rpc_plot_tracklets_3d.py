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
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import ROOT


def in_bounds(x, y, rect):
    """Check whether a point lies inside the detector rectangle."""
    x_min, x_max, y_min, y_max = rect
    return (x_min <= x <= x_max) and (y_min <= y <= y_max)


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


def track_angle_deg(p0, p1):
    """Angle with respect to the detector normal (z axis)."""
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    dz = p1[2] - p0[2]
    transverse = np.hypot(dx, dy)
    theta = np.degrees(np.arctan2(transverse, abs(dz) + 1e-12))
    return float(theta)


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
    ap.add_argument("--chamber-lw", type=float, default=1.0,
                    help="Line width for chamber contours (default: 1.0).")
    ap.add_argument("--clip-to-chamber", action="store_true",
                    help="Reject tracklets with any chamber point outside the detector rectangle.")
    ap.add_argument("--direction", choices=["down", "up"], default="down",
                    help="Arrow direction along the tracklet: 'down' for cosmic-like top-to-bottom, 'up' otherwise.")
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
            ax_fit, bx = fit_line(z_ref, x_ref)
            ay_fit, by = fit_line(z_ref, y_ref)
            x_pred = ax_fit * z_t + bx
            y_pred = ay_fit * z_t + by
        else:
            x_pred, y_pred = xA, yA

        dx = xT - x_pred
        dy = yT - y_pred
        dt = tT - tA

        if (abs(dx) > args.tol_x) or (abs(dy) > args.tol_y) or (abs(dt) > args.tol_t):
            continue

        order = np.argsort(zs)
        pts = [(hits[j][0], hits[j][1], zs[j]) for j in order]

        if args.clip_to_chamber and any(not in_bounds(px, py, args.chamber_rect) for px, py, _ in pts):
            continue

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

    fig = plt.figure(figsize=(10.2, 8.2))
    ax = fig.add_subplot(111, projection="3d")

    # --- Draw chamber contours and semi-transparent detector planes
    x_min, x_max, y_min, y_max = args.chamber_rect
    rect_x = [x_min, x_max, x_max, x_min, x_min]
    rect_y = [y_min, y_min, y_max, y_max, y_min]
    for zc in zs:
        rect_z = [zc] * len(rect_x)
        verts = [[(x_min, y_min, zc), (x_max, y_min, zc), (x_max, y_max, zc), (x_min, y_max, zc)]]
        plane = Poly3DCollection(verts, facecolors=(0.86, 0.90, 0.96, 0.18), edgecolors="none")
        ax.add_collection3d(plane)
        ax.plot(rect_x, rect_y, rect_z, args.chamber_style, linewidth=args.chamber_lw, alpha=0.95, zorder=1)

    thetas = [track_angle_deg(tr["pts"][0], tr["pts"][-1]) for tr in tracks]
    theta_max = max(5.0, np.percentile(thetas, 98))
    norm = colors.Normalize(vmin=0.0, vmax=theta_max)
    cmap = plt.get_cmap("viridis")

    for tr, theta in zip(tracks, thetas):
        pts = tr["pts"]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        zz = [p[2] for p in pts]

        c = cmap(norm(theta))  # blue/green = near-vertical tracks, yellow = more inclined tracks

        # Faint skeleton of the whole fitted/matched trajectory.
        ax.plot(xs, ys, zz, linewidth=0.75, alpha=0.18, color="0.2", solid_capstyle="round", zorder=2)

        p_low = pts[0]
        p_high = pts[-1]
        start = p_high if args.direction == "down" else p_low
        end = p_low if args.direction == "down" else p_high
        mid = (
            0.5 * (start[0] + end[0]),
            0.5 * (start[1] + end[1]),
            0.5 * (start[2] + end[2]),
        )

        # Arrow with head at the midpoint to convey direction without hiding the detector planes.
        ax.quiver(
            start[0], start[1], start[2],
            mid[0] - start[0], mid[1] - start[1], mid[2] - start[2],
            arrow_length_ratio=0.18, linewidth=1.35, color=c, alpha=0.98,
            normalize=False,
        )
        # Continuation after the arrowhead so the full tracklet remains visible.
        ax.plot([mid[0], end[0]], [mid[1], end[1]], [mid[2], end[2]],
                linewidth=1.35, alpha=0.98, color=c, solid_capstyle="round", zorder=3)

        xa, ya, za = tr["anchor"]
        xt, yt, zt = tr["target"]
        ax.scatter([xa], [ya], [za], s=18, color=c, edgecolors="white", linewidths=0.35,
                   depthshade=False, alpha=0.95, zorder=4)
        ax.scatter([xt], [yt], [zt], s=18, color=c, edgecolors="white", linewidths=0.35,
                   depthshade=False, alpha=0.95, zorder=4)

        if args.show_intermediate and len(pts) > 2:
            xi = xs[1:-1]
            yi = ys[1:-1]
            zi = zz[1:-1]
            ax.scatter(xi, yi, zi, s=8, color="0.55", alpha=0.38, depthshade=False, zorder=3)

    ax.set_xlabel("X", labelpad=10)
    ax.set_ylabel("Y", labelpad=10)
    ax.set_zlabel("Z", labelpad=8)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(min(zs), max(zs))
    try:
        ax.set_box_aspect((x_max - x_min, y_max - y_min, max(zs) - min(zs)))
    except AttributeError:
        pass

    ax.set_title(args.title or f"RPC tracklets (N={len(tracks)}) — {args.tag}", pad=16)
    ax.view_init(elev=args.elev, azim=args.azim)

    # Cleaner "publication-like" styling.
    ax.grid(True, alpha=0.22, linewidth=0.6)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
            axis.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
        except AttributeError:
            pass
    ax.tick_params(axis="both", which="major", labelsize=10, pad=3)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.06, shrink=0.84)
    cbar.set_label(r"Track inclination $\theta$ [deg]", rotation=90, labelpad=12)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    if args.pdf:
        fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Matched tracklets plotted: {len(tracks)} (scanned {scanned} events)")
    print(f"Arrow direction: {args.direction}")
    print(f"Saved: {out_png}")
    if args.pdf:
        print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
