#!/usr/bin/env python3
"""
rpc_wp_validation_plots.py

Fast WP-validation plotting utility for rpc-hv-efficiency outputs.
It produces:
  3. Delta x distribution
  4. Delta y distribution
  5. theta_x distribution
  6. theta_y distribution
  7. theta_x vs theta_y 2D map
  8. tracklet-level correlation matrix from the per-tracklet CSV

Expected ROOT histograms from the updated producer:
  h_dx, h_dy, h_theta_x, h_theta_y, h_theta_x_theta_y

Expected CSV columns, when --csv is given:
  xA, yA, x_pred, y_pred, xT, yT, dx_raw, dy_raw, dx, dy, dt,
  theta_x_deg, theta_y_deg, theta_deg

Example:
  python3 rpc_wp_validation_plots.py \
    --root out_dynamic_regime_36BX_diag/maps/tracklets_efficiency_HV6_beam.root \
    --csv  out_dynamic_regime_36BX_diag/maps/tracklets_full_HV6_beam.csv \
    --out-dir out_dynamic_regime_36BX_diag/plots/wp_validation/HV6_beam \
    --label $'Chamber 202, 36BX\n1.4 mm double gap iRPC\nFEB v2.3 Petiroc 2C\nthreshold ~ 40 fC' \
    --png --pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def cms_red_cmap():
    return LinearSegmentedColormap.from_list(
        "cms_white_red",
        ["#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"],
    )


def setup_style(dpi: int) -> None:
    plt.rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "axes.linewidth": 1.25,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "figure.constrained_layout.use": False,
    })


def add_cms_labels(ax: plt.Axes, label: str = "") -> None:
    ax.text(0.00, 1.015, "CMS", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=15, fontweight="bold")
    ax.text(0.115, 1.015, "Preliminary", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=13, style="italic")
    ax.text(1.00, 1.015, "GIF++", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=13)
    if label:
        ax.text(0.045, 0.93, label, transform=ax.transAxes,
                ha="left", va="top", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="none", alpha=0.78))


def save_all(fig: plt.Figure, out_base: Path, formats: Sequence[str]) -> List[Path]:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        out = out_base.with_suffix(f".{fmt}")
        fig.savefig(out, bbox_inches="tight")
        written.append(out)
    plt.close(fig)
    return written


def _list_with_uproot(root_path: Path) -> List[str]:
    import uproot  # type: ignore
    with uproot.open(root_path) as f:
        return [k.split(";")[0] for k in f.keys()]


def _read_th1_with_uproot(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray]:
    import uproot  # type: ignore
    with uproot.open(root_path) as f:
        if name not in f:
            keys = [k.split(";")[0] for k in f.keys()]
            raise KeyError(f"Histogram '{name}' not found in {root_path}. Available: {keys}")
        h = f[name]
        values = np.asarray(h.values(flow=False), dtype=float)
        edges = np.asarray(h.axes[0].edges(), dtype=float)
    return values, edges


def _read_th2_with_uproot(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    import uproot  # type: ignore
    with uproot.open(root_path) as f:
        if name not in f:
            keys = [k.split(";")[0] for k in f.keys()]
            raise KeyError(f"Histogram '{name}' not found in {root_path}. Available: {keys}")
        h = f[name]
        values = np.asarray(h.values(flow=False), dtype=float)
        x_edges = np.asarray(h.axes[0].edges(), dtype=float)
        y_edges = np.asarray(h.axes[1].edges(), dtype=float)
    return values, x_edges, y_edges


def _read_th1_with_pyroot(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray]:
    import ROOT  # type: ignore
    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")
    try:
        h = f.Get(name)
        if not h:
            keys = [k.GetName() for k in f.GetListOfKeys()]
            raise KeyError(f"Histogram '{name}' not found in {root_path}. Available: {keys}")
        n = int(h.GetNbinsX())
        ax = h.GetXaxis()
        edges = np.array([ax.GetBinLowEdge(i + 1) for i in range(n)] + [ax.GetBinUpEdge(n)], dtype=float)
        values = np.array([h.GetBinContent(i + 1) for i in range(n)], dtype=float)
    finally:
        f.Close()
    return values, edges


def _read_th2_with_pyroot(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    import ROOT  # type: ignore
    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")
    try:
        h = f.Get(name)
        if not h:
            keys = [k.GetName() for k in f.GetListOfKeys()]
            raise KeyError(f"Histogram '{name}' not found in {root_path}. Available: {keys}")
        nx = int(h.GetNbinsX())
        ny = int(h.GetNbinsY())
        ax = h.GetXaxis()
        ay = h.GetYaxis()
        x_edges = np.array([ax.GetBinLowEdge(i + 1) for i in range(nx)] + [ax.GetBinUpEdge(nx)], dtype=float)
        y_edges = np.array([ay.GetBinLowEdge(i + 1) for i in range(ny)] + [ay.GetBinUpEdge(ny)], dtype=float)
        values = np.zeros((nx, ny), dtype=float)
        for ix in range(nx):
            for iy in range(ny):
                values[ix, iy] = h.GetBinContent(ix + 1, iy + 1)
    finally:
        f.Close()
    return values, x_edges, y_edges


def read_th1(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray]:
    try:
        return _read_th1_with_pyroot(root_path, name)
    except (ModuleNotFoundError, ImportError):
        return _read_th1_with_uproot(root_path, name)


def read_th2(root_path: Path, name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        return _read_th2_with_pyroot(root_path, name)
    except (ModuleNotFoundError, ImportError):
        return _read_th2_with_uproot(root_path, name)


def bin_centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def plot_th1(root_path: Path, hname: str, xlabel: str, ylabel: str,
             title: str, out_base: Path, formats: Sequence[str], label: str,
             xlim: Optional[Tuple[float, float]] = None) -> List[Path]:
    values, edges = read_th1(root_path, hname)
    centers = bin_centers(edges)
    widths = np.diff(edges)

    fig, ax = plt.subplots(figsize=(7.2, 5.3))
    ax.step(edges[:-1], values, where="post", linewidth=1.9, color="black")
    ax.fill_between(edges[:-1], values, step="post", alpha=0.18, color="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=12)
    ax.grid(True, alpha=0.22, linewidth=0.7)
    if xlim is not None:
        ax.set_xlim(*xlim)

    integral = float(np.sum(values))
    mean = float(np.average(centers, weights=values)) if integral > 0 else float("nan")
    sigma = float(np.sqrt(np.average((centers - mean) ** 2, weights=values))) if integral > 0 else float("nan")
    stats = f"Entries = {integral:.0f}\nMean = {mean:.3g}\nRMS = {sigma:.3g}"
    ax.text(0.98, 0.88, stats, transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.7", alpha=0.85))
    add_cms_labels(ax, label)
    fig.tight_layout()
    return save_all(fig, out_base, formats)


def plot_th2(root_path: Path, hname: str, xlabel: str, ylabel: str,
             cbar_label: str, title: str, out_base: Path, formats: Sequence[str],
             label: str, cmap_name: str = "cms_white_red",
             vmin: Optional[float] = None, vmax: Optional[float] = None) -> List[Path]:
    values, x_edges, y_edges = read_th2(root_path, hname)

    cmap = cms_red_cmap() if cmap_name == "cms_white_red" else plt.get_cmap(cmap_name)
    fig, ax = plt.subplots(figsize=(7.2, 5.9))
    mesh = ax.pcolormesh(x_edges, y_edges, values.T, shading="auto", cmap=cmap,
                         vmin=vmin, vmax=vmax, rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=12)
    add_cms_labels(ax, label)
    fig.tight_layout()
    return save_all(fig, out_base, formats)


def plot_corr_matrix(csv_path: Path, out_base: Path, formats: Sequence[str], label: str) -> List[Path]:
    import pandas as pd

    df = pd.read_csv(csv_path)
    preferred_cols = [
        "xA", "yA", "x_pred", "y_pred", "xT", "yT",
        "dx_raw", "dy_raw", "dx", "dy", "dt",
        "theta_x_deg", "theta_y_deg", "theta_deg",
    ]
    cols = [c for c in preferred_cols if c in df.columns]
    if len(cols) < 2:
        raise SystemExit(f"Not enough columns for a correlation matrix. Available columns: {list(df.columns)}")

    corr = df[cols].corr(method="pearson")
    out_base.parent.mkdir(parents=True, exist_ok=True)
    corr.to_csv(out_base.with_suffix(".csv"))

    fig, ax = plt.subplots(figsize=(10.5, 9.0))
    im = ax.imshow(corr.values, vmin=-1.0, vmax=1.0, cmap="coolwarm")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pearson correlation coefficient")

    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_yticklabels(cols)
    ax.set_title("Tracklet-level correlation matrix near WP", pad=12)

    for i in range(len(cols)):
        for j in range(len(cols)):
            val = corr.values[i, j]
            color = "white" if abs(val) > 0.55 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)

    add_cms_labels(ax, label)
    fig.tight_layout()
    return save_all(fig, out_base, formats)


def main() -> None:
    ap = argparse.ArgumentParser(description="Make WP validation angle/residual/correlation plots.")
    ap.add_argument("--root", required=True, help="ROOT file closest to WP, e.g. tracklets_efficiency_HV6_beam.root")
    ap.add_argument("--csv", default=None, help="Per-tracklet CSV for the same HV, e.g. tracklets_full_HV6_beam.csv")
    ap.add_argument("--out-dir", default="wp_validation_plots")
    ap.add_argument("--label", default="", help="Text block placed inside the plots.")
    ap.add_argument("--png", action="store_true")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--dpi", type=int, default=250)
    ap.add_argument("--dx-lim", nargs=2, type=float, default=None, metavar=("MIN", "MAX"))
    ap.add_argument("--dy-lim", nargs=2, type=float, default=None, metavar=("MIN", "MAX"))
    ap.add_argument("--theta-x-lim", nargs=2, type=float, default=None, metavar=("MIN", "MAX"))
    ap.add_argument("--theta-y-lim", nargs=2, type=float, default=None, metavar=("MIN", "MAX"))
    args = ap.parse_args()

    formats: List[str] = []
    if args.png:
        formats.append("png")
    if args.pdf:
        formats.append("pdf")
    if not formats:
        formats = ["png"]

    setup_style(args.dpi)
    root_path = Path(args.root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: List[Path] = []
    written += plot_th1(root_path, "h_dx", r"$\Delta x$ [cm]", "Matched tracklets",
                        r"$\Delta x$ distribution", out_dir / "03_delta_x_distribution",
                        formats, args.label, tuple(args.dx_lim) if args.dx_lim else None)
    written += plot_th1(root_path, "h_dy", r"$\Delta y$ [cm]", "Matched tracklets",
                        r"$\Delta y$ distribution", out_dir / "04_delta_y_distribution",
                        formats, args.label, tuple(args.dy_lim) if args.dy_lim else None)
    written += plot_th1(root_path, "h_theta_x", r"$\theta_x$ [deg]", "Matched tracklets",
                        r"$\theta_x$ distribution", out_dir / "05_theta_x_distribution",
                        formats, args.label, tuple(args.theta_x_lim) if args.theta_x_lim else None)
    written += plot_th1(root_path, "h_theta_y", r"$\theta_y$ [deg]", "Matched tracklets",
                        r"$\theta_y$ distribution", out_dir / "06_theta_y_distribution",
                        formats, args.label, tuple(args.theta_y_lim) if args.theta_y_lim else None)
    written += plot_th2(root_path, "h_theta_x_theta_y", r"$\theta_x$ [deg]", r"$\theta_y$ [deg]",
                        "Matched tracklets", r"$\theta_x$ vs $\theta_y$", out_dir / "07_theta_x_vs_theta_y_2d",
                        formats, args.label)

    if args.csv:
        written += plot_corr_matrix(Path(args.csv), out_dir / "08_tracklet_correlation_matrix", formats, args.label)
    else:
        print("[warning] --csv not given: skipped 08_tracklet_correlation_matrix")

    print("Wrote:")
    for p in written:
        print(f"  {p}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        raise SystemExit(130)
