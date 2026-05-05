#!/usr/bin/env python3
"""
rpc_plot_efficiency_maps.py

Batch-plot ROOT efficiency maps produced by rpc-hv-efficiency.

Expected input ROOT files contain the standard histograms:
  - h_total     : denominator / eligible events
  - h_tracklet  : numerator / matched reconstructed tracklets
  - h_eff       : local efficiency map, h_tracklet / h_total

The script uses a ROOT/PyROOT reader when available, and falls back to uproot
when PyROOT is not available.  The plotting backend is forced to Agg so the
script can run safely on Wayland/headless machines and inside batch jobs.

Example:
  python3 src/rpc_plot_efficiency_maps.py \
    --input-dir out_dynamic_regime/maps \
    --out-dir out_dynamic_regime/map_plots \
    --png

To reproduce the older visual range shown in some legacy plots:
  python3 src/rpc_plot_efficiency_maps.py \
    --input-dir out_dynamic_regime/maps \
    --plot-range 0 60 0 160
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

import matplotlib

# Avoid Qt/Wayland problems when saving figures from scripts.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


HV_TO_KV = {
    "HV1": 6.0,
    "HV2": 6.5,
    "HV3": 6.7,
    "HV4": 6.8,
    "HV5": 6.9,
    "HV6": 7.0,
    "HV7": 7.1,
    "HV8": 7.2,
    "HV9": 7.3,
    "HV10": 7.4,
    "HV11": 7.5,
}


@dataclass(frozen=True)
class MapData:
    x_edges: np.ndarray
    y_edges: np.ndarray
    total: np.ndarray
    tracklet: np.ndarray
    eff: np.ndarray


@dataclass(frozen=True)
class FileMeta:
    hv_tag: str
    hv_kv: Optional[float]
    hv_label_for_file: str
    hv_label_for_title: str
    regime: str
    run_tag: str


def stamp(t0: float, msg: str) -> None:
    print(f"[{time.perf_counter() - t0:8.2f} s] {msg}")


def format_kv_for_file(kv: float) -> str:
    """Return compact file label: 7kV, 6p8kV, ..."""
    if math.isclose(kv, round(kv), abs_tol=1e-9):
        return f"{int(round(kv))}kV"
    return f"{kv:.1f}".replace(".", "p") + "kV"


def format_kv_for_title(kv: float) -> str:
    """Return plot-title label: 7 kV, 6.8 kV, ..."""
    if math.isclose(kv, round(kv), abs_tol=1e-9):
        return f"{int(round(kv))} kV"
    return f"{kv:.1f} kV"


def parse_file_meta(root_path: Path) -> FileMeta:
    stem = root_path.stem
    run_tag = stem
    prefix = "tracklets_efficiency_"
    if stem.startswith(prefix):
        run_tag = stem[len(prefix) :]

    hv_match = re.search(r"HV\d+", stem, flags=re.IGNORECASE)
    hv_tag = hv_match.group(0).upper() if hv_match else "HV"
    hv_kv = HV_TO_KV.get(hv_tag)

    lower = stem.lower()
    if "beam" in lower:
        regime = "beam"
    elif "normal" in lower:
        regime = "normal"
    else:
        regime = ""

    if hv_kv is None:
        hv_label_for_file = hv_tag
        hv_label_for_title = hv_tag
    else:
        hv_label_for_file = format_kv_for_file(hv_kv)
        hv_label_for_title = format_kv_for_title(hv_kv)

    return FileMeta(
        hv_tag=hv_tag,
        hv_kv=hv_kv,
        hv_label_for_file=hv_label_for_file,
        hv_label_for_title=hv_label_for_title,
        regime=regime,
        run_tag=run_tag,
    )


def _root_th2_to_numpy(hist) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    nx = int(hist.GetNbinsX())
    ny = int(hist.GetNbinsY())
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()

    x_edges = np.array(
        [xaxis.GetBinLowEdge(1 + i) for i in range(nx)]
        + [xaxis.GetBinUpEdge(nx)],
        dtype=float,
    )
    y_edges = np.array(
        [yaxis.GetBinLowEdge(1 + i) for i in range(ny)]
        + [yaxis.GetBinUpEdge(ny)],
        dtype=float,
    )

    values = np.zeros((nx, ny), dtype=float)
    for ix in range(1, nx + 1):
        for iy in range(1, ny + 1):
            values[ix - 1, iy - 1] = float(hist.GetBinContent(ix, iy))
    return values, x_edges, y_edges


def _read_with_pyroot(
    root_path: Path,
    h_total_name: str,
    h_tracklet_name: str,
    h_eff_name: str,
) -> MapData:
    import ROOT  # type: ignore

    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")

    try:
        h_total = f.Get(h_total_name)
        h_tracklet = f.Get(h_tracklet_name)
        h_eff = f.Get(h_eff_name)
        if not h_total:
            raise RuntimeError(f"Histogram '{h_total_name}' not found in {root_path}")
        if not h_tracklet:
            raise RuntimeError(f"Histogram '{h_tracklet_name}' not found in {root_path}")

        total, x_edges, y_edges = _root_th2_to_numpy(h_total)
        tracklet, x_edges_2, y_edges_2 = _root_th2_to_numpy(h_tracklet)

        if not (np.allclose(x_edges, x_edges_2) and np.allclose(y_edges, y_edges_2)):
            raise RuntimeError("h_total and h_tracklet have different binning")

        if h_eff:
            eff, x_edges_3, y_edges_3 = _root_th2_to_numpy(h_eff)
            if not (np.allclose(x_edges, x_edges_3) and np.allclose(y_edges, y_edges_3)):
                raise RuntimeError("h_eff has different binning from h_total")
        else:
            eff = np.divide(
                tracklet,
                total,
                out=np.zeros_like(tracklet, dtype=float),
                where=(total > 0),
            )
    finally:
        f.Close()

    return MapData(x_edges=x_edges, y_edges=y_edges, total=total, tracklet=tracklet, eff=eff)


def _read_with_uproot(
    root_path: Path,
    h_total_name: str,
    h_tracklet_name: str,
    h_eff_name: str,
) -> MapData:
    import uproot  # type: ignore

    with uproot.open(root_path) as f:
        if h_total_name not in f:
            raise RuntimeError(f"Histogram '{h_total_name}' not found in {root_path}")
        if h_tracklet_name not in f:
            raise RuntimeError(f"Histogram '{h_tracklet_name}' not found in {root_path}")

        h_total = f[h_total_name]
        h_tracklet = f[h_tracklet_name]

        total = np.asarray(h_total.values(flow=False), dtype=float)
        tracklet = np.asarray(h_tracklet.values(flow=False), dtype=float)
        x_edges = np.asarray(h_total.axes[0].edges(), dtype=float)
        y_edges = np.asarray(h_total.axes[1].edges(), dtype=float)

        if h_eff_name in f:
            h_eff = f[h_eff_name]
            eff = np.asarray(h_eff.values(flow=False), dtype=float)
        else:
            eff = np.divide(
                tracklet,
                total,
                out=np.zeros_like(tracklet, dtype=float),
                where=(total > 0),
            )

    return MapData(x_edges=x_edges, y_edges=y_edges, total=total, tracklet=tracklet, eff=eff)


def read_maps(
    root_path: Path,
    h_total_name: str,
    h_tracklet_name: str,
    h_eff_name: str,
) -> MapData:
    try:
        return _read_with_pyroot(root_path, h_total_name, h_tracklet_name, h_eff_name)
    except (ModuleNotFoundError, ImportError):
        return _read_with_uproot(root_path, h_total_name, h_tracklet_name, h_eff_name)


def setup_plot_style(dpi: int) -> None:
    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.grid": False,
            "figure.constrained_layout.use": False,
        }
    )


def draw_overlay_polygon(
    ax: plt.Axes,
    coords: Sequence[float],
    color: str,
    linewidth: float,
    alpha: float,
) -> None:
    if len(coords) < 6 or len(coords) % 2 != 0:
        raise ValueError("--overlay-poly must contain an even number of coordinates: x1 y1 x2 y2 ...")
    xy = np.asarray(coords, dtype=float).reshape(-1, 2)
    ax.plot(xy[:, 0], xy[:, 1], color=color, linewidth=linewidth, alpha=alpha)


def plot_one_map(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    out_base: Path,
    title: str,
    cbar_label: str,
    formats: Sequence[str],
    cmap: str,
    vmin: Optional[float],
    vmax: Optional[float],
    plot_range: Optional[Tuple[float, float, float, float]],
    overlay: bool,
    overlay_poly: Sequence[float],
    overlay_color: str,
    overlay_lw: float,
    overlay_alpha: float,
) -> List[Path]:
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    cm = plt.get_cmap(cmap)
    ax.set_facecolor(cm(0.0))

    # uproot/PyROOT convention here is values[ix, iy].  Matplotlib expects
    # pcolormesh data as [iy, ix], therefore the transpose.
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        values.T,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    ax.set_title(title)
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")

    if plot_range is None:
        ax.set_xlim(float(x_edges[0]), float(x_edges[-1]))
        ax.set_ylim(float(y_edges[0]), float(y_edges[-1]))
    else:
        xmin, xmax, ymin, ymax = plot_range
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    if overlay:
        draw_overlay_polygon(ax, overlay_poly, overlay_color, overlay_lw, overlay_alpha)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)

    fig.tight_layout()

    written: List[Path] = []
    for fmt in formats:
        out = out_base.with_suffix(f".{fmt}")
        fig.savefig(out, bbox_inches="tight")
        written.append(out)
    plt.close(fig)
    return written


def discover_root_files(args: argparse.Namespace) -> List[Path]:
    if args.files:
        return sorted(Path(p) for p in args.files)

    input_dir = Path(args.input_dir)
    if args.recursive:
        paths = sorted(input_dir.rglob(args.pattern))
    else:
        paths = sorted(input_dir.glob(args.pattern))
    return [p for p in paths if p.is_file()]


def build_name_tag(meta: FileMeta, legacy_names: bool) -> str:
    if legacy_names:
        # Historical style, e.g. Total_7kV_.pdf.  Avoid this when several
        # regimes are plotted into the same directory, because names collide.
        return f"{meta.hv_label_for_file}_"
    return f"{meta.hv_label_for_file}_{meta.run_tag}"


def output_subdir(base_out_dir: Path, meta: FileMeta, flat: bool) -> Path:
    if flat:
        return base_out_dir
    return base_out_dir / meta.run_tag


def finite_max(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(finite.max()) if finite.size else 1.0


def main() -> None:
    t0 = time.perf_counter()

    parser = argparse.ArgumentParser(
        description="Plot h_total, h_tracklet, and h_eff maps from rpc-hv-efficiency ROOT outputs."
    )
    parser.add_argument("--input-dir", default="out_dynamic_regime/maps", help="Directory containing ROOT map files.")
    parser.add_argument("--pattern", default="tracklets_efficiency_*.root", help="Glob pattern used inside --input-dir.")
    parser.add_argument("--files", nargs="*", default=None, help="Explicit ROOT files. Overrides --input-dir/--pattern.")
    parser.add_argument("--recursive", action="store_true", help="Search --input-dir recursively.")
    parser.add_argument("--out-dir", default=None, help="Output directory. Default: <input-dir>/../map_plots")
    parser.add_argument("--flat", action="store_true", help="Do not create one subfolder per HV/regime tag.")

    parser.add_argument("--h-total", default="h_total", help="Denominator TH2 name.")
    parser.add_argument("--h-tracklet", default="h_tracklet", help="Numerator TH2 name.")
    parser.add_argument("--h-eff", default="h_eff", help="Efficiency TH2 name. Recomputed if missing.")

    parser.add_argument("--png", action="store_true", help="Also write PNG files.")
    parser.add_argument("--no-pdf", action="store_true", help="Do not write PDF files.")
    parser.add_argument("--dpi", type=int, default=200, help="Output DPI for raster formats.")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap.")
    parser.add_argument("--count-vmax", type=float, default=None, help="Fixed vmax for Total/Tracklets maps.")
    parser.add_argument("--eff-vmax", type=float, default=1.0, help="Fixed vmax for Local Efficiency maps.")
    parser.add_argument(
        "--plot-range",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        default=None,
        help="Optional displayed axis range, e.g. --plot-range 0 60 0 160.",
    )

    parser.add_argument("--no-overlay", action="store_true", help="Do not draw the white active-area polygon.")
    parser.add_argument(
        "--overlay-poly",
        nargs="+",
        type=float,
        default=[0.0, 0.0, 0.0, 120.0, 30.0, 120.0, 55.0, 0.0, 0.0, 0.0],
        help="Active-area polygon as x1 y1 x2 y2 ... Default matches the legacy enclosure style.",
    )
    parser.add_argument("--overlay-color", default="white")
    parser.add_argument("--overlay-lw", type=float, default=1.2)
    parser.add_argument("--overlay-alpha", type=float, default=0.85)
    parser.add_argument("--legacy-names", action="store_true", help="Use old names such as Total_7kV_.pdf.")

    args = parser.parse_args()

    formats: List[str] = []
    if not args.no_pdf:
        formats.append("pdf")
    if args.png:
        formats.append("png")
    if not formats:
        raise SystemExit("Nothing to write: remove --no-pdf or add --png.")

    root_files = discover_root_files(args)
    if not root_files:
        raise SystemExit(
            f"No ROOT files found. Checked input-dir={args.input_dir!r}, pattern={args.pattern!r}."
        )

    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir) if args.out_dir else input_dir.parent / "map_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_plot_style(args.dpi)
    stamp(t0, f"Found {len(root_files)} ROOT file(s)")
    stamp(t0, f"Output directory: {out_dir}")

    plot_range = tuple(args.plot_range) if args.plot_range is not None else None
    overlay = not args.no_overlay
    total_written = 0

    for root_path in root_files:
        meta = parse_file_meta(root_path)
        stamp(t0, f"Reading {root_path}")
        maps = read_maps(root_path, args.h_total, args.h_tracklet, args.h_eff)

        target_dir = output_subdir(out_dir, meta, args.flat)
        target_dir.mkdir(parents=True, exist_ok=True)
        name_tag = build_name_tag(meta, args.legacy_names)

        regime_note = f" ({meta.regime})" if meta.regime else ""
        hv_note = meta.hv_label_for_title

        count_vmax = args.count_vmax
        if count_vmax is None:
            # Use a common count scale for denominator and numerator within
            # the same ROOT file, making them directly comparable.
            count_vmax = max(finite_max(maps.total), finite_max(maps.tracklet))

        map_specs = [
            (
                maps.total,
                "Total hits",
                f"Total_{name_tag}",
                target_dir / f"Total_{name_tag}",
                0.0,
                count_vmax,
            ),
            (
                maps.tracklet,
                "Reconstructed tracklets",
                f"Tracklets_{name_tag}",
                target_dir / f"Tracklets_{name_tag}",
                0.0,
                count_vmax,
            ),
            (
                maps.eff,
                "Local Efficiency",
                f"Local_Efficiency_{name_tag}",
                target_dir / f"Local_Efficiency_{name_tag}",
                0.0,
                args.eff_vmax,
            ),
        ]

        for values, title, cbar_label, out_base, vmin, vmax in map_specs:
            title_full = title
            # Keep the legacy title clean, but still record regime/HV in the
            # file and colorbar labels.  Uncomment the next line if desired.
            # title_full = f"{title} - {hv_note}{regime_note}"
            written = plot_one_map(
                values=values,
                x_edges=maps.x_edges,
                y_edges=maps.y_edges,
                out_base=out_base,
                title=title_full,
                cbar_label=cbar_label,
                formats=formats,
                cmap=args.cmap,
                vmin=vmin,
                vmax=vmax,
                plot_range=plot_range,
                overlay=overlay,
                overlay_poly=args.overlay_poly,
                overlay_color=args.overlay_color,
                overlay_lw=args.overlay_lw,
                overlay_alpha=args.overlay_alpha,
            )
            total_written += len(written)
            for out in written:
                stamp(t0, f"Wrote {out}")

        eligible = float(np.sum(maps.total))
        matched = float(np.sum(maps.tracklet))
        eff_global = matched / eligible if eligible > 0 else 0.0
        stamp(
            t0,
            f"{root_path.name}: HV={hv_note}, regime={meta.regime or 'n/a'}, "
            f"matched={matched:.0f}, eligible={eligible:.0f}, global eff={eff_global:.4f}",
        )

    stamp(t0, f"Done. Wrote {total_written} file(s).")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        raise SystemExit(130)
