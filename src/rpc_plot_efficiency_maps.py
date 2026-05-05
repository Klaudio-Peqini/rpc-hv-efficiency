#!/usr/bin/env python3
"""
rpc_plot_efficiency_maps_cms_style.py

Batch-plot ROOT efficiency maps produced by rpc-hv-efficiency, using a
CMS/GIF++-style visual convention similar to the conventional iRPC hit-profile
plots:

  - CMS Preliminary label at the top-left
  - GIF++ facility label at the top-right
  - inward ticks on all four sides, with minor ticks
  - ROOT/CMS-like white-to-red palette by default
  - active-area chamber outline overlaid in black by default
  - compact detector/run annotation inside the plotting frame

Expected input ROOT files contain the standard histograms:
  - h_total     : denominator / eligible events
  - h_tracklet  : numerator / matched reconstructed tracklets
  - h_eff       : local efficiency map, h_tracklet / h_total

The script uses PyROOT when available and falls back to uproot when PyROOT is
not available.  The plotting backend is forced to Agg so the script can run on
Wayland/headless machines and inside batch jobs.

Examples
--------
Default CMS/GIF++ style:

  python3 rpc_plot_efficiency_maps_cms_style.py \
    --input-dir out_dynamic_regime/maps \
    --png

Force the displayed active chamber range:

  python3 rpc_plot_efficiency_maps_cms_style.py \
    --input-dir out_dynamic_regime/maps \
    --plot-range 0 60 0 160 \
    --png

Use the older simple style with centered titles:

  python3 rpc_plot_efficiency_maps_cms_style.py \
    --input-dir out_dynamic_regime/maps \
    --plain-style
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

import matplotlib

# Avoid Qt/Wayland problems when saving figures from scripts.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import AutoMinorLocator  # noqa: E402


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

DEFAULT_CMS_ANNOTATION = (
    "GIF++ Test Beam November 2025\n"
    "1.4 mm double gap iRPC\n"
    "FEB v2.3 Petiroc 2C\n"
    "threshold ~ 40 fC\n\n"
    "GIF++ source off"
)

# Default outline chosen to match the conventional 60 cm x 160 cm display.
# If you need the older 0--120 cm outline, pass:
#   --overlay-poly 0 0 0 120 30 120 55 0 0 0
DEFAULT_CMS_OVERLAY_POLY = [
    0.0,
    0.0,
    0.0,
    146.0,
    25.0,
    146.0,
    30.0,
    140.0,
    54.0,
    5.0,
    51.0,
    0.0,
    0.0,
    0.0,
]


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


def setup_plot_style(dpi: int, cms_style: bool) -> None:
    base = {
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.grid": False,
        "figure.constrained_layout.use": False,
        "axes.linewidth": 0.9,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }

    if cms_style:
        base.update(
            {
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "xtick.minor.visible": True,
                "ytick.minor.visible": True,
                "xtick.major.size": 6,
                "ytick.major.size": 6,
                "xtick.minor.size": 3,
                "ytick.minor.size": 3,
                "xtick.major.width": 0.8,
                "ytick.major.width": 0.8,
                "xtick.minor.width": 0.7,
                "ytick.minor.width": 0.7,
            }
        )

    plt.rcParams.update(base)


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
    ax.plot(
        xy[:, 0],
        xy[:, 1],
        color=color,
        linewidth=linewidth,
        alpha=alpha,
        solid_joinstyle="miter",
        solid_capstyle="butt",
        zorder=5,
    )


def make_cmap(name: str, cms_style: bool):
    cmap = plt.get_cmap(name).copy()
    if cms_style:
        cmap.set_bad("white")
        cmap.set_under("white")
    return cmap


def maybe_mask_zeros(values: np.ndarray, mask_zeros: bool) -> np.ndarray:
    if not mask_zeros:
        return values
    return np.ma.masked_where(values <= 0.0, values)


def apply_cms_axes_style(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)


def draw_cms_header(
    ax: plt.Axes,
    cms_label: str,
    cms_extra: str,
    facility_label: str,
    fontsize_cms: float,
    fontsize_facility: float,
) -> None:
    left_text = rf"$\bf{{{cms_label}}}$"
    if cms_extra:
        left_text += rf" $\it{{{cms_extra}}}$"
    ax.text(
        0.0,
        1.012,
        left_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize_cms,
        clip_on=False,
    )
    if facility_label:
        ax.text(
            1.0,
            1.012,
            facility_label,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=fontsize_facility,
            clip_on=False,
        )


def draw_cms_annotation(
    ax: plt.Axes,
    text: str,
    fontsize: float,
    x: float,
    y: float,
) -> None:
    if not text.strip():
        return
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fontsize,
        linespacing=1.05,
        clip_on=False,
    )


def finalize_colorbar(cbar, cms_style: bool) -> None:
    if cms_style:
        cbar.ax.tick_params(direction="in", which="both", length=5, width=0.8)
        cbar.ax.minorticks_on()
        cbar.outline.set_linewidth(0.9)


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
    cms_style: bool,
    mask_zeros: bool,
    cms_label: str,
    cms_extra: str,
    facility_label: str,
    cms_annotation: str,
    cms_annotation_fontsize: float,
    cms_annotation_x: float,
    cms_annotation_y: float,
    cms_fontsize: float,
    facility_fontsize: float,
    show_plain_title_in_cms_style: bool,
) -> List[Path]:
    fig, ax = plt.subplots(figsize=(6.7, 6.1) if cms_style else (6.0, 5.2))
    cm = make_cmap(cmap, cms_style=cms_style)
    ax.set_facecolor("white" if cms_style else cm(0.0))

    plot_values = maybe_mask_zeros(values, mask_zeros=mask_zeros)

    # uproot/PyROOT convention here is values[ix, iy]. Matplotlib expects
    # pcolormesh data as [iy, ix], therefore the transpose.
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        plot_values.T,
        shading="auto",
        cmap=cm,
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    if cms_style:
        if show_plain_title_in_cms_style:
            ax.set_title(title, pad=16)
        draw_cms_header(
            ax=ax,
            cms_label=cms_label,
            cms_extra=cms_extra,
            facility_label=facility_label,
            fontsize_cms=cms_fontsize,
            fontsize_facility=facility_fontsize,
        )
        draw_cms_annotation(
            ax=ax,
            text=cms_annotation,
            fontsize=cms_annotation_fontsize,
            x=cms_annotation_x,
            y=cms_annotation_y,
        )
        apply_cms_axes_style(ax)
    else:
        ax.set_title(title)

    ax.set_xlabel("X [cm]", loc="right")
    ax.set_ylabel("Y [cm]", loc="top")

    if plot_range is None:
        ax.set_xlim(float(x_edges[0]), float(x_edges[-1]))
        ax.set_ylim(float(y_edges[0]), float(y_edges[-1]))
    else:
        xmin, xmax, ymin, ymax = plot_range
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    if overlay:
        draw_overlay_polygon(ax, overlay_poly, overlay_color, overlay_lw, overlay_alpha)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.055, fraction=0.055)
    cbar.set_label(cbar_label)
    finalize_colorbar(cbar, cms_style=cms_style)

    if cms_style:
        # Manual margins keep the CMS header and right colorbar label from being clipped.
        fig.subplots_adjust(left=0.115, right=0.86, bottom=0.105, top=0.90)
    else:
        fig.tight_layout()

    written: List[Path] = []
    for fmt in formats:
        out = out_base.with_suffix(f".{fmt}")
        fig.savefig(out, bbox_inches="tight", pad_inches=0.035)
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
        # Historical style, e.g. Total_7kV_.pdf. Avoid this when several
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


def bin_area_grid(x_edges: np.ndarray, y_edges: np.ndarray) -> np.ndarray:
    dx = np.diff(x_edges)
    dy = np.diff(y_edges)
    return dx[:, None] * dy[None, :]


def maybe_density(values: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray, enabled: bool) -> np.ndarray:
    if not enabled:
        return values
    area = bin_area_grid(x_edges, y_edges)
    return np.divide(values, area, out=np.zeros_like(values, dtype=float), where=(area > 0))


def make_parser() -> argparse.ArgumentParser:
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
    parser.add_argument("--dpi", type=int, default=250, help="Output DPI for raster formats.")

    parser.add_argument(
        "--cmap",
        default=None,
        help="Colormap used for all maps. If omitted, --count-cmap and --eff-cmap are used.",
    )
    parser.add_argument("--count-cmap", default="OrRd", help="Colormap for Total/Tracklets maps.")
    parser.add_argument("--eff-cmap", default="OrRd", help="Colormap for Local Efficiency maps.")
    parser.add_argument("--count-vmax", type=float, default=None, help="Fixed vmax for Total/Tracklets maps.")
    parser.add_argument("--eff-vmax", type=float, default=1.0, help="Fixed vmax for Local Efficiency maps.")
    parser.add_argument(
        "--density",
        action="store_true",
        help="Divide count maps by bin area. Use this if the colorbar should literally be [1/cm^2].",
    )
    parser.add_argument(
        "--count-cbar-label",
        default=None,
        help="Colorbar label for count maps. Default: Hits [1/cm^2] with --density, otherwise Hits.",
    )
    parser.add_argument(
        "--tracklet-cbar-label",
        default=None,
        help="Colorbar label for tracklet maps. Default: Tracklets [1/cm^2] with --density, otherwise Tracklets.",
    )
    parser.add_argument("--eff-cbar-label", default="Efficiency", help="Colorbar label for local efficiency map.")
    parser.add_argument(
        "--plot-range",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        default=None,
        help="Optional displayed axis range, e.g. --plot-range 0 60 0 160.",
    )

    parser.add_argument("--plain-style", action="store_true", help="Use the old simple style with centered titles.")
    parser.add_argument("--show-title", action="store_true", help="Also show the map title in CMS style.")
    parser.add_argument("--cms-label", default="CMS", help="Left header label.")
    parser.add_argument("--cms-extra", default="Preliminary", help="Italic text after the CMS label.")
    parser.add_argument("--facility-label", default="GIF++", help="Top-right facility label.")
    parser.add_argument("--cms-fontsize", type=float, default=17.0, help="CMS header font size.")
    parser.add_argument("--facility-fontsize", type=float, default=13.0, help="Facility header font size.")
    parser.add_argument(
        "--cms-annotation",
        default=DEFAULT_CMS_ANNOTATION,
        help="Multiline annotation inside the plot. Use quoted text with '\\n' for line breaks.",
    )
    parser.add_argument("--no-cms-annotation", action="store_true", help="Hide the inner detector/run annotation.")
    parser.add_argument("--cms-annotation-fontsize", type=float, default=10.5)
    parser.add_argument("--cms-annotation-x", type=float, default=0.965)
    parser.add_argument("--cms-annotation-y", type=float, default=0.955)

    parser.add_argument("--no-overlay", action="store_true", help="Do not draw the active-area polygon.")
    parser.add_argument(
        "--overlay-poly",
        nargs="+",
        type=float,
        default=DEFAULT_CMS_OVERLAY_POLY,
        help="Active-area polygon as x1 y1 x2 y2 ... Default matches the CMS/GIF++ enclosure style.",
    )
    parser.add_argument("--overlay-color", default="black", help="Active-area polygon color.")
    parser.add_argument("--overlay-lw", type=float, default=0.9, help="Active-area polygon line width.")
    parser.add_argument("--overlay-alpha", type=float, default=1.0, help="Active-area polygon alpha.")
    parser.add_argument(
        "--mask-zero-counts",
        action="store_true",
        help="Mask zero-valued Total/Tracklets bins so they are pure white.",
    )
    parser.add_argument(
        "--mask-zero-eff",
        action="store_true",
        help="Mask zero-valued efficiency bins so they are pure white.",
    )
    parser.add_argument("--legacy-names", action="store_true", help="Use old names such as Total_7kV_.pdf.")
    return parser


def main() -> None:
    t0 = time.perf_counter()
    parser = make_parser()
    args = parser.parse_args()

    cms_style = not args.plain_style

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

    setup_plot_style(args.dpi, cms_style=cms_style)
    stamp(t0, f"Found {len(root_files)} ROOT file(s)")
    stamp(t0, f"Output directory: {out_dir}")

    plot_range = tuple(args.plot_range) if args.plot_range is not None else None
    overlay = not args.no_overlay
    total_written = 0

    count_cmap = args.cmap if args.cmap is not None else args.count_cmap
    eff_cmap = args.cmap if args.cmap is not None else args.eff_cmap

    count_cbar_label = args.count_cbar_label
    tracklet_cbar_label = args.tracklet_cbar_label
    if count_cbar_label is None:
        count_cbar_label = r"Hits [1/cm$^{2}$]" if args.density else "Hits"
    if tracklet_cbar_label is None:
        tracklet_cbar_label = r"Tracklets [1/cm$^{2}$]" if args.density else "Tracklets"

    cms_annotation = "" if args.no_cms_annotation else args.cms_annotation.replace("\\n", "\n")

    for root_path in root_files:
        meta = parse_file_meta(root_path)
        stamp(t0, f"Reading {root_path}")
        maps = read_maps(root_path, args.h_total, args.h_tracklet, args.h_eff)

        total_values = maybe_density(maps.total, maps.x_edges, maps.y_edges, enabled=args.density)
        tracklet_values = maybe_density(maps.tracklet, maps.x_edges, maps.y_edges, enabled=args.density)
        eff_values = maps.eff

        target_dir = output_subdir(out_dir, meta, args.flat)
        target_dir.mkdir(parents=True, exist_ok=True)
        name_tag = build_name_tag(meta, args.legacy_names)

        hv_note = meta.hv_label_for_title

        count_vmax = args.count_vmax
        if count_vmax is None:
            # Use a common count scale for denominator and numerator within
            # the same ROOT file, making them directly comparable.
            count_vmax = max(finite_max(total_values), finite_max(tracklet_values))

        map_specs = [
            (
                total_values,
                "Total hits",
                count_cbar_label,
                target_dir / f"Total_{name_tag}",
                0.0,
                count_vmax,
                count_cmap,
                args.mask_zero_counts,
            ),
            (
                tracklet_values,
                "Reconstructed tracklets",
                tracklet_cbar_label,
                target_dir / f"Tracklets_{name_tag}",
                0.0,
                count_vmax,
                count_cmap,
                args.mask_zero_counts,
            ),
            (
                eff_values,
                "Local Efficiency",
                args.eff_cbar_label,
                target_dir / f"Local_Efficiency_{name_tag}",
                0.0,
                args.eff_vmax,
                eff_cmap,
                args.mask_zero_eff,
            ),
        ]

        for values, title, cbar_label, out_base, vmin, vmax, cmap, mask_zeros in map_specs:
            written = plot_one_map(
                values=values,
                x_edges=maps.x_edges,
                y_edges=maps.y_edges,
                out_base=out_base,
                title=title,
                cbar_label=cbar_label,
                formats=formats,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                plot_range=plot_range,
                overlay=overlay,
                overlay_poly=args.overlay_poly,
                overlay_color=args.overlay_color,
                overlay_lw=args.overlay_lw,
                overlay_alpha=args.overlay_alpha,
                cms_style=cms_style,
                mask_zeros=mask_zeros,
                cms_label=args.cms_label,
                cms_extra=args.cms_extra,
                facility_label=args.facility_label,
                cms_annotation=cms_annotation,
                cms_annotation_fontsize=args.cms_annotation_fontsize,
                cms_annotation_x=args.cms_annotation_x,
                cms_annotation_y=args.cms_annotation_y,
                cms_fontsize=args.cms_fontsize,
                facility_fontsize=args.facility_fontsize,
                show_plain_title_in_cms_style=args.show_title,
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
