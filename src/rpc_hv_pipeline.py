#!/usr/bin/env python3
"""
hv_pipeline.py (N-ch extended)

Runs the RPC HV-scan pipeline over folders HV1..HV11 where each HV folder contains data.root.

Supports two producer CLIs:
  (A) v2-style producer:
      tracklets_reconstruction_efficiency_v2.py
      args: --fileA --fileB --tree --tag --tol_x --tol_y --tol_t --hit_mode ...

  (B) optionC N-ch producer:
      tracklets_reconstruction_efficiency_nch_optionC.py
      args: --files f1 f2 ... --z z1 z2 ... --anchor-index --target-index --tree --tag
            --tol-x --tol-y --tol-t --hit-mode ...

This script supports:
  - 2-ch case: provide --ref-base and --test-base (legacy)
  - 3/4-ch case: provide --chambers BASE1 BASE2 BASE3 [BASE4]
      where each BASEk contains HV*/data.root

Downstream scripts (unchanged):
  - make_eff_summary.py
  - efficiency_vs_hv.py

HV mapping is embedded below (edit if needed).
"""

import argparse
import subprocess
import time
from pathlib import Path

# ---------------- HV mapping ----------------
HV_MAP = {
    "HV1":  6.0,
    "HV2":  6.5,
    "HV3":  6.7,
    "HV4":  6.8,
    "HV5":  6.9,
    "HV6":  7.0,
    "HV7":  7.1,
    "HV8":  7.2,
    "HV9":  7.3,
    "HV10": 7.4,
    "HV11": 7.5,
}

def run(cmd, logfile=None):
    print(" ".join(cmd))
    t0 = time.time()
    if logfile:
        with open(logfile, "w") as lf:
            p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    else:
        p = subprocess.run(cmd)
    dt = time.time() - t0
    if p.returncode != 0:
        raise RuntimeError(f"Command failed after {dt:.2f}s: {' '.join(cmd)}")
    print(f"  -> OK ({dt:.2f}s)")

def infer_producer_mode(script_path: str) -> str:
    name = Path(script_path).name.lower()
    if "optionc" in name or "nch" in name:
        return "optionc"
    return "v2"

def parse_z_list(z_list, n):
    if z_list is None:
        return None
    if len(z_list) != n:
        raise SystemExit(f"ERROR: --z must have exactly {n} numbers (got {len(z_list)})")
    return [str(x) for x in z_list]

def main():
    t_all = time.time()
    ap = argparse.ArgumentParser()

    # Inputs: either 2-ch legacy (ref/test) or N-ch (chambers list)
    ap.add_argument("--ref-base", help="Reference chamber base dir (2-ch legacy).")
    ap.add_argument("--test-base", help="Test chamber base dir (2-ch legacy).")
    ap.add_argument("--chambers", nargs="+",
                    help="Base directories for N chambers (3/4...). Each must contain HV*/data.root")

    ap.add_argument("--out-dir", default="out_hv")

    ap.add_argument("--python", default="python3")
    ap.add_argument("--producer", default="src/rpc_tracklets_efficiency.py",
                    help="Producer script (v2 or optionC).")
    ap.add_argument("--producer-mode", choices=["auto", "v2", "optionc"], default="auto",
                    help="How to call the producer script.")
    ap.add_argument("--summarizer", default="src/rpc_make_summary.py")
    ap.add_argument("--hvscript", default="src/rpc_plot_eff_vs_hv_cms.py")

    ap.add_argument("--tree", default="events")
    ap.add_argument("--roi-name", default="scint_overlap")
    ap.add_argument("--roi-rect", nargs=4, type=float, default=[0, 22, 0, 130],
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX"))

    # Matching configuration (mapped to each producer CLI)
    ap.add_argument("--tol-x", type=float, default=7.0)
    ap.add_argument("--tol-y", type=float, default=16.0)
    ap.add_argument("--tol-t", type=float, default=3.0)
    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])
    ap.add_argument("--require-ncluster-eq1", action="store_true",
                    help="Pass nCluster==1 requirement to producer when supported.")

    # OptionC-specific geometry (for N-ch producer)
    ap.add_argument("--z", nargs="+", type=float,
                    help="z positions for chambers (same count/order as --chambers or inferred 2-ch).")
    ap.add_argument("--z0", type=float, default=0.0, help="2-ch fallback z0 (if --z not given)")
    ap.add_argument("--z1", type=float, default=100.0, help="2-ch fallback z1 (if --z not given)")
    ap.add_argument("--anchor-index", type=int, default=0, help="optionC anchor-index (default 0)")
    ap.add_argument("--target-index", type=int, default=-1, help="optionC target-index (default -1=last)")

    ap.add_argument("--fit", action="store_true")

    args = ap.parse_args()

    mode = args.producer_mode
    if mode == "auto":
        mode = infer_producer_mode(args.producer)

    out_dir = Path(args.out_dir)
    map_dir = out_dir / "maps"
    sum_dir = out_dir / "summaries"
    plot_dir = out_dir / "plots"
    log_dir = out_dir / "logs"
    for d in (map_dir, sum_dir, plot_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Determine chamber bases
    chamber_bases = None
    if args.chambers:
        chamber_bases = [Path(p) for p in args.chambers]
        n_ch = len(chamber_bases)
        if n_ch < 2:
            raise SystemExit("ERROR: --chambers must have at least 2 dirs.")
    else:
        # 2-ch legacy
        if not args.ref_base or not args.test_base:
            raise SystemExit("ERROR: Provide either --chambers (N-ch) or both --ref-base and --test-base (2-ch).")
        chamber_bases = [Path(args.ref_base), Path(args.test_base)]
        n_ch = 2

    # Resolve z positions if optionC
    z_strs = None
    if mode == "optionc":
        if args.z is not None:
            z_strs = parse_z_list(args.z, n_ch)
        else:
            if n_ch == 2:
                z_strs = [str(args.z0), str(args.z1)]
            else:
                raise SystemExit("ERROR: For 3/4-ch optionC runs you must provide --z with one value per chamber.")

    # HV directories are taken from the FIRST chamber base (assumed complete set)
    hv_dirs = sorted([d for d in chamber_bases[0].glob("HV*") if d.is_dir()],
                     key=lambda x: (len(x.name), x.name))

    for hvdir in hv_dirs:
        hvtag = hvdir.name
        if hvtag not in HV_MAP:
            print(f"SKIP {hvtag}: not in HV_MAP")
            continue

        hv_value = HV_MAP[hvtag]
        print(f"\n=== {hvtag}  ({hv_value:.2f} kV) ===")

        # Build per-chamber file list for this HV
        hv_files = []
        missing = False
        for base in chamber_bases:
            f = base / hvtag / "data.root"
            if not f.exists():
                print(f"SKIP {hvtag}: missing {f}")
                missing = True
                break
            hv_files.append(f)
        if missing:
            continue

        # Producer outputs in CWD
        out_root_name = f"tracklets_efficiency_{hvtag}.root"
        out_csv = "tracklets_full.csv"
        out_ndjson = "tracklets.ndjson"

        map_out = map_dir / out_root_name
        csv_out = map_dir / f"tracklets_full_{hvtag}.csv"
        ndjson_out = map_dir / f"tracklets_{hvtag}.ndjson"
        sum_out = sum_dir / f"summary_{hvtag}.json"
        log_out = log_dir / f"producer_{hvtag}.log"

        # cleanup leftovers in CWD
        for f in (out_root_name, out_csv, out_ndjson):
            Path(f).unlink(missing_ok=True)

        # --- Step 1: producer ---
        if mode == "v2":
            if n_ch != 2:
                raise SystemExit("ERROR: v2-mode supports only 2 chambers. Use optionC producer for 3/4 chambers.")
            cmd = [
                args.python, args.producer,
                "--fileA", str(hv_files[0]),
                "--fileB", str(hv_files[1]),
                "--tree", args.tree,
                "--tag", hvtag,
                "--tol_x", str(args.tol_x),
                "--tol_y", str(args.tol_y),
                "--tol_t", str(args.tol_t),
                "--hit_mode", args.hit_mode,
            ]
            if args.require_ncluster_eq1:
                cmd += ["--require_ncluster_eq1"]

        else:  # optionC
            cmd = [
                args.python, args.producer,
                "--files", *[str(f) for f in hv_files],
                "--z", *z_strs,
                "--anchor-index", str(args.anchor_index),
                "--target-index", str(args.target_index),
                "--tree", args.tree,
                "--tag", hvtag,
                "--tol-x", str(args.tol_x),
                "--tol-y", str(args.tol_y),
                "--tol-t", str(args.tol_t),
                "--hit-mode", args.hit_mode,
            ]
            if args.require_ncluster_eq1:
                cmd += ["--require-ncluster-eq1"]

        run(cmd, logfile=log_out)

        if not Path(out_root_name).exists():
            raise RuntimeError(f"Producer did not produce {out_root_name}. See {log_out}")

        Path(out_root_name).replace(map_out)
        if Path(out_csv).exists():
            Path(out_csv).replace(csv_out)
        if Path(out_ndjson).exists():
            Path(out_ndjson).replace(ndjson_out)

        # --- Step 2: summarizer ---
        cmd_sum = [
            args.python, args.summarizer,
            "--root", str(map_out),
            "--hv", str(hv_value),
            "--roi-name", args.roi_name,
            "--roi-rect",
            str(args.roi_rect[0]), str(args.roi_rect[1]),
            str(args.roi_rect[2]), str(args.roi_rect[3]),
            "--out", str(sum_out),
        ]
        run(cmd_sum)

    # --- Step 3: HV analysis ---
    plot_out = plot_dir / f"eff_{args.roi_name}_vs_hv.png"
    fit_out = plot_dir / f"sigmoid_fit_{args.roi_name}.json"

    cmd_hv = [
        args.python, args.hvscript,
        "--glob", str(sum_dir / "summary_HV*.json"),
        "--out", str(plot_out),
    ]
    if args.fit:
        cmd_hv += ["--fit", "--fit-out", str(fit_out)]

    run(cmd_hv)

    print(f"\nALL DONE in {time.time() - t_all:.2f} s")
    print(f"Plot: {plot_out}")
    if args.fit:
        print(f"Fit : {fit_out}")

if __name__ == "__main__":
    main()
