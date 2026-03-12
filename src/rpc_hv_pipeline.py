#!/usr/bin/env python3
"""
rpc_hv_pipeline.py

Run the RPC HV-scan workflow with regime-aware outputs and optional per-HV
dynamic cut estimation.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HV_MAP = {
    "HV1": 6.0, "HV2": 6.5, "HV3": 6.7, "HV4": 6.8, "HV5": 6.9,
    "HV6": 7.0, "HV7": 7.1, "HV8": 7.2, "HV9": 7.3, "HV10": 7.4, "HV11": 7.5,
}

SCRIPT_DIR = Path(__file__).resolve().parent


def resolve_script(path_value, default_name):
    if path_value:
        return str(Path(path_value))
    return str(SCRIPT_DIR / default_name)


def run(cmd, logfile=None):
    print(" ".join(cmd))
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0

    combined = ""
    if p.stdout:
        combined += p.stdout
    if p.stderr:
        if combined and not combined.endswith("\n"):
            combined += "\n"
        combined += p.stderr

    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        Path(logfile).write_text(combined, encoding="utf-8")

    if combined.strip():
        print(combined, end="" if combined.endswith("\n") else "\n")

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


def load_cuts_json(path):
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"ERROR: --cuts-json file not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise SystemExit("ERROR: --cuts-json must contain a JSON object keyed by HV tag.")
    return data


def resolve_hv_cuts(hvtag, args, cuts_db):
    entry = cuts_db.get(hvtag, {}) if isinstance(cuts_db, dict) else {}
    return {
        "tol_x": float(entry.get("tol_x", args.tol_x)),
        "tol_y": float(entry.get("tol_y", args.tol_y)),
        "tol_t": float(entry.get("tol_t", args.tol_t)),
        "theta_max_deg": None if entry.get("theta_max_deg") is None else float(entry.get("theta_max_deg")),
        "source": "json" if entry else "defaults",
    }


def write_metadata(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def maybe_generate_cuts(args, mode, chamber_bases, z_strs, out_dir):
    cuts_script = args.cuts_script
    if args.auto_cuts and not cuts_script:
        cuts_script = str(Path(args.auto_cuts_script))
    if not cuts_script:
        return args.cuts_json

    cuts_out = Path(args.cuts_json) if args.cuts_json else (out_dir / f"cuts_by_hv_{args.regime}.json")
    cuts_out.parent.mkdir(parents=True, exist_ok=True)
    log_out = out_dir / "logs" / f"cuts_{args.regime}.log"

    cmd = [
        args.python, cuts_script,
        "--out", str(cuts_out),
        "--regime", args.regime,
        "--tree", args.tree,
        "--hit-mode", args.hit_mode,
        "--mode", mode,
        "--min-samples", str(args.cuts_min_samples),
        "--dx-nsigma", str(args.cuts_dx_nsigma),
        "--dy-nsigma", str(args.cuts_dy_nsigma),
        "--dt-nsigma", str(args.cuts_dt_nsigma),
        "--min-tol-x", str(args.cuts_min_tol_x),
        "--max-tol-x", str(args.cuts_max_tol_x),
        "--min-tol-y", str(args.cuts_min_tol_y),
        "--max-tol-y", str(args.cuts_max_tol_y),
        "--min-tol-t", str(args.cuts_min_tol_t),
        "--max-tol-t", str(args.cuts_max_tol_t),
        "--percentile", str(args.cuts_percentile),
        "--angle-max-deg", str(args.cuts_angle_max_deg),
    ]

    if args.require_ncluster_eq1:
        cmd += ["--require-ncluster-eq1"]

    if args.chambers:
        cmd += ["--chambers", *[str(p) for p in chamber_bases]]
    else:
        cmd += ["--ref-base", str(chamber_bases[0]), "--test-base", str(chamber_bases[1])]

    if mode == "optionc" and z_strs is not None:
        cmd += ["--z", *z_strs, "--anchor-index", str(args.anchor_index), "--target-index", str(args.target_index)]

    if args.roi_name:
        cmd += ["--roi-name", args.roi_name]
    if args.roi_rect:
        cmd += ["--roi-rect", *[str(x) for x in args.roi_rect]]

    run(cmd, logfile=log_out)
    if not cuts_out.exists():
        raise RuntimeError(f"Cuts script did not produce {cuts_out}. See {log_out}")
    return str(cuts_out)


def main():
    t_all = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-base")
    ap.add_argument("--test-base")
    ap.add_argument("--chambers", nargs="+")
    ap.add_argument("--out-dir", default="out_hv")

    ap.add_argument("--python", default=sys.executable if sys.executable else "python3")
    ap.add_argument("--producer", default=None)
    ap.add_argument("--producer-mode", choices=["auto", "v2", "optionc"], default="auto")
    ap.add_argument("--summarizer", default=None)
    ap.add_argument("--hvscript", default=None)
    ap.add_argument("--cuts-script", help="Optional cut-estimation script.")
    ap.add_argument("--cuts-json", help="Optional JSON with per-HV cuts.")
    ap.add_argument("--auto-cuts", action="store_true", help="Estimate per-HV cuts before running the HV loop.")
    ap.add_argument("--auto-cuts-script", default="rpc_estimate_tracklet_cuts.py",
                    help="Script used when --auto-cuts is enabled and --cuts-script is not provided.")

    ap.add_argument("--tree", default="events")
    ap.add_argument("--roi-name", default="scint_overlap")
    ap.add_argument("--roi-rect", nargs=4, type=float, default=[0, 22, 0, 130], metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    ap.add_argument("--regime", choices=["normal", "beam"], default="normal")

    ap.add_argument("--tol-x", type=float, default=7.0)
    ap.add_argument("--tol-y", type=float, default=20.0)
    ap.add_argument("--tol-t", type=float, default=3.0)
    ap.add_argument("--hit-mode", default="earliest", choices=["earliest", "centroid"])
    ap.add_argument("--require-ncluster-eq1", action="store_true")

    ap.add_argument("--z", nargs="+", type=float)
    ap.add_argument("--z0", type=float, default=0.0)
    ap.add_argument("--z1", type=float, default=100.0)
    ap.add_argument("--anchor-index", type=int, default=0)
    ap.add_argument("--target-index", type=int, default=-1)

    ap.add_argument("--cuts-min-samples", type=int, default=50)
    ap.add_argument("--cuts-dx-nsigma", type=float, default=3.0)
    ap.add_argument("--cuts-dy-nsigma", type=float, default=3.0)
    ap.add_argument("--cuts-dt-nsigma", type=float, default=3.0)
    ap.add_argument("--cuts-min-tol-x", type=float, default=3.0)
    ap.add_argument("--cuts-max-tol-x", type=float, default=15.0)
    ap.add_argument("--cuts-min-tol-y", type=float, default=6.0)
    ap.add_argument("--cuts-max-tol-y", type=float, default=35.0)
    ap.add_argument("--cuts-min-tol-t", type=float, default=1.0)
    ap.add_argument("--cuts-max-tol-t", type=float, default=8.0)
    ap.add_argument("--cuts-percentile", type=float, default=0.95)
    ap.add_argument("--cuts-angle-max-deg", type=float, default=25.0)

    ap.add_argument("--fit", action="store_true")
    args = ap.parse_args()

    args.producer = resolve_script(args.producer, "rpc_tracklet_efficiency.py")
    args.summarizer = resolve_script(args.summarizer, "rpc_make_eff_summary.py")
    args.hvscript = resolve_script(args.hvscript, "rpc_efficiency_vs_hv.py")
    if args.cuts_script:
        args.cuts_script = str(Path(args.cuts_script))
    else:
        args.auto_cuts_script = resolve_script(args.auto_cuts_script, "rpc_estimate_tracklet_cuts.py")

    mode = args.producer_mode if args.producer_mode != "auto" else infer_producer_mode(args.producer)

    out_dir = Path(args.out_dir)
    map_dir = out_dir / "maps"
    sum_dir = out_dir / "summaries"
    plot_dir = out_dir / "plots"
    log_dir = out_dir / "logs"
    meta_dir = out_dir / "metadata"
    for d in (map_dir, sum_dir, plot_dir, log_dir, meta_dir):
        d.mkdir(parents=True, exist_ok=True)

    if args.chambers:
        chamber_bases = [Path(p) for p in args.chambers]
        n_ch = len(chamber_bases)
        if n_ch < 2:
            raise SystemExit("ERROR: --chambers must have at least 2 dirs.")
    else:
        if not args.ref_base or not args.test_base:
            raise SystemExit("ERROR: Provide either --chambers or both --ref-base and --test-base.")
        chamber_bases = [Path(args.ref_base), Path(args.test_base)]
        n_ch = 2

    z_strs = None
    if mode == "optionc":
        if args.z is not None:
            z_strs = parse_z_list(args.z, n_ch)
        else:
            if n_ch == 2:
                z_strs = [str(args.z0), str(args.z1)]
            else:
                raise SystemExit("ERROR: For 3/4-ch optionC runs you must provide --z.")

    cuts_json_path = maybe_generate_cuts(args, mode, chamber_bases, z_strs, out_dir)
    cuts_db = load_cuts_json(cuts_json_path)

    hv_dirs = sorted([d for d in chamber_bases[0].glob("HV*") if d.is_dir()], key=lambda x: (len(x.name), x.name))

    for hvdir in hv_dirs:
        hvtag = hvdir.name
        if hvtag not in HV_MAP:
            print(f"SKIP {hvtag}: not in HV_MAP")
            continue
        hv_value = HV_MAP[hvtag]
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

        hv_cuts = resolve_hv_cuts(hvtag, args, cuts_db)
        print(f"=== {hvtag}  ({hv_value:.2f} kV) [{args.regime}] ===")
        print(f"Using cuts: tol_x={hv_cuts['tol_x']:.3f}, tol_y={hv_cuts['tol_y']:.3f}, tol_t={hv_cuts['tol_t']:.3f} [{hv_cuts['source']}]")
        if hv_cuts["theta_max_deg"] is not None:
            print(f"Optional theta_max_deg={hv_cuts['theta_max_deg']:.3f}")

        out_root_name = f"tracklets_efficiency_{hvtag}.root"
        out_csv = "tracklets_full.csv"
        out_ndjson = "tracklets.ndjson"

        map_out = map_dir / f"tracklets_efficiency_{hvtag}_{args.regime}.root"
        csv_out = map_dir / f"tracklets_full_{hvtag}_{args.regime}.csv"
        ndjson_out = map_dir / f"tracklets_{hvtag}_{args.regime}.ndjson"
        sum_out = sum_dir / f"summary_{hvtag}_{args.regime}.json"
        meta_out = meta_dir / f"metadata_{hvtag}_{args.regime}.json"
        log_out = log_dir / f"producer_{hvtag}_{args.regime}.log"

        for f in (out_root_name, out_csv, out_ndjson):
            Path(f).unlink(missing_ok=True)

        if mode == "v2":
            if n_ch != 2:
                raise SystemExit("ERROR: v2-mode supports only 2 chambers.")
            cmd = [
                args.python, args.producer,
                "--fileA", str(hv_files[0]),
                "--fileB", str(hv_files[1]),
                "--tree", args.tree,
                "--tag", hvtag,
                "--tol_x", str(hv_cuts["tol_x"]),
                "--tol_y", str(hv_cuts["tol_y"]),
                "--tol_t", str(hv_cuts["tol_t"]),
                "--hit_mode", args.hit_mode,
            ]
            if args.require_ncluster_eq1:
                cmd += ["--require_ncluster_eq1"]
        else:
            cmd = [
                args.python, args.producer,
                "--files", *[str(f) for f in hv_files],
                "--z", *z_strs,
                "--anchor-index", str(args.anchor_index),
                "--target-index", str(args.target_index),
                "--tree", args.tree,
                "--tag", hvtag,
                "--tol-x", str(hv_cuts["tol_x"]),
                "--tol-y", str(hv_cuts["tol_y"]),
                "--tol-t", str(hv_cuts["tol_t"]),
                "--hit-mode", args.hit_mode,
                "--regime", args.regime,
            ]
            if args.require_ncluster_eq1:
                cmd += ["--require-ncluster-eq1"]
            if hv_cuts["theta_max_deg"] is not None:
                cmd += ["--theta-max-deg", str(hv_cuts["theta_max_deg"])]

        run(cmd, logfile=log_out)
        if not Path(out_root_name).exists():
            raise RuntimeError(f"Producer did not produce {out_root_name}. See {log_out}")

        Path(out_root_name).replace(map_out)
        if Path(out_csv).exists():
            Path(out_csv).replace(csv_out)
        if Path(out_ndjson).exists():
            Path(out_ndjson).replace(ndjson_out)

        write_metadata(meta_out, {
            "hv_tag": hvtag,
            "hv_value_kv": hv_value,
            "regime": args.regime,
            "producer": args.producer,
            "producer_mode": mode,
            "tree": args.tree,
            "hit_mode": args.hit_mode,
            "require_ncluster_eq1": args.require_ncluster_eq1,
            "tol_x": hv_cuts["tol_x"],
            "tol_y": hv_cuts["tol_y"],
            "tol_t": hv_cuts["tol_t"],
            "theta_max_deg": hv_cuts["theta_max_deg"],
            "roi_name": args.roi_name,
            "roi_rect": list(args.roi_rect),
            "input_files": [str(f) for f in hv_files],
            "z": [float(z) for z in z_strs] if z_strs is not None else None,
            "anchor_index": args.anchor_index if mode == "optionc" else None,
            "target_index": args.target_index if mode == "optionc" else None,
            "cuts_source": hv_cuts["source"],
            "cuts_json": str(cuts_json_path) if cuts_json_path else None,
        })

        cmd_sum = [
            args.python, args.summarizer,
            "--root", str(map_out),
            "--hv", str(hv_value),
            "--roi-name", args.roi_name,
            "--roi-rect", str(args.roi_rect[0]), str(args.roi_rect[1]), str(args.roi_rect[2]), str(args.roi_rect[3]),
            "--out", str(sum_out),
        ]
        run(cmd_sum)

    plot_out = plot_dir / f"eff_{args.roi_name}_vs_hv_{args.regime}.png"
    fit_out = plot_dir / f"sigmoid_fit_{args.roi_name}_{args.regime}.json"
    cmd_hv = [args.python, args.hvscript, "--glob", str(sum_dir / f"summary_HV*_{args.regime}.json"), "--out", str(plot_out)]
    if args.fit:
        cmd_hv += ["--fit", "--fit-out", str(fit_out)]
    run(cmd_hv)

    print(f"ALL DONE in {time.time() - t_all:.2f} s")
    print(f"Regime: {args.regime}")
    if cuts_json_path:
        print(f"Cuts:   {cuts_json_path}")
    print(f"Plot:   {plot_out}")
    if args.fit:
        print(f"Fit:    {fit_out}")


if __name__ == "__main__":
    main()
