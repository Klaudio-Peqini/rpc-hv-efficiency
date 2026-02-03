#!/usr/bin/env python3
"""
make_synthetic_rpc_hvscan.py

Calibration-grade synthetic ROOT generator for RPC efficiency-vs-HV regression testing.

Main goal:
- Generate synthetic data with the SAME branch schema as your real data:
    Strip_X (vector<float>)
    Strip_Y (vector<float>)
    Strip_Tdc_diff (vector<float>)
    nCluster (Int_t scalar)
    Cluster_size (Int_t scalar)

- Generate 2/3/4 chambers with HV scan folder layout:
    OUT_BASE/
      ChamberSim1/HV1/data.root ... HV11/data.root
      ...
      ChamberSimN/HV*/data.root

- Enforce per-HV efficiency points and binomial errors:
    eff = matched/eligible
    err_binom = sqrt(eff*(1-eff)/eligible)

We construct integer (eligible, matched) pairs per HV that match the provided
targets to 4 decimal places, and then build the chamber hits deterministically
so your analysis recovers these values robustly.

Notes:
- Default is n-events = 60000 (matches your typical stats).
- Eligibility is defined via reference chambers having a hit (nCluster==1).
  Events not meant to be eligible have reference chambers set to nCluster=0.
- Target chamber hit decisions are enforced exactly among eligible events.
- Optional realism knobs exist (misalignment, noise-only, multicluster), but
  default to 0 so the calibration dataset is stable.

Requires:
- ROOT with PyROOT enabled (import ROOT).
"""

import argparse
import math
import random
from array import array
from pathlib import Path
from typing import Dict, Tuple, List

import ROOT

# ----------------------------------------------------------------------
# HV mapping (must match your analysis/pipeline)
# ----------------------------------------------------------------------
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
HV_ORDER = list(HV_MAP.keys())

# ----------------------------------------------------------------------
# Target efficiency points (your provided values)
# These are assumed to refer to the TARGET chamber efficiency among eligible events.
# ----------------------------------------------------------------------
TARGET_POINTS = {
    "HV1":  (0.0000, 0.0000),
    "HV2":  (0.0750, 0.0120),
    "HV3":  (0.5230, 0.0100),
    "HV4":  (0.7685, 0.0088),
    "HV5":  (0.8890, 0.0074),
    "HV6":  (0.9603, 0.0060),
    "HV7":  (0.9621, 0.0051),
    "HV8":  (0.9639, 0.0035),
    "HV9":  (0.9644, 0.0022),
    "HV10": (0.9649, 0.0010),
    "HV11": (0.9651, 0.0007),
}

# ----------------------------------------------------------------------
# Binomial error
# ----------------------------------------------------------------------
def binom_err(p: float, n: int) -> float:
    if n <= 0:
        return 0.0
    return math.sqrt(p * (1.0 - p) / n)

# ----------------------------------------------------------------------
# Find integer (eligible, matched) such that rounding matches targets.
# We require:
#   round(m/n, 4) == p_target
#   round(sqrt(p*(1-p)/n), 4) == err_target
# This makes the summary match the exact formatted values reliably.
# ----------------------------------------------------------------------
def find_counts_for_point(p_target: float, err_target: float) -> Tuple[int, int]:
    # Special case: p=0 and err=0
    if abs(p_target) < 1e-12 and abs(err_target) < 1e-12:
        # simplest: n=1, m=0 -> eff=0.0000 err=0.0000
        return 1, 0

    # Initial estimate from n ~ p(1-p)/err^2
    n_est = int(round(p_target * (1.0 - p_target) / (err_target ** 2)))
    n_est = max(n_est, 2)

    # Search around estimate
    # We widen the window for hard cases.
    lo = max(2, n_est - 50000)
    hi = n_est + 50000

    # Prefer smaller n if multiple solutions exist (but still realistic)
    best = None

    for n in range(lo, hi + 1):
        # ideal m
        m_float = p_target * n
        m = int(round(m_float))
        if m < 0: m = 0
        if m > n: m = n

        p = m / n
        e = binom_err(p, n)

        if round(p, 4) == round(p_target, 4) and round(e, 4) == round(err_target, 4):
            best = (n, m)
            break

    if best is None:
        # Fallback: try to match p exactly (rounded) and accept closest error rounding
        # (rare, but safer than failing)
        best_n, best_m = None, None
        best_score = 1e9
        for n in range(lo, hi + 1):
            m = int(round(p_target * n))
            m = min(max(m, 0), n)
            p = m / n
            if round(p, 4) != round(p_target, 4):
                continue
            e = binom_err(p, n)
            score = abs(round(e, 4) - round(err_target, 4))
            if score < best_score:
                best_score = score
                best_n, best_m = n, m
                if score == 0.0:
                    break
        if best_n is None:
            raise RuntimeError(f"Could not find counts for p={p_target}, err={err_target}")
        best = (best_n, best_m)

    return best

# ----------------------------------------------------------------------
# Truth track model: simple straight lines shared across all chambers
# ----------------------------------------------------------------------
def sample_truth_track(rng: random.Random, x0_range, y0_range, slope_sigma):
    x0 = rng.uniform(*x0_range)
    y0 = rng.uniform(*y0_range)
    tx = rng.gauss(0.0, slope_sigma)
    ty = rng.gauss(0.0, slope_sigma)
    return x0, y0, tx, ty

# ----------------------------------------------------------------------
# Create cluster vectors around a centroid (x,y,t) with smear
# ----------------------------------------------------------------------
def fill_cluster_vectors(
    rng: random.Random,
    vX, vY, vT,
    x_cm: float, y_cm: float, t_ns: float,
    cluster_size: int,
    smear_xy: float, smear_t: float
):
    vX.clear(); vY.clear(); vT.clear()
    for _ in range(cluster_size):
        vX.push_back(float(x_cm + rng.gauss(0.0, smear_xy)))
        vY.push_back(float(y_cm + rng.gauss(0.0, smear_xy)))
        vT.push_back(float(t_ns + rng.gauss(0.0, smear_t)))

def clear_vectors(vX, vY, vT):
    vX.clear(); vY.clear(); vT.clear()

# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-base", required=True, help="Output base directory.")
    ap.add_argument("--n-chambers", type=int, choices=[2, 3, 4], default=2)
    ap.add_argument("--n-events", type=int, default=60000, help="Total events per HV file (default 60000).")
    ap.add_argument("--tree", default="events")
    ap.add_argument("--seed", type=int, default=12345)

    # z positions (cm)
    ap.add_argument("--z", nargs="+", type=float, default=None,
                    help="Z positions for chambers (cm). Must have N values if provided.")
    ap.add_argument("--z0", type=float, default=0.0)
    ap.add_argument("--dz", type=float, default=50.0, help="Default spacing if --z not given.")

    # truth phase space (cm)
    ap.add_argument("--x0-min", type=float, default=5.0)
    ap.add_argument("--x0-max", type=float, default=17.0)
    ap.add_argument("--y0-min", type=float, default=10.0)
    ap.add_argument("--y0-max", type=float, default=120.0)
    ap.add_argument("--slope-sigma", type=float, default=0.01)

    # measurement smear
    ap.add_argument("--hit-smear-xy", type=float, default=0.2)
    ap.add_argument("--hit-smear-t", type=float, default=0.4)
    ap.add_argument("--t0-ns", type=float, default=0.0)

    # cluster size model
    ap.add_argument("--cluster-size-mean", type=float, default=2.0)
    ap.add_argument("--cluster-size-sigma", type=float, default=0.8)
    ap.add_argument("--cluster-size-min", type=int, default=1)
    ap.add_argument("--cluster-size-max", type=int, default=6)

    # optional realism knobs (default 0 for calibration stability)
    ap.add_argument("--p-multicluster", type=float, default=0.0,
                    help="Probability to add extra hits (nCluster=2) even when hit exists. Default 0.")
    ap.add_argument("--p-noise-only", type=float, default=0.0,
                    help="Probability of noise-only hit when otherwise miss. Default 0.")

    # misalignment offsets (cm), per chamber
    ap.add_argument("--dx", nargs="+", type=float, default=None)
    ap.add_argument("--dy", nargs="+", type=float, default=None)

    # which chamber is target (default last)
    ap.add_argument("--target-index", type=int, default=-1)

    args = ap.parse_args()
    rng = random.Random(args.seed)

    nC = args.n_chambers
    target_index = args.target_index if args.target_index >= 0 else (nC - 1)
    if not (0 <= target_index < nC):
        raise SystemExit("--target-index out of range")

    # z positions
    if args.z is not None:
        if len(args.z) != nC:
            raise SystemExit(f"--z must have {nC} values")
        zpos = list(args.z)
    else:
        zpos = [args.z0 + i * args.dz for i in range(nC)]

    # misalignment offsets
    dx = list(args.dx) if args.dx is not None else [0.0] * nC
    dy = list(args.dy) if args.dy is not None else [0.0] * nC
    if len(dx) != nC or len(dy) != nC:
        raise SystemExit("--dx/--dy must have N values if provided")

    out_base = Path(args.out_base)
    out_base.mkdir(parents=True, exist_ok=True)

    # Create chamber folders
    chamber_dirs = []
    for i in range(nC):
        chdir = out_base / f"ChamberSim{i+1}"
        chamber_dirs.append(chdir)
        for hvtag in HV_ORDER:
            (chdir / hvtag).mkdir(parents=True, exist_ok=True)

    # Pre-generate truth tracks (shared across all HV & chambers)
    x0_range = (args.x0_min, args.x0_max)
    y0_range = (args.y0_min, args.y0_max)
    tracks = [sample_truth_track(rng, x0_range, y0_range, args.slope_sigma) for _ in range(args.n_events)]

    # Precompute exact (eligible, matched) per HV for the target chamber
    counts: Dict[str, Tuple[int, int]] = {}
    for hvtag in HV_ORDER:
        p_t, e_t = TARGET_POINTS[hvtag]
        n, m = find_counts_for_point(p_t, e_t)
        counts[hvtag] = (n, m)

    print("[INFO] Using the following (eligible, matched) pairs for the TARGET chamber:")
    for hvtag in HV_ORDER:
        n, m = counts[hvtag]
        p = m / n
        e = binom_err(p, n)
        print(f"  {hvtag}: eligible={n:6d}, matched={m:6d}, eff={p:.4f}, err={e:.4f}")

    # Helper to create a ROOT file + branches
    def make_tree(file_path: Path):
        f = ROOT.TFile(str(file_path), "RECREATE")
        t = ROOT.TTree(args.tree, args.tree)

        vX = ROOT.std.vector('float')()
        vY = ROOT.std.vector('float')()
        vT = ROOT.std.vector('float')()

        nCluster = array('i', [0])
        cSize = array('i', [0])

        t.Branch("Strip_X", vX)
        t.Branch("Strip_Y", vY)
        t.Branch("Strip_Tdc_diff", vT)
        t.Branch("nCluster", nCluster, "nCluster/I")
        t.Branch("Cluster_size", cSize, "Cluster_size/I")

        return f, t, vX, vY, vT, nCluster, cSize

    # Generate files per HV
    for hvtag in HV_ORDER:
        hv_kv = HV_MAP[hvtag]
        eligible_target, matched_target = counts[hvtag]

        print(f"\n[GEN] {hvtag} ({hv_kv:.2f} kV)  target eligible={eligible_target} matched={matched_target}")

        # Make ROOT files for all chambers at this HV
        trees = []
        for ci in range(nC):
            file_path = chamber_dirs[ci] / hvtag / "data.root"
            f, t, vX, vY, vT, nCluster, cSize = make_tree(file_path)
            trees.append((f, t, vX, vY, vT, nCluster, cSize))

        # Decide which events are eligible (by reference hits)
        # We will choose exactly eligible_target events to be "eligible"
        # and enforce that ALL reference chambers have a hit there.
        indices = list(range(args.n_events))
        rng.shuffle(indices)
        eligible_set = set(indices[:eligible_target])

        # Among eligible events, choose exactly matched_target for target chamber hits
        eligible_list = list(eligible_set)
        rng.shuffle(eligible_list)
        matched_set = set(eligible_list[:matched_target])

        # Fill events
        for ievt in range(args.n_events):
            x0, y0, tx, ty = tracks[ievt]
            t_evt = args.t0_ns + rng.gauss(0.0, 0.2)

            is_eligible = (ievt in eligible_set)
            is_matched = (ievt in matched_set)

            for ci in range(nC):
                f, t, vX, vY, vT, nCluster, cSize = trees[ci]

                # reset
                clear_vectors(vX, vY, vT)
                nCluster[0] = 0
                cSize[0] = 0

                z = zpos[ci]
                x_true = x0 + tx * z
                y_true = y0 + ty * z

                # Determine chamber "hit"
                if ci == target_index:
                    # Target: only hit if eligible and selected as matched
                    hit = (is_eligible and is_matched)
                    # Allow optional noise-only hits when not matched (stress test)
                    if (not hit) and (args.p_noise_only > 0.0) and (rng.random() < args.p_noise_only):
                        hit = True
                        # put noise away from truth
                        x_true = rng.uniform(args.x0_min - 5.0, args.x0_max + 40.0)
                        y_true = rng.uniform(args.y0_min - 20.0, args.y0_max + 60.0)
                else:
                    # Reference chambers: hit iff eligible (to enforce denominator)
                    hit = is_eligible

                if hit:
                    # Apply misalignment offsets
                    x_meas = x_true + dx[ci]
                    y_meas = y_true + dy[ci]

                    # cluster size
                    cs = int(round(rng.gauss(args.cluster_size_mean, args.cluster_size_sigma)))
                    cs = max(args.cluster_size_min, min(args.cluster_size_max, cs))

                    fill_cluster_vectors(
                        rng, vX, vY, vT,
                        x_meas, y_meas, t_evt,
                        cs, args.hit_smear_xy, args.hit_smear_t
                    )

                    nCluster[0] = 1
                    cSize[0] = cs

                    # Optional multicluster stress (disabled by default)
                    if args.p_multicluster > 0.0 and rng.random() < args.p_multicluster:
                        # Inflate arrays and mark nCluster=2
                        nCluster[0] = 2
                        cs2 = int(round(rng.gauss(args.cluster_size_mean, args.cluster_size_sigma)))
                        cs2 = max(args.cluster_size_min, min(args.cluster_size_max, cs2))
                        # add noise hits far
                        for _ in range(cs2):
                            vX.push_back(float(rng.uniform(args.x0_min - 10.0, args.x0_max + 50.0)))
                            vY.push_back(float(rng.uniform(args.y0_min - 30.0, args.y0_max + 80.0)))
                            vT.push_back(float(t_evt + rng.gauss(5.0, 1.0)))
                        cSize[0] = cs + cs2

                # Fill
                t.Fill()

        # Write & close
        for f, t, *_ in trees:
            f.Write()
            f.Close()

    print("\nDone.")
    print(f"Output base: {out_base.resolve()}")
    print("Example file:")
    print(f"  {out_base}/ChamberSim1/HV6/data.root")

if __name__ == "__main__":
    main()

