#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_rpc_eff_hv.sh
#   Runs (1) HV-scan analysis and (2) CMS-style plotting
#   for 2, 3, or 4 chambers.
#
# Assumes repo layout:
#   src/rpc_hv_pipeline.py
#   src/rpc_tracklets_efficiency.py
#   src/rpc_plot_eff_vs_hv_cms.py
#
# Input chambers must each contain:
#   HV1/data.root ... HV11/data.root
# ============================================================

# --------- User configuration (edit these) ---------
NCH=2  # 2, 3, or 4

# Provide chamber base dirs (space-separated)
CHAMBERS_2="/path/to/ChA /path/to/ChB"
CHAMBERS_3="/path/to/ChA /path/to/ChB /path/to/ChC"
CHAMBERS_4="/path/to/ChA /path/to/ChB /path/to/ChC /path/to/ChD"

# Z positions in cm (must match order of chambers)
Z_2="0 10"
Z_3="0 50 100"
Z_4="0 40 80 120"

# Indices: anchor reference and target (0-based)
ANCHOR_INDEX=0
TARGET_INDEX_2=1
TARGET_INDEX_3=2
TARGET_INDEX_4=3

# Matching tolerances
TOL_X=7.0
TOL_Y=20.0
TOL_T=3.0

HIT_MODE="centroid"

# Output directory
OUT_DIR="out_hv_${NCH}ch"

# ROI plot (recommended); set to 0 for global plot
USE_ROI=1

# ---------------------------------------------------
PY=python3
PIPELINE="src/rpc_hv_pipeline.py"
PRODUCER="src/rpc_tracklets_efficiency.py"
PLOTTER="src/rpc_plot_eff_vs_hv_cms.py"

# --------- Select chamber list + z + target index ---------
if [[ "${NCH}" == "2" ]]; then
  CHAMBERS="${CHAMBERS_2}"
  Z="${Z_2}"
  TARGET_INDEX="${TARGET_INDEX_2}"
elif [[ "${NCH}" == "3" ]]; then
  CHAMBERS="${CHAMBERS_3}"
  Z="${Z_3}"
  TARGET_INDEX="${TARGET_INDEX_3}"
elif [[ "${NCH}" == "4" ]]; then
  CHAMBERS="${CHAMBERS_4}"
  Z="${Z_4}"
  TARGET_INDEX="${TARGET_INDEX_4}"
else
  echo "ERROR: NCH must be 2, 3, or 4"
  exit 1
fi

mkdir -p "${OUT_DIR}/plots"

# ============================================================
# (1) HV scan analysis
# ============================================================
echo "[1/2] Running HV scan analysis (${NCH} chambers) -> ${OUT_DIR}"
${PY} "${PIPELINE}" \
  --chambers ${CHAMBERS} \
  --producer "${PRODUCER}" \
  --producer-mode optionc \
  --z ${Z} \
  --anchor-index "${ANCHOR_INDEX}" \
  --target-index "${TARGET_INDEX}" \
  --tol-x "${TOL_X}" \
  --tol-y "${TOL_Y}" \
  --tol-t "${TOL_T}" \
  --hit-mode "${HIT_MODE}" \
  --require-ncluster-eq1 \
  --out-dir "${OUT_DIR}" \
  --fit

# ============================================================
# (2) CMS-style plotting
# ============================================================
echo "[2/2] Making CMS-style plot -> ${OUT_DIR}/plots"

if [[ "${USE_ROI}" == "1" ]]; then
  ${PY} "${PLOTTER}" \
    --use-roi \
    --series "ROI::${OUT_DIR}/summaries/summary_HV*.json" \
    --out "${OUT_DIR}/plots/eff_roi_cms_style.png"
else
  ${PY} "${PLOTTER}" \
    --series "Global::${OUT_DIR}/summaries/summary_HV*.json" \
    --out "${OUT_DIR}/plots/eff_global_cms_style.png"
fi

echo "Done."
echo "Outputs:"
echo "  - Summaries: ${OUT_DIR}/summaries/"
echo "  - ROOT maps: ${OUT_DIR}/maps/"
echo "  - Plots:     ${OUT_DIR}/plots/"

