#!/usr/bin/env bash
set -euo pipefail

python3 src/rpc_cms_style_eff_plot.py \
  --series 'Beam ROI::out_hv/summaries/summary_HV*_beam.json' \
  --use-roi \
  --out out_hv/plots/cms_style_beam_roi.png
