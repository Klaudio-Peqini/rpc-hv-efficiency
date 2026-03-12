#!/usr/bin/env bash
set -euo pipefail

python3 src/rpc_hv_pipeline.py \
  --chambers /path/to/Chamber201 /path/to/Chamber202 \
  --z 0 10 \
  --regime beam \
  --roi-name scint_overlap \
  --roi-rect 0 22 0 130 \
  --auto-cuts \
  --fit
