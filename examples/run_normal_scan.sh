#!/usr/bin/env bash
set -euo pipefail

python3 src/rpc_hv_pipeline.py \
  --chambers /path/to/Chamber201 /path/to/Chamber202 \
  --z 0 10 \
  --regime normal \
  --fit
