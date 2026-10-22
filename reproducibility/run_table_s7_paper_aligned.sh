#!/usr/bin/env bash
set -euo pipefail

# Statistical settings match smart_run_control("paper"):
# n = 250, 500, 1000; R = 20,000; base seed = 260906;
# quadrature order = 48.  Batch and jobs affect computation only.
python3 run_table_s7_exact_aligned.py \
  --n 250 500 1000 \
  --reps 20000 \
  --batch 50 \
  --jobs 3 \
  --quadrature-order 48 \
  --base-seed 260906 \
  --out table_s7_paper_aligned.csv \
  --raw-directory table_s7_raw
