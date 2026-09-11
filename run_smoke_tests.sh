#!/usr/bin/env bash
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# Current feedback-aware implementation.
fa="$root/simulation/feedback_aware"
python3 -m unittest discover -s "$fa/tests" -v
python3 "$fa/run_feedback_tuning.py" smoke --workers 1 --output "$fa/results/smoke"
python3 "$fa/review/gaussian_audit.py" --output "$fa/results/supplement_gaussian"

# Primary ADHD diagnostic when R and DTRlearn2 are available.
if command -v Rscript >/dev/null 2>&1 && \
   Rscript -e 'quit(save="no", status=if (requireNamespace("DTRlearn2", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
  (cd "$root/application/adhd" && Rscript run_adhd_sel_ma.R main)
else
  printf '%s\n' 'Skipping ADHD diagnostic: Rscript and DTRlearn2 are required.'
fi
