#!/usr/bin/env bash
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

printf '%s\n' '== Main recursive Q-learning smoke test =='
(
  cd "$root/simulation/main"
  Rscript smoke_test.R
)

printf '%s\n' '== Independent Python reference checks =='
for script in \
  reference_v2_1_validation.py \
  reference_v2_2_decomposition_validation.py \
  reference_v2_3_paired_decomposition.py \
  reference_v2_4_regret_validation.py \
  reference_primary_truth_quadrature.py
do
  python3 "$root/simulation/main/tests/$script"
done

if command -v Rscript >/dev/null 2>&1 && \
   Rscript -e 'quit(save="no", status=if (requireNamespace("DTRlearn2", quietly=TRUE)) 0 else 1)' >/dev/null 2>&1; then
  printf '%s\n' '== Simulated ADHD primary analysis =='
  (
    cd "$root/application/adhd"
    Rscript run_adhd_sel_ma.R main
  )
else
  printf '%s\n' 'Skipping ADHD analysis: Rscript and DTRlearn2 are required.'
fi

printf '%s\n' 'All v3.0.0 lightweight checks completed.'
