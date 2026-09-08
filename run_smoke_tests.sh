#!/usr/bin/env bash
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

Rscript "$root/simulation/main/smoke_test.R"
Rscript "$root/simulation/extensions/smoke_test.R"
Rscript "$root/simulation/extensions/smoke_test_E2pp.R"
Rscript "$root/application/adhd/smoke_test_stein.R"

cd "$root/simulation/main/tests"
for f in reference_*.py; do
  python3 "$f"
done
