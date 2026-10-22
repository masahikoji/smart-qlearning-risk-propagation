# Risk propagation under model selection and averaging in recursive Q-learning

This repository contains reproducibility materials for research on how downstream model selection and model averaging propagate through backward Q-learning via generated pseudo-outcomes.

The repository is intentionally journal-neutral. File names and documentation refer to scientific results rather than journal-specific table or section numbers so that the archived release remains valid if the manuscript is reformatted or submitted elsewhere.

## Current reproducibility materials

### Recursive Q-learning simulations

The primary simulation code is in `simulation/main/`.

It contains the audited two- and three-stage implementations used for the current risk-propagation analyses, including:

- hard information-criterion model selection and normalized exponential model averaging;
- local-signal experiments;
- smooth treatment-boundary and exact-tie scenarios;
- recursive risk-gap calculations;
- coefficient and criterion perturbation diagnostics;
- cross-stage decomposition calculations; and
- regime-value regret calculations.

The numerical-core version is recorded separately in `simulation/main/VERSION`. The repository release version and the simulation-core version are intentionally distinct.

A lightweight validation can be run with:

```bash
Rscript simulation/main/smoke_test.R
```

The production entry points are:

```bash
Rscript simulation/main/run_two_stage.R
Rscript simulation/main/run_three_stage.R
Rscript simulation/main/run_decomposition.R
Rscript simulation/main/summarize.R
```

The first three scripts default to the production (`paper`) profile when no profile argument is supplied. See `CODE_TO_RESULTS.md` for a result-oriented map.

### Recursive Gaussian-limit validation

The `reproducibility/` directory contains the independent finite-dimensional Gaussian recursion used to evaluate the local limiting risk gaps and its saved numerical outputs.

### Smooth-boundary and exact-tie validation

Exact conditional-integration code and archived outputs for the smooth-boundary and exact-tie calculations are stored with the reproducibility materials. Full-precision summary files and archived replicate-level arrays are included in the versioned release when they are used to support reported numerical values.

### Simulated ADHD SMART analysis

The application code is in `application/adhd/`.

Primary analysis:

```bash
cd application/adhd
Rscript run_adhd_sel_ma.R main
```

Candidate-library sensitivity analysis:

```bash
Rscript run_adhd_candidate_library_sensitivity.R all
```

The analysis uses the publicly distributed simulated `adhd` data from the `DTRlearn2` R package. The reported production analysis was originally run with an earlier release of `DTRlearn2`; the exact software environment used for archived runs should be documented with the corresponding output files rather than inferred from the current CRAN version.

## Validation and reproducibility

The repository includes smoke tests, reference outputs, numerical identity checks, and archived metadata for the principal computations. Large production simulations do not need to be rerun to inspect the implementation; saved summaries and, where stated in the research report, replicate-level arrays are included in the archived release.

A recommended lightweight check is to run the repository smoke tests and compare their outputs with the committed references before creating a new release.

## Historical materials

Earlier research versions are retained for traceability but are not part of the current analysis:

- `simulation/feedback_aware/`: materials associated with the v2.0.0 feedback-aware research version;
- `simulation/extensions/`: earlier exploratory materials;
- other explicitly marked legacy application scripts.

See `LEGACY_RELEASES.md` for details. Historical release tags and archived Zenodo records should not be overwritten.

## Citation and license

Use `CITATION.cff` for software citation metadata. Each archived release should be cited using its version-specific Zenodo DOI once minted.

The code in this repository is released under the MIT License. See `LICENSE`.
