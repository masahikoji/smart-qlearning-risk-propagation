# GitHub update steps for v3.0.0

## 1. Do not replace audited current analysis code unnecessarily

Keep the existing audited files under `simulation/main/` and the primary files under `application/adhd/` unless a file-level comparison identifies a real change. The principal simulation scripts were already verified against the audited project copies.

## 2. Replace repository-level metadata

Replace at repository root:

- `README.md`
- `CITATION.cff`
- `.zenodo.json`
- `VERSION`
- `CODE_TO_MANUSCRIPT.md` (replace with `CODE_TO_RESULTS.md` and remove the old file after confirming the new map is complete)

Add:

- `LICENSE`
- `LEGACY_RELEASES.md`

## 3. Add current reproducibility materials

Add the current archived materials that are not yet in the public GitHub repository, including as applicable:

- `reproducibility/compute_recursive_limit.py`
- recursive-limit CSV/raw output and metadata
- exact smooth-boundary conditional-integration code
- smooth-boundary full-precision CSV and metadata
- exact-tie validation outputs
- theoretical-constant calculation script/output
- figure-generation script under `manuscript_assets/`

Large raw replicate arrays may be uploaded directly to Zenodo instead of GitHub if repository size is a concern, but the README and `CODE_TO_RESULTS.md` should state that they are part of the archived release.

## 4. Mark older materials as legacy instead of deleting them

Retain:

- `simulation/feedback_aware/`
- `simulation/extensions/`
- historical application scripts

Add a short legacy notice to those directories if useful. Do not rewrite old Git tags or old Zenodo records.

## 5. Lightweight validation before release

Run the repository smoke tests and reference checks. There is no need to rerun all 20,000-replicate production simulations solely to create the release if the audited code and archived results are unchanged.

For the ADHD analysis, optionally verify the current CRAN release of `DTRlearn2` against the archived results. Do not claim a new package version was used unless that verification was actually performed.

## 6. Create the release

- Commit the journal-neutral v3.0.0 metadata and added reproducibility files.
- Merge to the default branch.
- Create Git tag/release `v3.0.0`.
- Let the Zenodo GitHub integration archive the release and mint a new version-specific DOI.
- Use the new version-specific DOI in the manuscript/report and submission system.

Because `.zenodo.json` and `CITATION.cff` are both present, verify `.zenodo.json` carefully before the release; Zenodo uses it in preference to `CITATION.cff` for GitHub archiving metadata.
