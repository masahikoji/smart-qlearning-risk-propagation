# Reproducibility upload set

These files are the confirmed final-run materials currently available for the boundary-validation part of the recursive Q-learning project.

Upload these files into the repository's `reproducibility/` directory on branch `release-v3.0.0`:

- `run_table_s7_exact_aligned.py`
- `run_table_s7_paper_aligned.sh`
- `table_s7_paper_aligned_final.csv`
- `table_s7_paper_aligned_final.csv.json`
- `run_table_s7_exact_round6.py`

The filenames retain their historical working names, but the repository-level documentation is journal-neutral.

Do not create the GitHub release yet. The following items still need to be added before release:

- recursive Gaussian-limit program and saved output;
- the raw replicate arrays for the completed 20,000-replicate smooth-boundary run (`table_s7_raw_n250.npz`, `table_s7_raw_n500.npz`, `table_s7_raw_n1000.npz`);
- the final exact-tie saved output, if available;
- any environment/session record that is promised in the manuscript/archive description.
