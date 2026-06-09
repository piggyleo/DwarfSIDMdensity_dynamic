# data/processed

These files are intermediate derivations; final per-galaxy inputs live in ../galaxies/.

## Contents

- `kinematics_ref3_jenkins2021_bootes1_leo4_leo5.csv` (file): Processed extraction from Jenkins et al. 2021 for Bootes I, Leo IV, and Leo V kinematics/probabilities. Approximate data rows: 327.
- `kinematics_ref9_ji2021_antlia2_crater2.csv` (file): Processed extraction from Ji et al. 2021 for Antlia 2 and Crater 2 kinematics/probabilities. Approximate data rows: 715.
- `star_observations_master.csv` (file): Intermediate master table of star observations used while assembling per-galaxy CSVs. Approximate data rows: 3735.

## Notes For Agents

- Keep this README synchronized if files are added, removed, or repurposed.

- `galaxy_structural_centers.csv`: Candidate structural center coordinates and PA/ellipticity provenance for the 27 Hayashi et al. 2023 galaxies. Read with `pd.read_csv(path, comment="#")`; inspect `status` and `shape_status` before modeling because some rows intentionally preserve Hayashi/source conflicts, upper limits, or unresolved PA provenance.
