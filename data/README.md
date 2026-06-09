# data

Top-level data directory. Read this together with PROJECT_MEMORY.md when resuming data work.

## Contents

- `galaxies` (directory): Per-galaxy standardized CSV files and summary tables for the 27 Hayashi targets.
- `processed` (directory): Intermediate processed tables derived from raw sources before final per-galaxy CSV assembly.
- `raw` (directory): Raw downloaded/source materials: VizieR snapshots, arXiv source packages, jsimon machine tables, and article snapshots.
- `validation` (directory): Audit logs, discrepancy reports, source checks, and sample-selection rule diagnostics.

## Notes For Agents

- Start with `../PROJECT_MEMORY.md` for current project state and scientific decisions.
- Do not overwrite raw evidence files unless explicitly asked.
- Per-galaxy modeling inputs should normally come from `galaxies/`, not directly from `raw/`.
