# data/raw

Raw materials should be treated as source evidence. Prefer adding derived files elsewhere rather than editing raw snapshots.

## Contents

- `arxiv` (directory): Downloaded arXiv source archives and selected extracted source directories.
- `jsimon` (directory): Local Simon & Geha / Josh Simon machine-readable stellar kinematic tables.
- `vizier` (directory): Local VizieR ASU TSV snapshots and source mapping manifests.
- `ads_2007ApJ670313_DATA.html` (file): ADS data page snapshot for Simon & Geha 2007.
- `iop_10.1086_521816.html` (file): IOP article/data page snapshot for Simon & Geha 2007.
- `SIMapj07b.pdf` (file): Supplementary/source PDF collected for Simon & Geha 2007 data tracing.
- `SIMapj07b_extracted.txt` (file): Text extracted from SIMapj07b.pdf for search and audit.
- `stacks_670_313.html` (file): IOP/ADS stack or article snapshot for ApJ 670, 313.

## Notes For Agents

- Treat files here as immutable evidence and record derived interpretations in `processed/`, `galaxies/`, or `validation/`.
