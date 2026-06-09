# data/raw/vizier

Files are ASU TSV responses with metadata comments. Parse carefully; multi-table responses contain several table blocks.

## Contents

- `galaxy_local_source_map.csv` (file): Mapping from Hayashi galaxies/ref IDs to local VizieR files and missing refs. Approximate data rows: 27.
- `J_AJ_165_55.tsv` (file): VizieR snapshot for Chiti et al. 2023, Tucana II velocity/abundance tables including multi-epoch table6.
- `J_ApJ_674_L81.tsv` (file): VizieR snapshot for Kuehn et al. 2008, CVn I RR Lyrae/distance reference.
- `J_ApJ_710_1664.tsv` (file): VizieR snapshot for de Jong et al. 2010, Leo IV/Leo V photometric source.
- `J_ApJ_733_46.tsv` (file): VizieR snapshot for Simon et al. 2011, Segue 1 spectroscopy and membership table.
- `J_ApJ_770_16.tsv` (file): VizieR snapshot for Kirby et al. 2013a, Segue 2 target/member table.
- `J_ApJ_779_102.tsv` (file): VizieR snapshot for Kirby et al. 2013b metallicity compilation; Hayashi Table 1 ref5.
- `J_ApJ_838_11.tsv` (file): VizieR snapshot for Simon et al. 2017, Tucana 3 multi-epoch velocity table.
- `J_ApJ_838_83.tsv` (file): VizieR snapshot for Kirby et al. 2017, Triangulum II spectroscopy tables.
- `J_ApJ_892_137.tsv` (file): VizieR snapshot for Simon et al. 2020, Grus 2/Tucana 4 velocity tables.
- `J_ApJ_920_92.tsv` (file): VizieR snapshot for Jenkins et al. 2021, Bootes I/Leo IV/Leo V member probabilities.
- `J_ApJ_921_32.tsv` (file): VizieR snapshot for Ji et al. 2021, Antlia 2 and Crater 2 probability/velocity tables.
- `J_ApJS_191_352.tsv` (file): VizieR snapshot for an auxiliary ApJS source used during source discovery/traceability.
- `J_MNRAS_397_1748.tsv` (file): VizieR snapshot for Belokurov et al. 2009, Segue 2 supplementary spectroscopy source.
- `J_MNRAS_488_2743.tsv` (file): VizieR snapshot for Torrealba et al. 2019, Antlia 2 spectroscopic/modeling source.
- `manifest.csv` (file): Manifest of downloaded VizieR sources, URLs, local file paths, and table counts. Approximate data rows: 13.

## Notes For Agents

- These files support provenance checks; changes here should be rare and source-preserving.
