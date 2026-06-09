# data/galaxies

Use `pd.read_csv(path, comment="#")` because CSVs may include bottom column-note comments.

## Contents

- `01_Antlia_2.csv` (file): Standardized per-galaxy table for Antlia 2; includes one global row plus preserved star rows. Approximate data rows: 730.
- `02_Bootes_I.csv` (file): Standardized per-galaxy table for Bootes I; includes one global row plus preserved star rows. Approximate data rows: 119.
- `03_Canes_Venatici_I.csv` (file): Standardized per-galaxy table for Canes Venatici I; includes one global row plus preserved star rows. Approximate data rows: 270.
- `04_Canes_Venatici_II.csv` (file): Standardized per-galaxy table for Canes Venatici II; includes one global row plus preserved star rows. Approximate data rows: 59.
- `05_Coma_Berenices.csv` (file): Standardized per-galaxy table for Coma Berenices; includes one global row plus preserved star rows. Approximate data rows: 103.
- `06_Crater_2.csv` (file): Standardized per-galaxy table for Crater 2; includes one global row plus preserved star rows. Approximate data rows: 208.
- `07_Draco_2.csv` (file): Standardized per-galaxy table for Draco 2; includes one global row plus preserved star rows. Approximate data rows: 52.
- `08_Eridanus_II.csv` (file): Standardized per-galaxy table for Eridanus II; includes one global row plus preserved star rows. Approximate data rows: 93.
- `09_Grus_1.csv` (file): Standardized per-galaxy table for Grus 1; includes one global row plus preserved star rows. Approximate data rows: 80.
- `10_Grus_2.csv` (file): Standardized per-galaxy table for Grus 2; includes one global row plus preserved star rows. Approximate data rows: 351.
- `11_Hercules.csv` (file): Standardized per-galaxy table for Hercules; includes one global row plus preserved star rows. Approximate data rows: 22.
- `12_Horologium_I.csv` (file): Standardized per-galaxy table for Horologium I; includes one global row plus preserved star rows. Approximate data rows: 18.
- `13_Hydra_II.csv` (file): Standardized per-galaxy table for Hydra II; includes one global row plus preserved star rows. Approximate data rows: 32.
- `14_Leo_IV.csv` (file): Standardized per-galaxy table for Leo IV; includes one global row plus preserved star rows. Approximate data rows: 105.
- `15_Leo_V.csv` (file): Standardized per-galaxy table for Leo V; includes one global row plus preserved star rows. Approximate data rows: 106.
- `16_Leo_T.csv` (file): Standardized per-galaxy table for Leo T; includes one global row plus preserved star rows. Approximate data rows: 62.
- `17_Pisces_II.csv` (file): Standardized per-galaxy table for Pisces II; includes one global row plus preserved star rows. Approximate data rows: 14.
- `18_Reticulum_II.csv` (file): Standardized per-galaxy table for Reticulum II; includes one global row plus preserved star rows. Approximate data rows: 26.
- `19_Segue_1.csv` (file): Standardized per-galaxy table for Segue 1; includes one global row plus preserved star rows. Approximate data rows: 523.
- `20_Segue_2.csv` (file): Standardized per-galaxy table for Segue 2; includes one global row plus preserved star rows. Approximate data rows: 1070.
- `21_Triangulum_II.csv` (file): Standardized per-galaxy table for Triangulum II; includes one global row plus preserved star rows. Approximate data rows: 35.
- `22_Tucana_2.csv` (file): Standardized per-galaxy table for Tucana 2; includes one global row plus preserved star rows. Approximate data rows: 61.
- `23_Tucana_3.csv` (file): Standardized per-galaxy table for Tucana 3; includes one global row plus preserved star rows. Approximate data rows: 152.
- `24_Tucana_4.csv` (file): Standardized per-galaxy table for Tucana 4; includes one global row plus preserved star rows. Approximate data rows: 287.
- `25_Ursa_Major_I.csv` (file): Standardized per-galaxy table for Ursa Major I; includes one global row plus preserved star rows. Approximate data rows: 107.
- `26_Ursa_Major_II.csv` (file): Standardized per-galaxy table for Ursa Major II; includes one global row plus preserved star rows. Approximate data rows: 109.
- `27_Willman_1.csv` (file): Standardized per-galaxy table for Willman 1; includes one global row plus preserved star rows. Approximate data rows: 46.
- `_coverage_summary.csv` (file): Per-galaxy row/member/source coverage summary. Approximate data rows: 27.
- `_field_completeness_summary.csv` (file): Per-galaxy completeness summary for required fields, member flags, and source refs. Approximate data rows: 27.
- `_jsimon_velocity_fill_report.csv` (file): Report from filling velocity fields using Simon & Geha/jsimon local tables. Approximate data rows: 7.
- `_kinematic_velocity_gaps.csv` (file): Current list of per-galaxy rows still missing LOS velocity values. Approximate data rows: 3.
- `_velocity_missing_root_cause.csv` (file): Root-cause notes for remaining velocity gaps. Approximate data rows: 302.

## Notes For Agents

- Preserve one CSV per Hayashi target and keep all source-listed stars unless the user explicitly approves a derived table.
- `member_flag` semantics: `1` member, `0` non-member, `2` uncertain/no hard member flag.
- Use validation logs before deciding which rows form a Hayashi-equivalent modeling sample.
