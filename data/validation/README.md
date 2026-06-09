# data/validation

Validation logs are the main place to check current audit status before changing data.

## Contents

- `csv_column_notes_index_2026-05-06.csv` (file): Index of bottom-of-file column-note blocks appended to project CSV files. Approximate data rows: 1207.
- `hayashi_equivalent_sample_rules_2026-05-07.csv` (file): Current per-galaxy candidate rules for deriving Hayashi-equivalent modeling samples. Approximate data rows: 27.
- `jsimon_full_rebuild_report_2026-05-06.csv` (file): Report for rebuilding Simon & Geha/jsimon-based galaxy CSVs from complete local tables. Approximate data rows: 6.
- `member_count_discrepancy_log_2026-05-06.csv` (file): Machine-readable member-count discrepancy audit against Hayashi Table 1 n. Approximate data rows: 27.
- `member_count_discrepancy_log_2026-05-06.md` (file): Human-readable member-count discrepancy audit and follow-up notes.
- `member_flag_normalization_by_source_2026-05-06.csv` (file): Per-source summary of member_flag normalization results. Approximate data rows: 32.
- `member_flag_normalization_changes_2026-05-06.csv` (file): Row-level or source-level changes made during member_flag normalization. Approximate data rows: 1746.
- `member_flag_normalization_summary_2026-05-06.csv` (file): Summary report for the member_flag normalization pass. Approximate data rows: 27.
- `negative_member_count_missing_member_rows_2026-05-06.csv` (file): Rows suspected missing during negative member-count audit; currently expected to be empty/diagnostic after rebuild. Approximate data rows: 0.
- `negative_member_count_omission_audit_2026-05-06.csv` (file): Machine-readable audit of galaxies where current member count was below Hayashi n. Approximate data rows: 1.
- `negative_member_count_omission_audit_2026-05-06.md` (file): Human-readable omission audit, including Reticulum II member-definition finding.
- `random_source_audit_2026-05-04.csv` (file): Earlier random source audit comparing selected CSV rows with source data. Approximate data rows: 10.
- `random_source_audit_2026-05-07.csv` (file): Latest fixed-seed random source audit; 12/12 sampled rows passed source-value checks. Approximate data rows: 12.
- `random_source_audit_2026-05-07.md` (file): Human-readable version of the latest random source audit.
- `step2_ultrawork_audit_2026-05-07.md` (file): Human-readable Step 2 ultrawork data-side audit and sample-rule summary.

## Notes For Agents

- Prefer the latest dated audit when conflicts appear.
- `hayashi_equivalent_sample_rules_2026-05-07.csv` is the compact current status table for sample-selection rules.
