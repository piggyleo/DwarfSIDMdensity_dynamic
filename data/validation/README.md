# data/validation

Validation logs are the main place to check current audit status before changing data.

## Contents

- `all_27_galaxy_sample_selection_issues_2026-06-16.md` (file): Current all-galaxy summary of stellar-sample selection rules, known issues, and required follow-up for all 27 Hayashi targets.
- `csv_column_notes_index_2026-05-06.csv` (file): Index of bottom-of-file column-note blocks appended to project CSV files. Approximate data rows: 1207.
- `frozen_13_model_member_flag_precheck_2026-06-16.md` (file): Precheck and migration audit for the 13 frozen Hayashi-equivalent galaxies, including whether scientific re-screening was needed.
- `galaxy_image_record_2026-05-15.md` (file): Per-galaxy image record with structural centers and links to locally stored source-figure screenshots.
- `hayashi_equivalent_sample_rules_2026-05-07.csv` (file): Current per-galaxy candidate rules for deriving Hayashi-equivalent modeling samples. Approximate data rows: 27.
- `hercules_full_source_rebuild_2026-06-14.csv` (file): Row-level audit of the 86-target Hercules rebuild, including previous overlap status, hard membership, kinematics, metallicity provenance, and action.
- `jsimon_full_rebuild_report_2026-05-06.csv` (file): Report for rebuilding Simon & Geha/jsimon-based galaxy CSVs from complete local tables. Approximate data rows: 6.
- `member_count_discrepancy_log_2026-05-06.csv` (file): Machine-readable member-count discrepancy audit against Hayashi Table 1 n. Approximate data rows: 27.
- `member_count_discrepancy_log_2026-05-06.md` (file): Human-readable member-count discrepancy audit and follow-up notes.
- `member_flag_normalization_by_source_2026-05-06.csv` (file): Per-source summary of member_flag normalization results. Approximate data rows: 32.
- `member_flag_normalization_changes_2026-05-06.csv` (file): Row-level or source-level changes made during member_flag normalization. Approximate data rows: 1746.
- `member_flag_normalization_summary_2026-05-06.csv` (file): Summary report for the member_flag normalization pass. Approximate data rows: 27.
- `model_member_flag_application_2026-06-16.csv` (file): Row-level audit of `model_member_flag` assignment for the 14 unified-trial galaxies.
- `model_member_flag_frozen_13_2026-06-16.csv` (file): Row-level audit of `model_member_flag` assignment for the 13 frozen Hayashi-equivalent galaxies.
- `negative_member_count_missing_member_rows_2026-05-06.csv` (file): Rows suspected missing during negative member-count audit; currently expected to be empty/diagnostic after rebuild. Approximate data rows: 0.
- `negative_member_count_omission_audit_2026-05-06.csv` (file): Machine-readable audit of galaxies where current member count was below Hayashi n. Approximate data rows: 1.
- `negative_member_count_omission_audit_2026-05-06.md` (file): Human-readable omission audit, including Reticulum II member-definition finding.
- `random_source_audit_2026-05-04.csv` (file): Earlier random source audit comparing selected CSV rows with source data. Approximate data rows: 10.
- `random_source_audit_2026-05-07.csv` (file): Latest fixed-seed random source audit; 12/12 sampled rows passed source-value checks. Approximate data rows: 12.
- `random_source_audit_2026-05-07.md` (file): Human-readable version of the latest random source audit.
- `segue2_cross_source_match_audit_2026-06-14.csv` (file): All 39 reciprocal-nearest-neighbour matches between the two Segue 2 source catalogues, including identity and velocity-source resolution.
- `source_member_semantics_repair_2026-06-14.csv` (file): Row-level audit of source membership labels repaired after checking their source-specific meanings.
- `step2_ultrawork_audit_2026-05-07.md` (file): Human-readable Step 2 ultrawork data-side audit and sample-rule summary.
- `structural_center_search_log_2026-05-14.md` (file): Per-galaxy structural-center coordinates, source references, verification status, and unresolved center issues.
- `trial_classification_issue_log_2026-06-14.csv` (file): Machine-readable issue register for the 14 source-aware trial membership classifications.
- `trial_membership_classification_2026-06-14.md` (file): Human-readable rules, results, source decisions, unresolved issues, and verification for the trial classifications.

## Notes For Agents

- Prefer the latest dated audit when conflicts appear.
- For a project-wide overview of stellar-sample screening issues, start with `all_27_galaxy_sample_selection_issues_2026-06-16.md`.
- For the current unified member/multi-epoch trial, start with `trial_membership_classification_2026-06-14.md` and its issue-log CSV.
- `hayashi_equivalent_sample_rules_2026-05-07.csv` is the compact current status table for sample-selection rules.
