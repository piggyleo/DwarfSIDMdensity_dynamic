# Step 2 Ultrawork Data-Side Audit (2026-05-07)

Scope: continue project Step 2 under Hayashi et al. 2023 Table 1 reference policy; preserve one CSV per galaxy, all source-listed stars, asymmetric error columns, and member_flag semantics 1/0/2.

## Artifacts

- `data/validation/random_source_audit_2026-05-07.csv`
- `data/validation/random_source_audit_2026-05-07.md`
- `data/validation/hayashi_equivalent_sample_rules_2026-05-07.csv`

## Random Source Audit

- Fixed seed: `20260507`.
- Sample size: 12 randomly selected source rows from locally traceable VizieR/jsimon sources.
- Result: 12/12 PASS after source-row matching by ID plus coordinate/velocity fallback.
- Fields checked: RA, Dec, LOS velocity, LOS velocity error down/up, and member flag mapping.

Sampled galaxies: Antlia 2, Ursa Major I, Leo V, Tucana 3, Ursa Major II, Bootes I, Triangulum II, Tucana 2, Tucana 4, Crater 2, Canes Venatici I, Canes Venatici II.

## Candidate Hayashi-Equivalent Sample Rules

### candidate_unconfirmed

- `Antlia 2`: source_ref_id=9 (Ji table4); memprob>=0.5 and Bin/binary_flag == N; exclude Torrealba rows (probability_threshold; Hayashi n=283). Gives 283 in source-side audit; current CSV lacks explicit binary flag column, so implementation may need source Bin column or derived column.
- `Bootes I`: member_flag==1 and memprob>=0.95 (probability_threshold; Hayashi n=66). Jenkins subjective Member=1 has 73; threshold gives 66.
- `Crater 2`: memprob>0 (equivalently memprob>=0.1 in current table) and velocity/error non-null (probability_threshold; Hayashi n=141). Ji table has 141 nonzero/high-probability rows; member_flag remains 2 because no hard member flag.
- `Draco 2`: member_flag==1 and memprob>=0.85 and velocity/error non-null (probability_threshold; Hayashi n=9). Exact Hayashi n; stricter than source Member=Y 14 stars.
- `Leo IV`: member_flag==1 and memprob>=0.75 (probability_threshold; Hayashi n=18). Jenkins Member=1 has 20; threshold gives 18.
- `Leo V`: member_flag==1 and memprob>=0.95 (probability_threshold; Hayashi n=7). Threshold gives 7; two subjective members lack memprob and are excluded.
- `Segue 2`: member_flag==1 and velocity/error non-null; deduplicate across ref28/ref29 by coordinates (deduplicate_cross_source; Hayashi n=26). Raw flag1=27 but coordinate-dedup gives 26; decide duplicate velocity priority later.

### ready_low_risk

- `Canes Venatici I`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=214). Count matches Hayashi after full jsimon rebuild; source_ref_id now set to Hayashi Table 1 ref5, with jsimon/S&G 2007 retained as upstream machine-table trace.
- `Canes Venatici II`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=25). Count matches Hayashi after full jsimon rebuild; source_ref_id now set to Hayashi Table 1 ref5, with jsimon/S&G 2007 retained as upstream machine-table trace.
- `Coma Berenices`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=59). Count matches Hayashi after full jsimon rebuild; source_ref_id now set to Hayashi Table 1 ref5, with jsimon/S&G 2007 retained as upstream machine-table trace.
- `Eridanus II`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=92). Source table recorded as kinematic-analysis final selection.
- `Horologium I`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=5). Count matches Hayashi.
- `Hydra II`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=13). Count matches Hayashi.
- `Leo T`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=19). Count matches Hayashi after full jsimon rebuild; source_ref_id now set to Hayashi Table 1 ref5, with jsimon/S&G 2007 retained as upstream machine-table trace.
- `Pisces II`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=7). Count matches Hayashi.
- `Ursa Major I`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=39). Count matches Hayashi after full jsimon rebuild; note one nonmember velocity gap in older summary should not affect member likelihood.
- `Ursa Major II`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=20). Count matches Hayashi after full jsimon rebuild.
- `Willman 1`: member_flag==1 and velocity/error non-null (clean_default; Hayashi n=40). Count matches Hayashi.

### candidate_needs_velocity_merge

- `Grus 1`: member_flag==1; group by unique star_id/coordinates before likelihood (deduplicate_multiepoch; Hayashi n=8). 16 member rows are 8 unique stars with repeated observations.
- `Tucana 2`: source_catalog==J/AJ/165/55/table6; group all 60 velocity rows by unique star_id/coordinates (deduplicate_multiepoch_all_rows; Hayashi n=19). No hard member column; table6 has 19 unique stars.
- `Tucana 3`: member_flag==1; group by unique star_id/coordinates (deduplicate_multiepoch; Hayashi n=26). 44 member rows are 26 unique stars.
- `Tucana 4`: source_catalog==J/ApJ/892/137/table3 and member_flag==1; group by unique star_id/coordinates; exclude table2 flag2 (deduplicate_multiepoch_exclude_table2; Hayashi n=11). Table3 unique members=11; table2 AAT rows have no hard member flag.

### unresolved

- `Grus 2`: source_catalog==J/ApJ/892/137/table3 and member_flag==1; group by unique star_id/coordinates; then identify 2 exclusions (partial_deduplicate_unresolved; Hayashi n=19). Dedup gives 21, still 2 above Hayashi; source field insufficient to safely remove 2.
- `Hercules`: No safe rule yet; current member_flag==1 gives 21 (unresolved_source_mismatch; Hayashi n=18). Current rows use jsimon ref33 but Hayashi refs are 2;5;19; need locate strict 18-star sample.

### unresolved_user_decision

- `Reticulum II`: Two explicit branches: confirmed_members=member_flag==1 gives 18; candidate_rows=all 25 source rows with velocity gives Hayashi n (member_definition_mismatch; Hayashi n=25). No row omission; Hayashi n=25 conflicts with source-confirmed 18 members.

### ready_with_flag2_exclusion

- `Segue 1`: member_flag==1 and velocity/error non-null; exclude member_flag==2 (clean_default_exclude_flag2; Hayashi n=71). flag1 matches Hayashi; flag2 are blank/non-hard members.
- `Triangulum II`: member_flag==1 and velocity/error non-null; exclude member_flag==2 (clean_default_exclude_flag2; Hayashi n=13). flag2 row has no HRV and special source flag.

## Issues Found By Global Consistency Audit

- Six jsimon-rebuilt CSVs store `member_flag` as `0.0/1.0` strings in raw CSV text rather than integer-style `0/1`; semantics are correct but formatting should be normalized later.
- Resolved after audit: off-Table-1 jsimon-based rows for CVn I/II, Coma Berenices, Hercules, and Leo T now use Hayashi Table 1 `source_ref_id=5` (Kirby et al. 2013b); `source_catalog/notes` retain Simon & Geha (2007) jsimon machine tables as the upstream kinematic trace. Ursa Major I/II keep `source_ref_id=33` because Hayashi Table 1 explicitly lists ref 33 for them.
- `source_ref_id` is often serialized as float-like values such as `33.0`; downstream loaders should normalize to integer strings defensively.

## Do Not Modify Yet

- Do not overwrite original star CSVs with candidate sample selections. Preserve all source-listed rows.
- Do not force Reticulum II to 25 members until the user chooses confirmed-member vs candidate-row modeling.
- Do not trim Hercules or Grus 2 automatically; both require further source-level evidence.