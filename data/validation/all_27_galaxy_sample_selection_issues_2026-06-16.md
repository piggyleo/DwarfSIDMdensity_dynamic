# All 27 Galaxy Stellar-Sample Selection Issues (2026-06-16)

This note summarizes, galaxy by galaxy, whether the current stellar-table
selection has known issues. It focuses on membership screening, duplicate or
multi-epoch handling, binary/velocity-variable handling, and source-definition
conflicts.

Primary upstream project files:

- `data/validation/hayashi_equivalent_sample_rules_2026-05-07.csv`
- `data/validation/member_count_discrepancy_log_2026-05-06.csv`
- `data/validation/trial_classification_issue_log_2026-06-14.csv`
- `data/validation/trial_membership_classification_2026-06-14.md`
- `data/processed/membership/trial_membership_summary_2026-06-14.csv`
- `PROJECT_MEMORY.md`

Terminology:

- **Member count** means the number of unique stars classified as members by
  the currently adopted source-aware rule.
- **Baseline dynamical count** means confirmed members with a selected usable
  velocity after excluding known binaries, velocity variables, RR Lyrae stars,
  or unresolved velocity-source conflicts.
- **Frozen Hayashi-equivalent set** means the 13 galaxies whose current default
  modeling sample still uses the Hayashi-matching rule `member_flag==1` with
  valid velocity and velocity errors.
- **Unified trial set** means the 14 galaxies reclassified with the newer rule:
  prefer the source paper's membership semantics and do not force the member
  count to match Hayashi Table 1.

## Issue-Type Summary

| Issue type | Galaxies |
|---|---|
| No current stellar-sample selection problem | Canes Venatici I, Canes Venatici II, Coma Berenices, Eridanus II, Horologium I, Hydra II, Leo T, Pisces II, Ursa Major II |
| No current problem, but minor implementation caution | Ursa Major I, Willman 1 |
| Uncertain rows excluded by default | Segue 1, Triangulum II |
| Probability-threshold membership rule | Antlia 2, Crater 2 |
| Source catalogue versus source-paper conflict already repaired | Bootes I |
| Single-epoch or insufficient binary screening | Draco 2, Hercules |
| Borderline or candidate-member sensitivity | Hercules, Reticulum II |
| Multi-epoch velocity-variable or binary handling | Grus 1, Grus 2, Tucana 4 |
| Cross-source identity and velocity provenance | Segue 2 |
| Source-text internal count inconsistency | Tucana 2 |
| Tidally disrupting / very small dispersion systematic | Tucana 3 |

## Per-Galaxy Notes

### 1. Antlia 2

- Current rule: use Ji et al. table4, select stars with `P_mem > 0.95`, and
  prefer the Ji consistently reprocessed velocity product over older overlapping
  Torrealba rows.
- Current result in unified trial set: 283 members, 276 baseline dynamical
  velocities after excluding 7 confirmed binary candidates.
- Problem: the `P_mem > 0.95` threshold is a project rule, not a hard member
  boundary explicitly defined by the source paper. The table stores many
  probability rows rather than a simple hard member flag.
- Required follow-up: keep the threshold configurable and run sensitivity tests
  around the adopted probability cut.

### 2. Bootes I

- Current rule: use Jenkins et al. final paper-level `M` member list rather than
  trusting the raw VizieR `Member=1` column blindly.
- Current result in unified trial set: 69 members, 64 baseline dynamical
  velocities after excluding four binary candidates and one RR Lyrae star.
- Problem: four VizieR rows, `Boo1_9`, `Boo1_35`, `Boo1_106`, and `Boo1_113`,
  are marked `Member=1` in the catalogue but are absent from the paper's final
  `M`/`VCNM` Table 2 list. They were changed to `member_flag=2` with an audit
  trail.
- Current status: no unresolved membership-selection problem, but retain the
  catalogue-versus-paper conflict audit.

### 3. Canes Venatici I

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 214 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem. The table was rebuilt
  from the Simon & Geha / Josh Simon machine table, while `source_ref_id=5`
  preserves the Hayashi Table 1 reference chain through Kirby et al. 2013b.
- Implementation caution: this provenance split should not be mistaken for a
  membership conflict.

### 4. Canes Venatici II

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 25 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.
- Implementation caution: as for CVn I, the kinematic machine-table provenance
  is recorded separately from the Hayashi Table 1 primary reference chain.

### 5. Coma Berenices

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 59 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.
- Implementation caution: none specific to membership selection.

### 6. Crater 2

- Current rule: use Ji et al. table5 and select stars with `P_mem > 0.95`.
- Current result in unified trial set: 141 members, 138 baseline dynamical
  velocities after excluding 3 source binary candidates.
- Problem: the source table provides membership probabilities rather than a
  hard source member flag. The adopted threshold is a project baseline choice.
- Required follow-up: keep the probability threshold configurable and perform
  threshold sensitivity tests.

### 7. Draco 2

- Current rule: use Longeard et al. final hard member classification.
- Current result in unified trial set: 14 members and 14 baseline dynamical
  velocities.
- Problem: the current table provides no multi-epoch binary screen for the
  selected members. The stars are therefore single-epoch untested, not proven
  non-variable.
- Hayashi note: Hayashi Table 1 lists 9 stars, but the unified rule currently
  follows the source hard classification and does not force the count to 9.

### 8. Eridanus II

- Current rule: frozen Hayashi-equivalent default, using the source-recorded
  kinematic-analysis final selection.
- Current result: 92 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.
- Implementation caution: this source has not yet gone through the same
  detailed binary-sensitivity audit as the 14 unified trial galaxies.

### 9. Grus 1

- Current rule: use Chiti et al. confirmed `M` members; source `CM` rows are
  candidate members and are not included in the baseline.
- Current result in unified trial set: 8 members, 6 baseline dynamical
  velocities after excluding 2 binary candidates.
- Problem: the table contains repeated observations. One of the two binary
  candidate decisions is recoverable only after using external M2FS information
  not present as rows in the current galaxy CSV.
- Required follow-up: preserve the source binary label; later ingest external
  epochs if a joint multi-instrument velocity-offset model is implemented.

### 10. Grus 2

- Current rule: use Simon et al. 21 hard members, prefer IMACS table3 for the
  baseline velocity product, and do not average AAT and IMACS velocities.
- Current result in unified trial set: 21 members, 19 baseline dynamical
  velocities.
- Problem: two confirmed members are clear velocity variables and are excluded:
  `J220352.01-462446.5` with `p = 1.0e-8`, and `J220433.75-462639.8` with
  `p = 4.6e-5`.
- Additional sensitivity issue: three confirmed members have weaker evidence
  for velocity variability but remain above the current exclusion threshold
  `p < 0.01`: `J220420.12-462341.3` with `p = 0.020`,
  `J220437.47-461955.0` with `p = 0.030`, and
  `J220403.28-462520.0` with `p = 0.083`.
- Required follow-up: keep the 19-star baseline, and run a sensitivity sample
  excluding the weak-variable candidates.

### 11. Hercules

- Current rule: use the complete Simon & Geha machine table hard `Mem`
  classification, based on velocity, CMD position, Na I equivalent width,
  spatial position, fitted spectral type, metallicity, and visual spectral
  inspection.
- Current result in unified trial set: 30 members and 30 baseline dynamical
  velocities from 86 total source targets.
- Problem 1: no Hercules-specific multi-epoch binary screen is available. These
  are single-epoch untested velocities.
- Problem 2: source member `308_734`, with velocity near `30.7 km/s`, is
  explicitly described by the source paper as near the edge of the velocity,
  Na I, and metallicity membership ranges. It remains in the baseline because
  the project does not override source hard classifications.
- Problem 3: star `309_528` has an implausible upstream Simon & Geha
  metallicity `[Fe/H]=+1.12`, while the Kirby overlap value is
  `-2.97 +/- 0.16`. The project table preserves the Kirby value and records
  the conflict.
- Required follow-up: run sensitivity tests excluding `308_734`; do not claim
  the 30 velocities are binary-clean.

### 12. Horologium I

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 5 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.

### 13. Hydra II

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 13 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.

### 14. Leo IV

- Current rule: use Jenkins et al. final subjective member list, not a
  probability threshold chosen to match Hayashi.
- Current result in unified trial set: 20 members, 18 baseline dynamical
  velocities after excluding one binary candidate and one RR Lyrae star.
- Problem: no unresolved member-screening problem. The important detail is that
  the paper's final hard classification is used instead of the earlier
  Hayashi-matching probability threshold.

### 15. Leo V

- Current rule: use Jenkins et al. final subjective member list.
- Current result in unified trial set: 11 members, 8 baseline dynamical
  velocities after excluding two binary candidates and one RR Lyrae star.
- Problem: no unresolved member-screening problem. As with Leo IV, the current
  rule uses the source paper's final member classification instead of forcing
  the Hayashi count.

### 16. Leo T

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 19 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.
- Implementation caution: as with the other Simon & Geha / Kirby-chain tables,
  provenance and primary Hayashi reference chain are both recorded and should
  not be confused.

### 17. Pisces II

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 7 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.

### 18. Reticulum II

- Current rule: use Koposov et al. confirmed members only.
- Current result in unified trial set: 18 members and 18 baseline dynamical
  velocities.
- Problem: the source table contains 25 candidate/source rows, and Hayashi
  Table 1 also lists 25, but the source paper identifies only 18 confirmed
  members. The remaining 7 rows are not treated as baseline members.
- Required follow-up: keep the 7 non-confirmed candidate rows outside the
  baseline, and include them only in an explicit sensitivity sample if desired.

### 19. Segue 1

- Current rule: frozen Hayashi-equivalent default, use only `member_flag==1`
  rows with valid velocity and velocity errors; exclude `member_flag==2`.
- Current result: 71 modeling members, matching Hayashi Table 1.
- Problem: 129 additional rows have `member_flag=2`; these are blank or
  non-hard source membership states and are not assumed to be members.
- Required follow-up: no immediate issue for the default Hayashi-equivalent
  sample, but do not include `member_flag=2` in the likelihood unless an
  explicit alternative membership model is defined.

### 20. Segue 2

- Current rule: merge Belokurov and Kirby source catalogues with a
  reciprocal-nearest-neighbour identity audit. Prefer Kirby classification and
  velocity where Kirby has an effective observation; otherwise retain the
  Belokurov-only confirmed member.
- Current result in unified trial set: 27 members, 26 baseline dynamical
  velocities after excluding the Kirby RR Lyrae.
- Problem 1: the two source catalogues contain 39 reciprocal cross-source target
  matches within 0.2 arcsec. Exact-coordinate matching misses real duplicates,
  while coarse coordinate rounding can incorrectly merge unrelated nearby
  targets.
- Problem 2: one Kirby row contains a velocity backfilled from Belokurov. That
  value must not be treated as an independent Kirby epoch.
- Current severity: moderate but manageable, because the project now has a
  dedicated cross-source match audit and an explicit velocity-provenance rule.
- Required follow-up: retain the full match audit and never count backfilled
  values as independent measurements.

### 21. Triangulum II

- Current rule: frozen Hayashi-equivalent default, use the 13
  `member_flag==1` rows with valid velocity and velocity errors; exclude the
  one `member_flag==2` row.
- Current result: 13 modeling members, matching Hayashi Table 1.
- Problem: the sole `member_flag=2` row, `[MIC2016] 25`, has a special source
  flag and no radial velocity. It cannot enter a dynamical likelihood.
- Current status: no unresolved issue for the default sample.

### 22. Tucana 2

- Current rule: Chiti et al. Table 6 is treated as a compilation of velocity
  measurements of 19 confirmed Tucana II members. Measurements are grouped by
  unique star, and the source paper's explicit preferred instrument per star
  is used before velocity merging.
- Current result in unified trial set: 19 members, 17 baseline dynamical
  velocities after excluding `TucII-078` and `TucII-309` as known binary
  candidates.
- Problem: the source text gives an internal count conflict. The explicit
  member list contains 19 members. Excluding the two named binary candidates
  leaves 17 nonbinary stars. The explicit per-star preferred-instrument list
  also enumerates 17 stars, but the paper states elsewhere that the dynamical
  sample contains 16 stars.
- Current project decision: do not delete an unidentified 17th star just to
  match the text. Use the auditable 17-star rule provisionally and record the
  `16/17` inconsistency as unresolved.
- Required follow-up: locate any hidden quality cut or erratum that identifies
  the missing 16th-star exclusion; otherwise treat 16 as a likely text/count
  inconsistency and report a 17-star baseline plus sensitivity note.

### 23. Tucana 3

- Current rule: use Simon et al. 26 hard members and merge same-system repeated
  velocities after velocity-constancy checks.
- Current result in unified trial set: 26 members and 26 baseline dynamical
  velocities.
- Problem: no source-identified member binary is excluded. However, Tucana III
  is tidally disrupting and has an extremely small velocity dispersion, so
  undetected binaries and tidal streaming can be important physical
  systematics even if the stellar sample is source-reproduced.
- Required follow-up: retain all 26 in the source-reproduced baseline and run
  binary/tidal sensitivity analyses where possible.

### 24. Tucana 4

- Current rule: use Simon et al. table3 hard members; table2 AAT rows without
  hard membership are retained as provenance but not used in the baseline.
- Current result in unified trial set: 11 members, 9 baseline dynamical
  velocities after excluding two clear binaries.
- Problem: a third member has a velocity-variability probability around
  `p = 0.03`. The source baseline excludes only the two clear binaries, so this
  weaker candidate remains in the current baseline.
- Required follow-up: keep the 9-star source baseline and run a sensitivity
  sample excluding the `p = 0.03` star.

### 25. Ursa Major I

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 39 modeling members, matching Hayashi Table 1.
- Problem: no current member-sample selection problem. Older velocity-gap logs
  mention one nonmember or metallicity-source row without a velocity; this does
  not affect the 39-member likelihood sample.
- Implementation caution: the table uses Simon & Geha / Josh Simon machine
  data with the Hayashi Table 1 reference chain retained.

### 26. Ursa Major II

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 20 modeling members, matching Hayashi Table 1.
- Problem: no current stellar-sample selection problem.

### 27. Willman 1

- Current rule: frozen Hayashi-equivalent default, `member_flag==1` with valid
  velocity and velocity errors.
- Current result: 40 modeling members, matching Hayashi Table 1.
- Problem: no current membership-screening problem. However, coordinate-based
  de-duplication must be careful: two stars can collide under coarse coordinate
  rounding, so source IDs should be preserved for identity.
- Implementation caution: use source star IDs rather than rounded coordinates
  alone when constructing the likelihood sample.

