# Trial Membership Classification Audit (2026-06-14)

## Scope

This audit implements a source-aware trial classification for the 14 galaxies
outside the frozen Hayashi-equivalent backup set. It does not overwrite the
raw source snapshots. Semantic repairs are applied only to the project galaxy
CSVs, with every changed row recorded in
`source_member_semantics_repair_2026-06-14.csv`.

Work continued after the initial pause. The two original pause conditions are
now carried as explicit issues rather than silently forcing a numerical
answer. The machine-readable issue list is
`trial_classification_issue_log_2026-06-14.csv`.

Scripts:

- `scripts/repair_source_membership_semantics.py`
- `scripts/rebuild_hercules_full_catalog.py`
- `scripts/build_trial_membership_catalogs.py`

Derived outputs:

- `data/processed/membership/*_trial_classified_observations.csv`
- `data/processed/membership/*_trial_unique_stars.csv`
- `data/processed/membership/trial_membership_summary_2026-06-14.csv`

## Critical Semantic Repairs

1. Grus 1: four source rows labelled `CM` were changed from
   `member_flag=0` to `member_flag=2`. Chiti et al. (2022) define `CM` as a
   candidate member with compatible velocity but no metallicity measurement.
2. Segue 2: four Kirby et al. (2013) rows labelled `B` were changed from
   `member_flag=2` to `member_flag=1`. The VizieR ReadMe defines `B` as a
   confirmed horizontal-branch member.
3. Bootes I: `Boo1_9`, `Boo1_35`, `Boo1_106`, and `Boo1_113` were changed
   from `member_flag=1` to `member_flag=2`. They have `Member=1` in the
   published VizieR catalogue but are absent from both the final `M` and
   `VCNM` lists in Jenkins et al. (2021) Table 2. This makes the final member
   count 69, matching the paper.
4. Tucana 2: all 60 Table 6 observation rows were changed from
   `member_flag=2` to `member_flag=1`. Chiti et al. (2023) explicitly describe
   the appendix table as a compilation of velocity measurements of Tucana II
   members. The rows represent 19 unique member stars.
5. Hercules: the previous 21-row Kirby metallicity-overlap convenience table
   was replaced by all 86 rows of `jsimon_data/Herc_feh.dat`: 30 source hard
   members and 56 nonmembers. Kirby metallicities are preserved for the 21
   overlap stars; the remaining metallicities retain explicit Simon & Geha
   upstream provenance. The row-by-row reconstruction is recorded in
   `hercules_full_source_rebuild_2026-06-14.csv`.

## Trial Rules

- Probability-only source tables use the strict rule `P_mem > 0.95`.
  Lower-probability rows are labelled `not_selected_probability`, not forced
  to hard nonmembers.
- Hard source classifications map to confirmed, nonmember, or unclassified.
- Membership and velocity usability are independent.
- Source-identified binaries and RR Lyrae stars remain members but their
  velocities are excluded from the baseline dynamical sample.
- Same-system multi-epoch velocities are combined by inverse-variance
  weighting only if a constant-velocity chi-square test gives `p >= 0.01`.
- Measurements from multiple catalogues or multiple instruments are not
  averaged until their velocity zero points are resolved.
- The current identity key uses coordinates to seven decimal degrees. This is
  adequate for exact repeated coordinates but is not the final cross-source
  matcher.

## Trial Results

| Galaxy | Unique stars | Confirmed | Confirmed velocity exclusions | Baseline dynamical |
|---|---:|---:|---:|---:|
| Antlia 2 | 726 | 283 | 7 | 276 |
| Bootes I | 118 | 69 | 5 | 64 |
| Crater 2 | 207 | 141 | 3 | 138 |
| Draco 2 | 51 | 14 | 0 | 14 |
| Grus 1 | 69 | 8 | 2 | 6 |
| Grus 2 | 247 | 21 | 2 | 19 |
| Hercules | 86 | 30 | 0 | 30 |
| Leo IV | 104 | 20 | 2 | 18 |
| Leo V | 105 | 11 | 3 | 8 |
| Reticulum II | 25 | 18 | 0 | 18 |
| Segue 2 | 960 | 27 | 1 | 26 |
| Tucana 2 | 19 | 19 | 2 | 17 |
| Tucana 3 | 106 | 26 | 0 | 26 |
| Tucana 4 | 192 | 11 | 2 | 9 |

The baseline counts are provisional wherever cross-source or cross-instrument
measurements remain unresolved. The Jenkins validation cases reproduce the
paper exactly:

- Bootes I: 69 total members, 64 used for velocity dispersion.
- Leo IV: 20 total members, 18 used for velocity dispersion.
- Leo V: 11 total members, 8 used for velocity dispersion.

## Unresolved and Recorded Issues

### Tucana II: Source Text Gives 16 but the Explicit Selection Gives 17

The source paper establishes 19 confirmed members and excludes two binary
candidates, `TucII-078` and `TucII-309`. This leaves 17 nonbinary members.
The paper's explicit instrument-priority list also enumerates 17 stars, but
the next sentence states that the dynamical sample contains 16 stars.

No star was removed to force the stated count. The provisional baseline now
implements the paper's explicit instrument priority on a star-by-star basis,
giving 17 usable nonbinary members. The statement that the sample contains 16
stars remains recorded as an internal publication inconsistency.

Primary source:

- https://arxiv.org/abs/2205.01740

### Segue 2: Cross-Source Identity and Backfilled Velocity Provenance

The two source catalogues contain 39 reciprocal nearest-neighbour target
matches within 0.2 arcsec, including four pairs classified as members in both
catalogues. Exact coordinate matching misses the pairs, while globally
rounding coordinates to five decimals incorrectly merges unrelated nearby
targets.

The implemented rule uses reciprocal nearest-neighbour matching, gives the
Kirby et al. classification and velocity priority when that study contains an
effective observation, and otherwise retains the Belokurov et al. result.
This exception is required for Belokurov member `064`: its matching Kirby row
was not observed and contains a velocity backfilled from Belokurov. The
backfilled value is not counted as an independent epoch.

The complete 39-pair audit is in
`segue2_cross_source_match_audit_2026-06-14.csv`.

Primary catalogues:

- https://vizier.cds.unistra.fr/viz-bin/cat/J/MNRAS/397/1748
- https://vizier.cds.unistra.fr/viz-bin/cat/J/ApJ/770/16

### Hercules: Complete Source Table Restored, Remaining Systematics Recorded

The former project table contained only the 21 stars shared with the Kirby
metallicity catalogue. The complete Simon & Geha machine table contains 86
observed targets and provides a hard `Mem` classification: 30 members and 56
nonmembers. All 86 targets are now retained.

Simon & Geha classified membership using velocity, CMD position, Na I
equivalent width, spatial position, fitted spectral type, metallicity, and
individual spectral inspection. Their source text explicitly describes one
included Hercules member at approximately `30.7 km/s` as lying near the edge
of the velocity, Na I, and metallicity ranges. This is source row `308_734`
in the machine table. It remains a member in the baseline because the project
does not override a hard source classification, but a sensitivity run that
excludes it is required.

Further recorded limitations:

- The source table contains one representative velocity per target and does
  not provide a Hercules-specific multi-epoch binary screen.
- Source star `309_528` has an implausible upstream `[Fe/H]=+1.12`, while its
  Kirby overlap metallicity is `-2.97 +/- 0.16`. The rebuilt table preserves
  the Kirby value and records the provenance conflict.
- Hayashi Table 1 lists 18 stars for Hercules. The unified project sample is
  intentionally source-defined rather than count-matched to Hayashi, so the
  30-member table is not trimmed to reproduce 18.

Primary source:

- https://arxiv.org/abs/0706.0516

## Verification

- All three scripts pass Python byte-code compilation.
- The semantic repair script is repeatable; a second run validates repaired
  rows without changing their flags.
- Assertions passed for Jenkins et al. (2021): `69/64`, `20/18`, and `11/8`.
- Tucana 2 has 60 observation rows, 19 unique confirmed members, and two
  source-identified binary candidates.
- Grus 1 reproduces the source dynamical selection of 8 members minus 2 binary
  candidates, leaving 6.
- Grus 2 reproduces 21 members minus 2 clear binaries, leaving 19.
- Tucana 4 reproduces 11 members minus 2 clear binaries, leaving 9.
- Segue 2 has 27 unique confirmed members across both sources and 26 usable
  velocities after excluding the RR Lyrae.
- Hercules contains all 86 source rows, with exactly 30 hard members and 56
  nonmembers; its rebuild script succeeds on repeated runs.
