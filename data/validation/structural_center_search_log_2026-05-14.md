# Structural Center Search Log (2026-05-14)

Scope: start collecting galaxy structural centers needed to replace the current loader's temporary member-star-median origin. The target list is the 27 galaxies in Hayashi et al. 2023 Table 1.

## Outputs

- Created `data/processed/galaxy_structural_centers.csv`.
- Each row records center coordinates, center source reference id, source citation/table, source URL or local path, status, and notes.
- Column notes are appended at the bottom of the CSV; read with `pd.read_csv(path, comment="#")`.

## Extraction Policy

- Prefer a Hayashi Table 1 reference for each galaxy.
- For galaxies with Hayashi ref2, use Muñoz et al. 2018 because Hayashi Table 1 cites it as the compiled structural parameter source.
- Use Muñoz Table 1 for centers and Muñoz Table 3 Plummer rows for position angle and ellipticity when available.
- For objects without ref2, use the relevant discovery/structure reference listed by Hayashi Table 1 and mark rows `needs_primary_table_audit` until the original table has been directly extracted.
- Keep Hayashi `qprime` in a separate column. Do not overwrite it with source-derived `q=1-e` until conflicts are resolved.

## Current Status Counts

After the continuation pass on the initially flagged objects:

- `primary_source_extracted`: 21 galaxies.
- `center_extracted_q_mismatch`: 1 galaxy, Canes Venatici I.
- `center_extracted_minor_q_rounding`: 1 galaxy, Coma Berenices.
- `center_primary_shape_hayashi_ref9_extracted`: 1 galaxy, Crater 2.
- `primary_source_extracted_upper_limit`: 1 galaxy, Grus 2.
- `center_primary_shape_source_unresolved`: 1 galaxy, Tucana 3.
- `primary_source_extracted_q_conflict_with_hayashi`: 1 galaxy, Tucana 4.

Shape-status counts mirror the same risk classes:

- `primary_source_extracted`: 21 galaxies.
- `primary_source_extracted_q_mismatch`: Canes Venatici I.
- `primary_source_extracted_minor_q_rounding`: Coma Berenices.
- `hayashi_ref9_shape_extracted`: Crater 2.
- `ellipticity_upper_limit_no_pa`: Grus 2.
- `hayashi_qprime_consistent_but_source_unresolved`: Tucana 3.
- `primary_source_extracted_q_conflict_with_hayashi`: Tucana 4.

## Source Notes

- Muñoz et al. 2018 PDF downloaded to `data/raw/Munoz_2018_ApJ_860_66.pdf`; extracted text is `data/raw/Munoz_2018_ApJ_860_66_pypdf.txt`.
- Antlia 2 values came from local VizieR table `data/raw/vizier/J_MNRAS_488_2743.tsv`, Hayashi ref1.
- Reticulum II and Horologium I also have values in Koposov et al. 2015b, but Hayashi includes ref2 and the Hayashi `qprime` values match Muñoz Plummer ellipticities, so the current adopted center source is ref2.
- Crater 2 now separates center and shape provenance. The center is from Torrealba et al. 2016, while the Hayashi-equivalent ellipticity/PA is from Ji et al. 2021 because Torrealba reports only a low-ellipticity upper limit and Ji's `e=0.12±0.02`, `PA=135±4 deg` matches Hayashi `qprime=0.88`.
- Draco 2 is now directly tied to Laevens et al. 2015a Table 1: center `15:52:47.6`, `+64:33:55`, ellipticity `0.24`, PA `70 deg`, implying `q=0.76`, consistent with Hayashi.
- Grus 2 is retained as Hayashi-equivalent round in the numeric fields. The source morphology gives an ellipticity upper limit (`e<0.2`) and no PA, so `ellipticity=0.0` should be interpreted as a modeling convention matching Hayashi `qprime=1.00`, not as a measured nonzero shape.
- Tucana 2 is now tied to Koposov et al. 2015a: center `RA=342.9796 deg`, `Dec=-58.5689 deg`, ellipticity `0.39`, PA `107 deg`, consistent with Hayashi `qprime=0.61`.
- Tucana 3 remains unresolved for shape provenance. The center follows Drlica-Wagner et al. 2015, and the numeric ellipticity is set to `e=0.20` only because Hayashi `qprime=0.80`. Drlica-Wagner et al. 2015 states that the Tucana III main body was required to be azimuthally symmetric in the MCMC fit and the linear tidal feature was analyzed separately, so ref16 does not provide a normal main-body PA/ellipticity row. The exact Hayashi-reference source for PA and ellipticity still needs to be found.
- Tucana 4 currently has an explicit conflict: Simon et al. 2020 structural values give `e=0.39` and therefore `q=0.61`, while Hayashi Table 1 lists `qprime=0.40`. This must be resolved before modeling Tucana 4.
- Canes Venatici I has a real source/Hayashi mismatch: Muñoz et al. 2018 Plummer `e=0.44` gives `q=0.56`, while Hayashi Table 1 lists `qprime=0.61` even though `b_star=452 pc` matches the Muñoz Plummer radius. Keep both values separate.
- Coma Berenices has only a minor mismatch: Muñoz Plummer `e=0.37` gives `q=0.63`, while Hayashi lists `qprime=0.62`; this is likely rounding or convention-level.

## Next Audit Targets

1. Resolve Tucana 3 shape provenance and PA. Current row has Hayashi-equivalent `e=0.20` but no confirmed source table for PA/ellipticity.
2. Resolve Tucana 4 `qprime` conflict before modeling. Current source-derived `q=0.61` disagrees with Hayashi `qprime=0.40`.
3. Decide how to handle Canes Venatici I before precision reproduction. Current source-derived `q=0.56` disagrees with Hayashi `qprime=0.61`.
4. Treat Grus 2 as round only if the modeling goal is Hayashi-equivalent reproduction; otherwise the source only supports an ellipticity upper limit and no PA.
5. For strict code ingestion, allow separate center and shape provenance fields: `center_source_ref_id` and `shape_source_ref_id`.
