# Frozen 13 Galaxies: Precheck and `model_member_flag` Migration (2026-06-16)

This note records whether the 13 frozen Hayashi-equivalent galaxies need
scientific re-screening before they are migrated to the new `model_member_flag`
schema.

Migration status: completed on 2026-06-16. The scientific member samples were
not changed.

Important distinction:

- **Schema migration needed** means adding `model_member_flag`,
  `model_selection_reason`, and `model_velocity_rule` so the modeling pipeline
  can read the same column for all galaxies.
- **Scientific re-screening needed** means changing the actual selected member
  sample under the current standard.

Current loader behavior after the 2026-06-16 update:

- Default modeling uses `model_member_flag in [1, 2]`.
- If `model_member_flag` does not exist, the loader fails fast. This is
  intentional: every modeling-ready galaxy must be migrated before sampling.

## Summary

| Galaxy | Current selected members | Scientific re-screening needed? | Reason |
|---|---:|---|---|
| Canes Venatici I | 214 | No | `member_flag==1` matches Hayashi and no uncertain member rows are present. |
| Canes Venatici II | 25 | No | `member_flag==1` matches Hayashi and no uncertain member rows are present. |
| Coma Berenices | 59 | No | `member_flag==1` matches Hayashi and no uncertain member rows are present. |
| Eridanus II | 92 | No | Source table is recorded as the kinematic-analysis final selection. |
| Horologium I | 5 | No | Hard source members match Hayashi. |
| Hydra II | 13 | No | Hard source members match Hayashi. |
| Leo T | 19 | No | `member_flag==1` matches Hayashi and no uncertain member rows are present. |
| Pisces II | 7 | No | Hard source members match Hayashi. |
| Segue 1 | 71 | Minor confirmation only | `member_flag==1` matches Hayashi; 129 `member_flag==2` rows are blank/non-hard source membership states and should become `model_member_flag=3`. |
| Triangulum II | 13 | Minor confirmation only | `member_flag==1` matches Hayashi; the single `member_flag==2` row has no radial velocity and should become `model_member_flag=3`. |
| Ursa Major I | 39 | No | Member sample matches Hayashi; an older velocity-gap note concerns a nonmember/metallicity row and does not affect the 39-member likelihood sample. |
| Ursa Major II | 20 | No | `member_flag==1` matches Hayashi and no uncertain member rows are present. |
| Willman 1 | 40 | No, but preserve source IDs | Member sample matches Hayashi; avoid coarse coordinate-only de-duplication because two rows can collide under rounded coordinates. |

## Recommended Migration Actions

1. Completed: `model_member_flag=1` was added to the current
   `member_flag==1` rows with valid velocities.
2. Completed: `model_member_flag=3` was added to all nonmembers and uncertain
   rows.
3. Completed: Segue 1 and Triangulum II `member_flag==2` rows are excluded
   from the default model sample with `model_member_flag=3`.
4. Completed: Willman 1 rows were not coordinate-deduplicated; source
   `star_id` is preserved.
5. No frozen galaxy required changing its default member count under the
   present standard.

Row-level migration audit:

- `data/validation/model_member_flag_frozen_13_2026-06-16.csv`
