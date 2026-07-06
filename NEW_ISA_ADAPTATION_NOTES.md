# Ep_ISA_NEW adaptation notes

This folder is an isolated adaptation of `Ep_ISA` to the updated ISA logic in:

`F:\phd\Drophila\3Model_motif_discovering\ISA\ISA_NEW\deepISA`

It keeps the Ep_ISA Fi-NeMo input wrapper and uses a separate Python package name, `Ep_ISA_NEW`, so it does not collide with the existing `Ep_ISA` package or old result folders.

## First stage: preflight motif-pair audit

`EpQuickStart.run_isa()` now starts from `preflight_audit` by default.

The audit reads `Data/motif_locs.csv` before any ISA model scoring and writes:

- `Data/preflight_motif_pair_audit.csv`
- `Data/preflight_motif_pair_audit_by_region.csv`
- `Data/preflight_overlap_or_abutting_pairs.csv`

The pair-count rule is matched to `make_pairs_for_region()`:

`same region, sorted by start_rel, distance = start2_rel - end1_rel, 1 <= distance <= receptive_field`

This answers whether motif-pair counts are already lower than single-motif counts at the motif-location input level, before single ISA filtering, null generation, or interaction scoring.

`same-region pairs` means all theoretical interval combinations inside each region: `n * (n - 1) / 2`.
Only `receptive_field_pairs` are kept by deepISA pair scoring.
The remaining skipped same-region pairs are split into:

- `overlapping_or_adjacent_pairs`: skipped because `distance < 1`;
- `too_far_pairs`: skipped because `distance > receptive_field`.

The overlap detail table records exact TF labels, coarse family labels derived from the prefix before `/`, interval coordinates, distance, overlap size, and whether the skipped pair is exact same TF, same coarse family, or different family.

## Updated ISA behavior

Compared with the previous Ep_ISA implementation:

- single null default is `8192`, matching the updated deepISA default;
- motif locations are deduplicated by `chrom,start,end,region`, keeping the highest score;
- single ISA filtering defaults to `positive_all_tracks`, retaining motifs only when the requested target track passes the positive single-null threshold;
- pair null default is `8192`;
- pair-null pseudo-site length is the median motif length from observed motif pairs, matching updated deepISA;
- raw pair scoring no longer writes final `interaction_t*` directly during `run_combi_isa()`;
- `aggregate_isa` first calls `add_interaction()`, which writes normalized `interaction_t*`;
- normalized interaction uses `tau` and sets unqualified pair rows to `NaN`;
- `calc_coop_score()` now consumes the normalized interaction table and pair null table.

Ep_ISA_NEW keeps the Ep_ISA-specific fix that `pred_orig.csv` is computed for both motif and non-motif regions, so null ISA rows from non-motif regions can be scored safely.

The positive single-filter default is important for the updated `add_interaction()` stage. The updated normalized interaction gate requires both motif singles to pass the target-track positive single-null threshold. If single ISA is filtered with the older `any_tails` rule, many negative or wrong-track motifs can enter `motif_combi_isa.csv`, and `interaction_t*` will become mostly `NaN`.

## Existing-data preflight audit

A dry audit of the current completed result folders was written to:

`F:\phd\Drophila\3Model_motif_discovering\ISA\result\new_ep_isa_preflight_audit`

Summary with `receptive_field=255`:

| Task | single motif loci after new dedup | all same-region pairs | receptive-field pairs | RF pair / single |
|---|---:|---:|---:|---:|
| CAGE | 29,459 | 26,082 | 21,523 | 0.731 |
| DEV | 34,807 | 33,517 | 28,919 | 0.831 |
| HK | 33,385 | 32,492 | 26,482 | 0.793 |

This shows that pair candidates are already fewer than single motif loci at the input motif-location level. The lower pair count is therefore not primarily caused by downstream ISA output handling.

Skipped overlap/abutting pairs are also audited because they are the main reason `same-region pairs` differ from `receptive_field_pairs` in the current data.
