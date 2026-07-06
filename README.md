# Ep_ISA_NEW

**In Silico Ablation for TF cooperativity, adapted to consume Fi-NeMo scan results.**

Ep_ISA_NEW takes Fi-NeMo motif hits as input and runs the full deepISA ablation pipeline (single ISA → combinatorial ISA → null model → cooperativity scoring) to quantify transcription factor synergy, redundancy, and independence.

Built on [deepISA](https://github.com/anderssonlab/deepISA) (by Xuening He).

---

## Quick Start (Google Colab)

```python
# 1. Clone and install
!git clone https://github.com/<your-org>/Ep_ISA_NEW.git
%cd Ep_ISA_NEW
!pip install -e .

# 2. Import
from Ep_ISA_NEW.quickstart import EpQuickStart

# 3. Prepare region table (maps Fi-NeMo peak_id to genomic coords)
import pandas as pd
df_regions = pd.DataFrame({
    'peak_id': [0, 1, 2, ...],
    'chrom':   ['chr2R', 'chr2R', 'chr3L', ...],
    'start':   [37800, 38100, 52000, ...],
    'end':     [38049, 38349, 52249, ...],
})

# 4. Initialize
qs = EpQuickStart(
    results_dir='./results',
    fasta_path='dm6.fa',
    df_regions=df_regions,
)

# 5. Define model (train new OR load checkpoint)
qs.define_model(model_config={
    'seq_len': 249,  # match Fi-NeMo input length
    'ks': [15, 9, 9, 9, 9],
    'cs': [64, 64, 64, 64, 64],
    'ds': [1, 2, 4, 8, 16],
    'dropout': 0.1,
})

# Option A: Train from scratch
qs.train(bw_paths=['signal.bw'])

# Option B: Load pretrained
# qs.load_checkpoint('best')

# 6. Load Fi-NeMo hits (replaces map_motifs)
qs.load_finemo(
    hits_tsv_path='finemo_scans/CAGE_NEW/hits.tsv',
    finemo_h5_path='IC_Trimmed_MetaClusters_for_Finemo.h5',
    auto_threshold_percentile=50,  # auto-compute score threshold
)

# 7. Run ISA pipeline
qs.run_isa(isa_config={
    'tracks': [0],
    'null_percentile': 80,
    'min_count': 10,
    'q_val_thresh': 0.1,
})

# 8. Generate reports
qs.report()
```

---

## Architecture

```
Ep_ISA_NEW/
├── src/Ep_ISA_NEW/
│   ├── quickstart.py        EpQuickStart orchestrator
│   ├── adapter/
│   │   └── finemo_io.py     Fi-NeMo hits.tsv → motif_locs.csv conversion
│   ├── modeling/             CNN model (from deepISA)
│   ├── scoring/              ISA core algorithm (from deepISA, minus mapper.py)
│   ├── plotting/             Visualization (from deepISA)
│   ├── exploring/            Downstream analysis (from deepISA)
│   └── utils.py              Utilities (from deepISA)
├── data/                     Reference data (JASPAR, TF families, PPI, etc.)
└── notebooks/
    └── Ep_ISA_NEW_tutorial.ipynb
```

### Key difference from deepISA

| Aspect | deepISA | Ep_ISA_NEW |
|--------|---------|--------|
| Motif discovery | `map_motifs` (JASPAR bigBED scan) | `load_finemo` (Fi-NeMo hits.tsv) |
| Input coordinate system | Genomic absolute | Fi-NeMo relative → auto-converts to absolute via `df_regions` |
| Score threshold | Default 500 (JASPAR integer) | Auto-suggested percentile of `hit_coefficient` (float) |
| ISA stages | 5 (map → single → combi → null → aggregate) | 4 (skip map_motifs) |

---

## API Reference

### `EpQuickStart`

#### `__init__(results_dir, fasta_path, df_regions, device=None)`

- `df_regions`: DataFrame with columns `peak_id`, `chrom`, `start`, `end`.
  `peak_id` must match Fi-NeMo hits `peak_id` column. If absent, index is used.

#### `define_model(model_config=None, model_obj=None, mode='dual')`

- `model_config`: dict with `seq_len`, `ks`, `cs`, `ds`, `dropout`.
- `model_obj`: pass a pre-instantiated model (e.g. AlphaGenome) instead.
- `mode`: `'dual'` (regression + classification) or `'regression'`.

#### `train(trainer_config=None, bw_paths=None, target_reg_col='target_reg', rc_aug=True)`

Compiles training data from `df_regions` + BigWig signals, then trains.

#### `load_checkpoint(suffix='best')`

Loads `model_{suffix}.pt` from `results_dir/Models/`.

#### `load_finemo(hits_tsv_path, finemo_h5_path=None, score_col='hit_coefficient', score_threshold=None, similarity_threshold=None, auto_threshold_percentile=None)`

Reads Fi-NeMo `hits.tsv`, maps `motif_name → TF` via H5 attributes, converts
coordinates to genomic absolute using `df_regions`, and writes `motif_locs.csv`
+ `non_motif_locs.csv`.

- `score_col`: which Fi-NeMo score column to use. Default `'hit_coefficient'`.
- `auto_threshold_percentile`: if set (e.g. 50), auto-compute `score_threshold`
  from the score distribution. Recommended because Fi-NeMo scores are float.
- `similarity_threshold`: minimum `hit_similarity` (cosine). None = no filter.
- `finemo_h5_path`: H5 motif database for TF annotation. If None, motif_name
  is used as-is (no TF mapping).

#### `run_isa(isa_config, start_from='single_isa')`

Runs the 4-stage ISA pipeline. Supports `start_from` checkpointing.

`isa_config` keys:
```python
{
    'tracks': [0],                  # prediction track indices
    'null_percentile': 80,          # null distribution percentile threshold
    'receptive_field': 255,         # max motif pair distance (default: model.rf)
    'num_regions_per_batch': 200,   # batch size for region processing
    'pred_batch_size': 1024,        # GPU inference batch size
    'min_count': 10,                # min samples for cooperativity scoring
    'q_val_thresh': 0.1,            # FDR q-value threshold
}
```

#### `report()`

Generates ~17 plots per track: ISA distributions, cooperativity heatmaps,
TF family analysis, PPI validation, cell-type specificity, etc.

---

### `adapter.finemo_io`

Standalone functions for Fi-NeMo I/O (usable without EpQuickStart):

| Function | Purpose |
|----------|---------|
| `load_finemo_hits(path)` | Load raw 15-column hits.tsv |
| `load_motif_annotation(h5_path)` | Build {motif_name: {MC_ID, TF_Name}} from H5 |
| `annotate_hits(df, h5_path)` | Add MC_ID + TF_Name columns to hits |
| `prepare_region_map(df_regions)` | Build {peak_id: {chrom, start, end, region}} |
| `hits_to_motif_locs(df, region_map, ...)` | Convert hits → deepISA motif_locs format |
| `compute_non_motif_regions(df, region_map)` | Subtract motifs from regions |
| `suggest_score_threshold(df, percentile)` | Auto-suggest score threshold |
| `load_finemo_scan(hits_path, region_df, ...)` | High-level: full pipeline I/O |

---

## Fi-NeMo hits.tsv Format (15 columns)

| Column | Type | Description |
|--------|------|-------------|
| `chr` | str | Chromosome (`NA` if no peaks provided) |
| `start` | int | Trimmed motif start (0-based half-open) |
| `end` | int | Trimmed motif end (0-based half-open) |
| `start_untrimmed` | int | Untrimmed start |
| `end_untrimmed` | int | Untrimmed end |
| `motif_name` | str | `pos_patterns.pattern_{N}` |
| `hit_coefficient` | float | Local contribution weight (**recommended score**) |
| `hit_coefficient_global` | float | Globally normalized coefficient |
| `hit_similarity` | float | Cosine similarity to motif CWM [0, 1] |
| `hit_correlation` | float | Pearson correlation |
| `hit_importance` | float | Fractional importance in sequence |
| `hit_importance_sq` | float | Squared importance |
| `strand` | str | `+` or `-` |
| `peak_name` | str | Peak name from BED |
| `peak_id` | int | 0-based peak index |

**Coordinate note**: If Fi-NeMo was called without `-p <peaks>`, `chr=NA` and
`start/end` are relative to the input region. Ep_ISA_NEW converts these to absolute
genomic coordinates using `df_regions`.

---

## Output Files

All outputs are written to `results_dir/`:

```
results/
├── Data/
│   ├── motif_locs.csv           Fi-NeMo hits → ISA format
│   ├── non_motif_locs.csv       Non-motif intervals (for null model)
│   ├── pred_orig.csv            Original predictions per region
│   ├── motif_single_isa.csv     Single motif ISA scores
│   ├── null_isa.csv             Null ISA distribution
│   ├── motif_combi_isa.csv      Pairwise ISA + interaction scores
│   ├── null_interaction.csv     Null interaction distribution
│   ├── tf_importance.csv        Per-TF importance (KS test)
│   ├── coop_tf_pair_t{N}.csv    TF pair cooperativity scores
│   └── coop_tf_t{N}.csv         Single TF cooperativity scores
├── Models/
│   ├── model_best.pt            Best checkpoint
│   └── model_config.json        Model architecture
└── Plots/
    └── *.png                    ~17 plots per track
```

### coop_tf_pair columns

| Column | Description |
|--------|-------------|
| `tf_pair` | `TF1\|TF2` |
| `coop_score` | Direction-weighted score [-1, +1] |
| `mw_p` / `mw_q` | Mann-Whitney U p-value / BH-FDR q-value |
| `cooperativity` | `Synergistic` / `Redundant` / `Intermediate` / `Independent` |
| `median_distance` | Typical motif pair distance (bp) |

---

## Dependencies

- Python >= 3.9
- PyTorch >= 2.0
- pandas, numpy, scipy, scikit-learn
- bioframe, biopython, pyBigWig
- matplotlib, seaborn
- h5py
- loguru, tqdm
- statsmodels (FDR correction)

---

## Citation

If you use Ep_ISA_NEW, please cite:
- deepISA: He, X. "deepISA: A unified Python framework for mapping transcription factor cooperativity using deep learning and in silico ablation."
- Fi-NeMo: Relevant Fi-NeMo citation.

## License

Same as deepISA.
