# serialNull Script Pipeline

This directory contains a clean, script-based serialNull workflow for both CellEdgeR and GraphCompass.

## What it does

- Builds one shared sample manifest from `analysis/data/raw/serialNull/*.h5`.
- Defines two groupings from the same ordered sample list:
  - `group_odd_even` (null expectation)
  - `group_first_half_second_half` (non-null expectation)
- Runs CellEdgeR and GraphCompass analyses using the same manifest.
- Measures and appends step timings for every run.

## Main files

- `01_prepare_manifest.R`: writes `config/sample_splits.csv`.
- `02_run_celledger.R`: runs CellEdgeR with cached graph/motif stages.
- `03_run_graphcompass.py`: runs GraphCompass from the same manifest.
- `04_summarize_timings.R`: writes run timing summaries.
- `run_serialnull_pipeline.sh`: runs all steps in sequence.
- `serialnull_utils.R`: shared R helpers (I/O, timing, decoding, stats).

## Outputs

- Shared manifest:
  - `config/sample_splits.csv`
  - `config/sample_df_odd_even.csv`
  - `config/sample_df_first_half_second_half.csv`
- Timings:
  - `results/timings.csv` (append-only)
  - `results/timings_summary_latest.csv`
  - `results/timings_totals_latest.csv`
  - `results/timings_totals_history.csv` (append-only)
- CellEdgeR:
  - `results/celledger/pvalues.csv`
  - `results/celledger/uniformity_summary.csv`
  - `results/celledger/sample_manifest_used.csv`
  - `cache/celledger/split_stats_*.rds` (compact per-split cache used for reruns)
- GraphCompass:
  - `results/graphcompass/pvalues.csv`
  - `results/graphcompass/uniformity_summary.csv`
  - `results/graphcompass/sample_manifest_used.csv`
  - `cache/graphcompass/adata_serialnull.h5ad` (cached AnnData for faster reruns)
  - `cache/graphcompass/adata_manifest_used.csv` (manifest guard for AnnData cache)

## Run

```bash
bash analysis/run_serialnull_main.sh
```

Equivalent direct entrypoint:

```bash
bash analysis/serialnull/run_serialnull_pipeline.sh
```

Optional environment variables:

- `SERIALNULL_RUN_ID` (default: UTC timestamp)
- `SERIALNULL_N_CORES` (default: `4`)
- `SERIALNULL_RECOMPUTE_GRAPH=true|false`
- `SERIALNULL_RECOMPUTE_MOTIFS=true|false`
- `SERIALNULL_RECOMPUTE_FIT=true|false`
- `SERIALNULL_WRITE_FIT_CACHE=true|false` (default `false`; enables large `fit_*.rds` caches)
- `SERIALNULL_RECOMPUTE_GRAPHCOMPASS=true|false` (default `false`; recomputes AnnData + GraphCompass metric caches)
- `SERIALNULL_GRAPHCOMPASS_METHODS=portrait,filtration,wl` (default `portrait,filtration,wl`; diffusion is skipped by default)
- `SERIALNULL_GRAPHCOMPASS_N_JOBS=4` (default: `SERIALNULL_N_CORES` or `4`)
- `SERIALNULL_GRAPHCOMPASS_MAX_CELLS_PER_SAMPLE=250000` (optional per-sample downsampling cap for GraphCompass)
- `SERIALNULL_GRAPHCOMPASS_SAMPLE_SEED=1` (default `1`; used when downsampling)
- `PYTHON_BIN=/path/to/python`
- `R_BIN=Rscript`
