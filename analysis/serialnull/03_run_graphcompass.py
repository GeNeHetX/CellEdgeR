#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import sys
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional


def _configure_third_party_warnings() -> None:
    warnings.filterwarnings(
        "ignore",
        message="The legacy Dask DataFrame implementation is deprecated.*",
        category=FutureWarning,
        module=r"dask\.dataframe",
    )
    warnings.filterwarnings(
        "ignore",
        message="pkg_resources is deprecated as an API.*",
        category=UserWarning,
        module=r"xarray_schema",
    )
    warnings.filterwarnings(
        "ignore",
        message="nopython is set for njit and is ignored",
        category=RuntimeWarning,
        module=r"numba\.core\.decorators",
    )
    warnings.filterwarnings(
        "ignore",
        message=r".*__version__.*anndata.*",
        category=FutureWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="Importing read_text from `anndata` is deprecated.*",
        category=FutureWarning,
        module=r"anndata",
    )


_configure_third_party_warnings()

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import binomtest, kstest, mannwhitneyu


REQUIRED_MANIFEST_COLUMNS = [
    "sample",
    "position",
    "group_odd_even",
    "group_first_half_second_half",
]

VALID_GRAPHCOMPASS_METRICS = ("portrait", "diffusion", "filtration", "wl")
DEFAULT_GRAPHCOMPASS_METRICS = ("portrait", "filtration", "wl")


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def resolve_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "DESCRIPTION").exists() and (candidate / "analysis").is_dir():
            return candidate
    raise RuntimeError(f"Could not resolve repository root from {start}")


def parse_positive_int_env(name: str, default: int) -> int:
    raw = os.environ.get(name, str(default)).strip()
    try:
        value = int(raw)
    except ValueError as exc:
        raise RuntimeError(f"Environment variable {name} must be a positive integer, got: {raw!r}") from exc
    if value < 1:
        raise RuntimeError(f"Environment variable {name} must be >= 1, got: {value}")
    return value


def parse_optional_positive_int_env(name: str) -> Optional[int]:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return None
    try:
        value = int(raw)
    except ValueError as exc:
        raise RuntimeError(f"Environment variable {name} must be a positive integer or empty, got: {raw!r}") from exc
    if value < 1:
        raise RuntimeError(f"Environment variable {name} must be >= 1, got: {value}")
    return value


def parse_metrics_env(name: str, default: tuple[str, ...]) -> list[str]:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return list(default)

    values = [x.strip().lower() for x in raw.split(",") if x.strip()]
    unknown = [x for x in values if x not in VALID_GRAPHCOMPASS_METRICS]
    if unknown:
        raise RuntimeError(
            f"Environment variable {name} contains unknown metric(s): {unknown}. "
            f"Allowed: {list(VALID_GRAPHCOMPASS_METRICS)}"
        )

    ordered_unique = []
    seen = set()
    for value in values:
        if value not in seen:
            ordered_unique.append(value)
            seen.add(value)

    if not ordered_unique:
        raise RuntimeError(f"Environment variable {name} is empty after parsing.")
    return ordered_unique


def build_cache_suffix(max_cells_per_sample: Optional[int], sample_seed: int) -> str:
    max_cells_tag = "allcells" if max_cells_per_sample is None else f"maxcells{max_cells_per_sample}"
    return f"{max_cells_tag}_seed{sample_seed}"


@dataclass
class Paths:
    repo_root: Path
    base_dir: Path
    data_dir: Path
    manifest_csv: Path
    cache_dir: Path
    results_dir: Path
    timings_csv: Path


class TimingLogger:
    def __init__(self, timings_csv: Path, run_id: str):
        self.timings_csv = timings_csv
        self.run_id = run_id
        self.timings_csv.parent.mkdir(parents=True, exist_ok=True)

    def append(self, engine: str, split: str, step: str, status: str, started_at: str, finished_at: str, elapsed_sec: float, details: str = "") -> None:
        write_header = not self.timings_csv.exists()
        with self.timings_csv.open("a", newline="") as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow([
                    "run_id",
                    "engine",
                    "split",
                    "step",
                    "status",
                    "started_at",
                    "finished_at",
                    "elapsed_sec",
                    "details",
                ])
            writer.writerow([
                self.run_id,
                engine,
                split,
                step,
                status,
                started_at,
                finished_at,
                f"{elapsed_sec:.6f}",
                details,
            ])

    def timed(self, engine: str, split: str, step: str, fn: Callable[[], object]) -> object:
        print(f"[{engine}] START split={split} step={step}", flush=True)
        started_at = utc_now_iso()
        started_ts = time.time()
        status = "ok"
        details = ""
        result = None
        try:
            result = fn()
            return result
        except Exception as exc:  # noqa: BLE001
            status = "error"
            details = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            finished_at = utc_now_iso()
            elapsed = time.time() - started_ts
            self.append(engine, split, step, status, started_at, finished_at, elapsed, details)
            print(f"[{engine}] END   split={split} step={step} status={status} elapsed={elapsed:.2f}s", flush=True)


def decode_class_vector(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr)

    if a.dtype.kind in {"S", "a", "U"}:
        a_str = np.char.decode(a, "utf-8", errors="ignore") if a.dtype.kind in {"S", "a"} else a.astype(str)
        if a_str.ndim > 1:
            return np.asarray(["".join(row.tolist()).strip() for row in a_str], dtype=str)
        return a_str.astype(str)

    if a.dtype == object:
        out = []
        for value in a:
            if isinstance(value, (bytes, bytearray, np.bytes_)):
                out.append(value.decode("utf-8", errors="ignore"))
            else:
                out.append(str(value))
        return np.asarray(out, dtype=str)

    if a.ndim == 2 and a.dtype.kind in {"i", "u", "f"}:
        out = []
        for row in a:
            byte_vals = [int(v) for v in row if int(v) > 0]
            out.append(bytes(byte_vals).decode("utf-8", errors="ignore"))
        return np.asarray(out, dtype=str)

    return a.astype(str)


def read_serial_h5(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        x = np.asarray(f["x"][:], dtype=float)
        y = np.asarray(f["y"][:], dtype=float)
        cls = decode_class_vector(f["class"][:])

    if not (len(x) == len(y) == len(cls)):
        raise ValueError(f"Length mismatch in {path.name}: x={len(x)}, y={len(y)}, class={len(cls)}")

    return pd.DataFrame({"x": x, "y": y, "label": cls})


def mw_pvalue(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def pair_class(sample_a: str, sample_b: str, sample_to_group: dict[str, str], group_a: str, group_b: str) -> str:
    ga = sample_to_group[sample_a]
    gb = sample_to_group[sample_b]
    if ga == group_a and gb == group_a:
        return f"within_{group_a}"
    if ga == group_b and gb == group_b:
        return f"within_{group_b}"
    return "between"


def pvals_from_pairwise(df: pd.DataFrame, value_col: str, method: str, sample_to_group: dict[str, str], group_a: str, group_b: str) -> pd.DataFrame:
    sub = df.copy()
    sub["pair_class"] = [
        pair_class(a, b, sample_to_group, group_a, group_b)
        for a, b in zip(sub["sample_a"], sub["sample_b"])
    ]

    records = []
    for cell_type, chunk in sub.groupby("cell_type"):
        wa = chunk.loc[chunk["pair_class"] == f"within_{group_a}", value_col].to_numpy()
        wb = chunk.loc[chunk["pair_class"] == f"within_{group_b}", value_col].to_numpy()
        between = chunk.loc[chunk["pair_class"] == "between", value_col].to_numpy()
        within_pooled = np.concatenate([wa, wb]) if (len(wa) + len(wb)) else np.array([])

        records.append(
            {
                "strategy": f"{method}__within_{group_a}_vs_within_{group_b}",
                "feature": str(cell_type),
                "p_value": mw_pvalue(wa, wb),
            }
        )
        records.append(
            {
                "strategy": f"{method}__between_vs_within_pooled",
                "feature": str(cell_type),
                "p_value": mw_pvalue(between, within_pooled),
            }
        )

    return pd.DataFrame(records)


def pvals_from_filtration(df: pd.DataFrame, group_a: str, group_b: str) -> pd.DataFrame:
    records = []
    for (cell_type, weight), chunk in df.groupby(["cell_type", "weight"]):
        a = chunk.loc[chunk["condition"] == group_a, "value"].to_numpy()
        b = chunk.loc[chunk["condition"] == group_b, "value"].to_numpy()
        records.append(
            {
                "strategy": f"filtration_curves__{group_a}_vs_{group_b}",
                "feature": f"{cell_type}@w={float(weight):.5g}",
                "p_value": mw_pvalue(a, b),
            }
        )
    return pd.DataFrame(records)


def pvals_from_wl(df: pd.DataFrame, sample_to_group: dict[str, str], group_a: str, group_b: str) -> pd.DataFrame:
    sub = df.copy()
    sub["pair_class"] = [
        pair_class(a, b, sample_to_group, group_a, group_b)
        for a, b in zip(sub["sample_a"], sub["sample_b"])
    ]

    wa = sub.loc[sub["pair_class"] == f"within_{group_a}", "distance"].to_numpy()
    wb = sub.loc[sub["pair_class"] == f"within_{group_b}", "distance"].to_numpy()
    between = sub.loc[sub["pair_class"] == "between", "distance"].to_numpy()
    within_pooled = np.concatenate([wa, wb]) if (len(wa) + len(wb)) else np.array([])

    return pd.DataFrame(
        [
            {
                "strategy": f"wlkernel__within_{group_a}_vs_within_{group_b}",
                "feature": "global_graph",
                "p_value": mw_pvalue(wa, wb),
            },
            {
                "strategy": "wlkernel__between_vs_within_pooled",
                "feature": "global_graph",
                "p_value": mw_pvalue(between, within_pooled),
            },
        ]
    )


def uniformity_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for strategy, chunk in df.groupby("strategy"):
        p = chunk["p_value"].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        n = len(p)
        if n < 2:
            rows.append(
                {
                    "strategy": strategy,
                    "n": n,
                    "mean_p": np.nan,
                    "median_p": np.nan,
                    "frac_p_lt_0_05": np.nan,
                    "ks_pvalue": np.nan,
                    "binom_pvalue": np.nan,
                }
            )
            continue

        ks = kstest(p, "uniform")
        bt = binomtest(int(np.sum(p < 0.05)), n=n, p=0.05)
        rows.append(
            {
                "strategy": strategy,
                "n": n,
                "mean_p": float(np.mean(p)),
                "median_p": float(np.median(p)),
                "frac_p_lt_0_05": float(np.mean(p < 0.05)),
                "ks_pvalue": float(ks.pvalue),
                "binom_pvalue": float(bt.pvalue),
            }
        )

    return pd.DataFrame(rows)


def manifest_matches_cache(manifest: pd.DataFrame, manifest_csv: Path) -> bool:
    if not manifest_csv.exists():
        return False
    try:
        cached = pd.read_csv(manifest_csv)
    except Exception:  # noqa: BLE001
        return False
    if not all(c in cached.columns for c in REQUIRED_MANIFEST_COLUMNS):
        return False
    cached = cached[REQUIRED_MANIFEST_COLUMNS].copy()
    cached["sample"] = cached["sample"].astype(str)
    manifest_std = manifest[REQUIRED_MANIFEST_COLUMNS].copy()
    manifest_std["sample"] = manifest_std["sample"].astype(str)
    return cached.reset_index(drop=True).equals(manifest_std.reset_index(drop=True))


def _write_adata_cache(adata: ad.AnnData, adata_h5ad: Path, manifest_csv: Path, manifest: pd.DataFrame) -> None:
    adata_h5ad.parent.mkdir(parents=True, exist_ok=True)
    manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(adata_h5ad, compression="lzf")
    manifest.to_csv(manifest_csv, index=False)


def main() -> None:
    repo_root = resolve_repo_root(Path(__file__).resolve().parent)

    graphcompass_src = repo_root / "analysis" / "external" / "graphcompass" / "src"
    if graphcompass_src.exists() and str(graphcompass_src) not in sys.path:
        sys.path.insert(0, str(graphcompass_src))

    import graphcompass as gc  # noqa: PLC0415

    paths = Paths(
        repo_root=repo_root,
        base_dir=repo_root / "analysis" / "serialnull",
        data_dir=repo_root / "analysis" / "data" / "raw" / "serialNull",
        manifest_csv=repo_root / "analysis" / "serialnull" / "config" / "sample_splits.csv",
        cache_dir=repo_root / "analysis" / "serialnull" / "cache" / "graphcompass",
        results_dir=repo_root / "analysis" / "serialnull" / "results" / "graphcompass",
        timings_csv=repo_root / "analysis" / "serialnull" / "results" / "timings.csv",
    )

    paths.cache_dir.mkdir(parents=True, exist_ok=True)
    paths.results_dir.mkdir(parents=True, exist_ok=True)

    run_id = os.environ.get("SERIALNULL_RUN_ID", datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    recompute = os.environ.get("SERIALNULL_RECOMPUTE_GRAPHCOMPASS", "false").lower() == "true"
    enabled_metrics = parse_metrics_env("SERIALNULL_GRAPHCOMPASS_METHODS", DEFAULT_GRAPHCOMPASS_METRICS)
    default_n_jobs = parse_positive_int_env("SERIALNULL_N_CORES", default=4)
    n_jobs = parse_positive_int_env("SERIALNULL_GRAPHCOMPASS_N_JOBS", default=default_n_jobs)
    max_cells_per_sample = parse_optional_positive_int_env("SERIALNULL_GRAPHCOMPASS_MAX_CELLS_PER_SAMPLE")
    sample_seed = parse_positive_int_env("SERIALNULL_GRAPHCOMPASS_SAMPLE_SEED", default=1)
    cache_suffix = build_cache_suffix(max_cells_per_sample=max_cells_per_sample, sample_seed=sample_seed)

    # graphcompass.tl.distance uses this env variable to set joblib parallelism.
    os.environ["GRAPHCOMPASS_N_JOBS"] = str(n_jobs)

    print(
        f"[graphcompass] config metrics={enabled_metrics} n_jobs={n_jobs} "
        f"max_cells_per_sample={max_cells_per_sample} sample_seed={sample_seed}",
        flush=True,
    )

    logger = TimingLogger(paths.timings_csv, run_id)

    manifest = logger.timed(
        engine="graphcompass",
        split="all",
        step="load_sample_manifest",
        fn=lambda: pd.read_csv(paths.manifest_csv),
    )

    missing_cols = [c for c in REQUIRED_MANIFEST_COLUMNS if c not in manifest.columns]
    if missing_cols:
        raise RuntimeError(f"Manifest missing required columns: {missing_cols}")

    manifest = manifest[REQUIRED_MANIFEST_COLUMNS].copy()
    manifest["sample"] = manifest["sample"].astype(str)

    samples = manifest["sample"].tolist()
    if len(samples) != 56:
        raise RuntimeError(f"Expected 56 samples in manifest, found {len(samples)}")

    if len(set(samples)) != len(samples):
        raise RuntimeError("Duplicate sample names in manifest")

    logger.timed(
        engine="graphcompass",
        split="all",
        step="validate_manifest_splits",
        fn=lambda: _validate_manifest_groups(manifest),
    )

    adata_cache_h5ad = paths.cache_dir / f"adata_serialnull_{cache_suffix}.h5ad"
    adata_manifest_used_csv = paths.cache_dir / f"adata_manifest_used_{cache_suffix}.csv"
    pairwise_manifest_used_csv = paths.cache_dir / f"sample_manifest_used_{cache_suffix}.csv"

    legacy_adata_cache_h5ad = paths.cache_dir / "adata_serialnull.h5ad"
    legacy_adata_manifest_used_csv = paths.cache_dir / "adata_manifest_used.csv"
    legacy_pairwise_manifest_used_csv = paths.cache_dir / "sample_manifest_used.csv"

    # Preserve compatibility with existing cache names for the default config.
    if max_cells_per_sample is None and sample_seed == 1:
        if (not adata_cache_h5ad.exists()) and legacy_adata_cache_h5ad.exists():
            adata_cache_h5ad = legacy_adata_cache_h5ad
        if (not adata_manifest_used_csv.exists()) and legacy_adata_manifest_used_csv.exists():
            adata_manifest_used_csv = legacy_adata_manifest_used_csv
        if (not pairwise_manifest_used_csv.exists()) and legacy_pairwise_manifest_used_csv.exists():
            pairwise_manifest_used_csv = legacy_pairwise_manifest_used_csv

    use_adata_cache = (
        (not recompute)
        and adata_cache_h5ad.exists()
        and manifest_matches_cache(manifest, adata_manifest_used_csv)
    )

    adata = None
    if use_adata_cache:
        try:
            adata = logger.timed(
                engine="graphcompass",
                split="all",
                step="load_cached_anndata",
                fn=lambda: ad.read_h5ad(adata_cache_h5ad),
            )
        except Exception as exc:  # noqa: BLE001
            print(
                f"[graphcompass] Cached AnnData load failed ({exc}); rebuilding cache.",
                flush=True,
            )

    if adata is None:
        adata = logger.timed(
            engine="graphcompass",
            split="all",
            step="build_anndata_from_manifest",
            fn=lambda: _build_adata(
                paths.data_dir,
                samples,
                max_cells_per_sample=max_cells_per_sample,
                sample_seed=sample_seed,
            ),
        )
        logger.timed(
            engine="graphcompass",
            split="all",
            step="write_anndata_cache",
            fn=lambda: _write_adata_cache(adata, adata_cache_h5ad, adata_manifest_used_csv, manifest),
        )

    odd_even_map = dict(zip(manifest["sample"], manifest["group_odd_even"]))
    adata.obs["condition"] = pd.Categorical(adata.obs["sample"].map(odd_even_map), categories=["odd", "even"])

    portrait_csv = paths.cache_dir / f"pairwise_portrait_{cache_suffix}.csv"
    diffusion_csv = paths.cache_dir / f"pairwise_diffusion_{cache_suffix}.csv"
    filtration_csv = paths.cache_dir / f"filtration_long_{cache_suffix}.csv"
    wl_pairs_csv = paths.cache_dir / f"wl_pairwise_{cache_suffix}.csv"

    legacy_portrait_csv = paths.cache_dir / "pairwise_portrait.csv"
    legacy_diffusion_csv = paths.cache_dir / "pairwise_diffusion.csv"
    legacy_filtration_csv = paths.cache_dir / "filtration_long.csv"
    legacy_wl_pairs_csv = paths.cache_dir / "wl_pairwise.csv"

    if max_cells_per_sample is None and sample_seed == 1:
        if (not portrait_csv.exists()) and legacy_portrait_csv.exists():
            portrait_csv = legacy_portrait_csv
        if (not diffusion_csv.exists()) and legacy_diffusion_csv.exists():
            diffusion_csv = legacy_diffusion_csv
        if (not filtration_csv.exists()) and legacy_filtration_csv.exists():
            filtration_csv = legacy_filtration_csv
        if (not wl_pairs_csv.exists()) and legacy_wl_pairs_csv.exists():
            wl_pairs_csv = legacy_wl_pairs_csv

    pairwise_manifest_valid = (not recompute) and manifest_matches_cache(manifest, pairwise_manifest_used_csv)

    pairwise_portrait: Optional[pd.DataFrame] = None
    pairwise_diffusion: Optional[pd.DataFrame] = None
    filtration_long: Optional[pd.DataFrame] = None
    wl_pairwise: Optional[pd.DataFrame] = None

    wrote_metric_cache = False
    spatial_graphs_ready = False

    if "portrait" in enabled_metrics:
        if pairwise_manifest_valid and portrait_csv.exists():
            pairwise_portrait = logger.timed(
                "graphcompass",
                "all",
                "load_cached_pairwise_portrait",
                lambda: pd.read_csv(portrait_csv),
            )
        else:
            pairwise_portrait = logger.timed(
                "graphcompass",
                "all",
                "compute_pairwise_portrait",
                lambda: _compute_pairwise(
                    adata,
                    gc,
                    method="portrait",
                    compute_spatial_graphs=not spatial_graphs_ready,
                ),
            )
            spatial_graphs_ready = True
            pairwise_portrait.to_csv(portrait_csv, index=False)
            wrote_metric_cache = True

    if "diffusion" in enabled_metrics:
        if pairwise_manifest_valid and diffusion_csv.exists():
            pairwise_diffusion = logger.timed(
                "graphcompass",
                "all",
                "load_cached_pairwise_diffusion",
                lambda: pd.read_csv(diffusion_csv),
            )
        else:
            pairwise_diffusion = logger.timed(
                "graphcompass",
                "all",
                "compute_pairwise_diffusion",
                lambda: _compute_pairwise(
                    adata,
                    gc,
                    method="diffusion",
                    compute_spatial_graphs=not spatial_graphs_ready,
                ),
            )
            spatial_graphs_ready = True
            pairwise_diffusion.to_csv(diffusion_csv, index=False)
            wrote_metric_cache = True

    if "filtration" in enabled_metrics:
        if pairwise_manifest_valid and filtration_csv.exists():
            filtration_long = logger.timed(
                "graphcompass",
                "all",
                "load_cached_filtration",
                lambda: pd.read_csv(filtration_csv),
            )
        else:
            filtration_long = logger.timed(
                "graphcompass",
                "all",
                "compute_filtration_curves",
                lambda: _compute_filtration_long(
                    adata,
                    gc,
                    odd_even_map,
                    compute_spatial_graphs=not spatial_graphs_ready,
                ),
            )
            spatial_graphs_ready = True
            filtration_long.to_csv(filtration_csv, index=False)
            wrote_metric_cache = True

    if "wl" in enabled_metrics:
        if pairwise_manifest_valid and wl_pairs_csv.exists():
            wl_pairwise = logger.timed(
                "graphcompass",
                "all",
                "load_cached_wl",
                lambda: pd.read_csv(wl_pairs_csv),
            )
        else:
            wl_pairwise = logger.timed(
                "graphcompass",
                "all",
                "compute_wl_pairwise",
                lambda: _compute_wl_pairwise(
                    adata,
                    gc,
                    compute_spatial_graphs=not spatial_graphs_ready,
                ),
            )
            spatial_graphs_ready = True
            wl_pairwise.to_csv(wl_pairs_csv, index=False)
            wrote_metric_cache = True

    if wrote_metric_cache:
        manifest.to_csv(pairwise_manifest_used_csv, index=False)

    split_defs = [
        ("odd_even", "group_odd_even", "odd", "even"),
        ("first_half_second_half", "group_first_half_second_half", "first_half", "second_half"),
    ]

    pval_frames = []
    summary_frames = []

    for split_name, col, group_a, group_b in split_defs:
        split_pvals = logger.timed(
            "graphcompass",
            split_name,
            "compute_split_pvalues",
            lambda split_col=col, ga=group_a, gb=group_b: _compute_split_pvalues(
                manifest,
                split_col,
                ga,
                gb,
                pairwise_portrait=pairwise_portrait,
                pairwise_diffusion=pairwise_diffusion,
                filtration_long=filtration_long,
                wl_pairwise=wl_pairwise,
            ),
        )
        split_pvals["split"] = split_name
        split_pvals["run_id"] = run_id
        pval_frames.append(split_pvals)

        split_summary = logger.timed(
            "graphcompass",
            split_name,
            "compute_uniformity_summary",
            lambda df=split_pvals: uniformity_summary(df),
        )
        split_summary["split"] = split_name
        split_summary["run_id"] = run_id
        summary_frames.append(split_summary)

    pvals_out = pd.concat(pval_frames, axis=0, ignore_index=True)
    pvals_out = pvals_out[np.isfinite(pvals_out["p_value"])].copy()

    summary_out = pd.concat(summary_frames, axis=0, ignore_index=True)

    logger.timed(
        "graphcompass",
        "all",
        "write_results",
        lambda: _write_graphcompass_results(paths.results_dir, manifest, pvals_out, summary_out),
    )

    print(f"GraphCompass results written to: {paths.results_dir}")
    print(f"p-values rows: {len(pvals_out)}")
    print(f"uniformity rows: {len(summary_out)}")


def _validate_manifest_groups(manifest: pd.DataFrame) -> None:
    oe = manifest["group_odd_even"].astype(str).tolist()
    hs = manifest["group_first_half_second_half"].astype(str).tolist()

    expected_oe = ["odd" if (i % 2 == 1) else "even" for i in range(1, len(oe) + 1)]
    expected_hs = ["first_half" if i <= (len(hs) // 2) else "second_half" for i in range(1, len(hs) + 1)]

    if oe != expected_oe:
        raise RuntimeError("Manifest odd/even labels do not match alternating assignment")
    if hs != expected_hs:
        raise RuntimeError("Manifest first/second-half labels do not match expected split")


def _build_adata(
    data_dir: Path,
    samples: list[str],
    max_cells_per_sample: Optional[int] = None,
    sample_seed: int = 1,
) -> ad.AnnData:
    obs_parts = []
    xy_parts = []

    for sample_idx, sample in enumerate(samples):
        path = data_dir / f"{sample}.h5"
        if not path.exists():
            raise FileNotFoundError(f"Missing .h5 for sample: {sample}")

        df = read_serial_h5(path)
        n = len(df)
        if max_cells_per_sample is not None and n > max_cells_per_sample:
            print(
                f"[graphcompass] downsampling sample={sample} n={n} -> {max_cells_per_sample}",
                flush=True,
            )
            rng = np.random.default_rng(sample_seed + sample_idx)
            keep_idx = np.sort(rng.choice(n, size=max_cells_per_sample, replace=False))
            df = df.iloc[keep_idx].reset_index(drop=True)
            n = len(df)
        idx = [f"{sample}_{i}" for i in range(n)]

        obs = pd.DataFrame(
            {
                "sample": sample,
                "cell_type": df["label"].astype(str).values,
            },
            index=idx,
        )
        obs_parts.append(obs)
        xy_parts.append(df[["x", "y"]].to_numpy())

    obs_all = pd.concat(obs_parts, axis=0)
    obs_all["sample"] = pd.Categorical(obs_all["sample"], categories=samples, ordered=True)
    obs_all["cell_type"] = obs_all["cell_type"].astype("category")

    x_onehot = pd.get_dummies(obs_all["cell_type"], dtype=float)
    adata = ad.AnnData(X=sparse.csr_matrix(x_onehot.to_numpy()), obs=obs_all)
    adata.obsm["spatial"] = np.vstack(xy_parts)
    return adata


def _compute_pairwise(adata: ad.AnnData, gc, method: str, compute_spatial_graphs: bool) -> pd.DataFrame:
    gc.tl.distance.compare_conditions(
        adata,
        library_key="sample",
        cluster_key="cell_type",
        method=method,
        compute_spatial_graphs=compute_spatial_graphs,
        kwargs_spatial_neighbors={"coord_type": "generic", "spatial_key": "spatial", "delaunay": True},
    )
    return adata.uns["pairwise_similarities"].copy()


def _compute_filtration_long(
    adata: ad.AnnData,
    gc,
    sample_to_group: dict[str, str],
    compute_spatial_graphs: bool,
) -> pd.DataFrame:
    gc.tl.filtration_curves.compare_conditions(
        adata,
        library_key="sample",
        cluster_key="cell_type",
        condition_key="condition",
        compute_spatial_graphs=compute_spatial_graphs,
    )

    curves = adata.uns["filtration_curves"]["curves"]
    long_parts = []
    for sample, df in curves.items():
        tmp = df.copy()
        tmp["sample"] = sample
        tmp["condition"] = sample_to_group[sample]
        value_cols = [
            c for c in tmp.columns if c not in {"graph_label", "weight", "sample", "condition"}
        ]
        tmp = tmp.melt(
            id_vars=["sample", "condition", "weight"],
            value_vars=value_cols,
            var_name="cell_type",
            value_name="value",
        )
        long_parts.append(tmp)

    return pd.concat(long_parts, axis=0, ignore_index=True)


def _compute_wl_pairwise(adata: ad.AnnData, gc, compute_spatial_graphs: bool) -> pd.DataFrame:
    gc.tl.wlkernel.compare_conditions(
        adata,
        library_key="sample",
        cluster_key="cell_type",
        compute_spatial_graphs=compute_spatial_graphs,
    )
    wl_mat = adata.uns["wl_kernel"]["wasserstein_distance"].copy()

    rows = []
    samples = wl_mat.index.tolist()
    for i, sample_a in enumerate(samples):
        for j in range(i + 1, len(samples)):
            sample_b = samples[j]
            rows.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "distance": float(wl_mat.loc[sample_a, sample_b]),
                }
            )
    return pd.DataFrame(rows)


def _compute_split_pvalues(
    manifest: pd.DataFrame,
    split_col: str,
    group_a: str,
    group_b: str,
    pairwise_portrait: Optional[pd.DataFrame] = None,
    pairwise_diffusion: Optional[pd.DataFrame] = None,
    filtration_long: Optional[pd.DataFrame] = None,
    wl_pairwise: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    sample_to_group = dict(zip(manifest["sample"], manifest[split_col]))

    pval_parts = []
    if pairwise_portrait is not None:
        pval_parts.append(
            pvals_from_pairwise(
                pairwise_portrait,
                "similarity_score",
                "distance_portrait",
                sample_to_group,
                group_a,
                group_b,
            )
        )
    if pairwise_diffusion is not None:
        pval_parts.append(
            pvals_from_pairwise(
                pairwise_diffusion,
                "similarity_score",
                "distance_diffusion",
                sample_to_group,
                group_a,
                group_b,
            )
        )
    if filtration_long is not None:
        filt = filtration_long.copy()
        filt["condition"] = filt["sample"].map(sample_to_group)
        pval_parts.append(pvals_from_filtration(filt, group_a, group_b))
    if wl_pairwise is not None:
        pval_parts.append(pvals_from_wl(wl_pairwise, sample_to_group, group_a, group_b))

    if not pval_parts:
        raise RuntimeError("No GraphCompass metrics were enabled; no p-values can be computed.")

    pvals = pd.concat(pval_parts, axis=0, ignore_index=True)
    return pvals


def _write_graphcompass_results(results_dir: Path, manifest: pd.DataFrame, pvals_out: pd.DataFrame, summary_out: pd.DataFrame) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(results_dir / "sample_manifest_used.csv", index=False)
    pvals_out.to_csv(results_dir / "pvalues.csv", index=False)
    summary_out.to_csv(results_dir / "uniformity_summary.csv", index=False)


if __name__ == "__main__":
    main()
