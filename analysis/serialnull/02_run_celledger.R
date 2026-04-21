#!/usr/bin/env Rscript

source(file.path("analysis", "serialnull", "serialnull_utils.R"))

repo_root <- resolve_repo_root()
local_rlib <- file.path(repo_root, ".Rlib")
if (dir.exists(local_rlib)) {
  .libPaths(c(normalizePath(local_rlib), .libPaths()))
}

suppressPackageStartupMessages({
  library(Matrix)
})
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE)
} else if (requireNamespace("CellEdgeR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(CellEdgeR))
} else {
  stop("Install CellEdgeR or run this script from the repo with devtools available.")
}

paths <- serialnull_paths(repo_root)

celledger_cache_dir <- file.path(paths$cache_dir, "celledger")
celledger_results_dir <- file.path(paths$results_dir, "celledger")
ensure_dir(celledger_cache_dir)
ensure_dir(celledger_results_dir)

run_id <- Sys.getenv("SERIALNULL_RUN_ID", unset = format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))
n_cores <- as.integer(Sys.getenv("SERIALNULL_N_CORES", unset = "4"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

recompute_graph <- identical(tolower(Sys.getenv("SERIALNULL_RECOMPUTE_GRAPH", unset = "false")), "true")
recompute_motifs <- identical(tolower(Sys.getenv("SERIALNULL_RECOMPUTE_MOTIFS", unset = "false")), "true")
recompute_fit <- identical(tolower(Sys.getenv("SERIALNULL_RECOMPUTE_FIT", unset = "false")), "true")
write_fit_cache_requested <- identical(tolower(Sys.getenv("SERIALNULL_WRITE_FIT_CACHE", unset = "false")), "true")
if (write_fit_cache_requested) {
  message("SERIALNULL_WRITE_FIT_CACHE is ignored; split_stats caches are used to avoid storing multi-GB fitted cellgraph objects.")
}
read_fit_cache <- FALSE
write_fit_cache <- FALSE

max_edge_len <- NA_real_
include_wedge <- TRUE
celledger_offset_version <- "volume_chung_lu_v2"
celledger_stats_version <- "volume_filter_v3"

is_current_celledger_motifs <- function(obj) {
  is.list(obj) &&
    is.list(obj$parameters) &&
    identical(obj$parameters$offset_version, celledger_offset_version) &&
    identical(obj$parameters$offset_modes, "volume") &&
    is.list(obj$offsets) &&
    "volume" %in% names(obj$offsets)
}

manifest <- run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "celledger",
  split = "all",
  step = "load_sample_manifest",
  fn = function() {
    if (!file.exists(paths$manifest_csv)) {
      stop("Missing manifest CSV: ", paths$manifest_csv, ". Run 01_prepare_manifest.R first.")
    }
    utils::read.csv(paths$manifest_csv, stringsAsFactors = FALSE)
  }
)

required_cols <- c("sample", "position", "group_odd_even", "group_first_half_second_half")
if (!all(required_cols %in% names(manifest))) {
  stop("Manifest is missing required columns: ", paste(setdiff(required_cols, names(manifest)), collapse = ", "))
}

run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "celledger",
  split = "all",
  step = "validate_manifest_against_data",
  fn = function() {
    expected_manifest <- build_sample_manifest(paths$data_dir, n_expected = 56L)
    expected_manifest <- expected_manifest[, required_cols, drop = FALSE]
    manifest_cmp <- manifest[, required_cols, drop = FALSE]
    expected_manifest$position <- as.integer(expected_manifest$position)
    manifest_cmp$position <- as.integer(manifest_cmp$position)
    expected_manifest$sample <- as.character(expected_manifest$sample)
    manifest_cmp$sample <- as.character(manifest_cmp$sample)
    expected_manifest$group_odd_even <- as.character(expected_manifest$group_odd_even)
    manifest_cmp$group_odd_even <- as.character(manifest_cmp$group_odd_even)
    expected_manifest$group_first_half_second_half <- as.character(expected_manifest$group_first_half_second_half)
    manifest_cmp$group_first_half_second_half <- as.character(manifest_cmp$group_first_half_second_half)
    if (!identical(expected_manifest, manifest_cmp)) {
      stop(
        "Manifest does not match current data ordering/splits. ",
        "Re-run 01_prepare_manifest.R."
      )
    }
    invisible(NULL)
  }
)

cellgraph_motifs <- NULL
cellgraph_motifs_ready <- FALSE
ensure_cellgraph_motifs <- function() {
  if (cellgraph_motifs_ready) {
    return(cellgraph_motifs)
  }

  cellgraph_rds <- file.path(celledger_cache_dir, "cellgraph_base.rds")
  motif_rds <- file.path(celledger_cache_dir, "cellgraph_with_motifs.rds")

  if (file.exists(motif_rds) && !recompute_motifs) {
    cached_motifs <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = "all",
      step = "load_cached_motifs",
      fn = function() readRDS(motif_rds)
    )
    if (is_current_celledger_motifs(cached_motifs)) {
      cellgraph_motifs <<- cached_motifs
      cellgraph_motifs_ready <<- TRUE
      return(cellgraph_motifs)
    }
    cached_version <- if (is.list(cached_motifs$parameters) && !is.null(cached_motifs$parameters$offset_version)) {
      cached_motifs$parameters$offset_version
    } else {
      "<missing>"
    }
    message("Ignoring stale motif cache with offset_version: ", cached_version)
    rm(cached_motifs)
    gc(verbose = FALSE)
  }

  if (file.exists(cellgraph_rds) && !recompute_graph) {
    cellgraph <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = "all",
      step = "load_cached_cellgraph",
      fn = function() readRDS(cellgraph_rds)
    )
  } else {
    cells_by_sample <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = "all",
      step = "read_h5_slides",
      fn = function() {
        samples <- manifest$sample
        files <- file.path(paths$data_dir, paste0(samples, ".h5"))
        missing <- files[!file.exists(files)]
        if (length(missing)) {
          stop("Missing .h5 files for manifest samples: ", paste(basename(missing), collapse = ", "))
        }
        out <- lapply(files, read_serial_h5)
        stats::setNames(out, samples)
      }
    )
    cellgraph <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = "all",
      step = "build_cellgraph_delaunay",
      fn = function() {
        obj <- build_cell_graphs(cells_by_sample, n_cores = n_cores, verbose = TRUE)
        saveRDS(obj, cellgraph_rds)
        obj
      }
    )
    rm(cells_by_sample)
    gc(verbose = FALSE)
  }

  cellgraph_motifs <<- run_timed_step(
    timings_csv = paths$timings_csv,
    run_id = run_id,
    engine = "celledger",
    split = "all",
    step = "count_motifs_graphs",
    fn = function() {
      obj <- count_motifs_graphs(
        cellgraph,
        max_edge_len = max_edge_len,
        include_wedge = include_wedge,
        n_cores = n_cores,
        verbose = TRUE
      )
      saveRDS(obj, motif_rds)
      obj
    }
  )
  rm(cellgraph)
  gc(verbose = FALSE)

  cellgraph_motifs_ready <<- TRUE
  cellgraph_motifs
}

split_defs <- list(
  odd_even = list(
    split_col = "group_odd_even",
    expected_levels = c("odd", "even")
  ),
  first_half_second_half = list(
    split_col = "group_first_half_second_half",
    expected_levels = c("first_half", "second_half")
  )
)

all_pvals <- list()
all_uniformity <- list()

for (split_name in names(split_defs)) {
  split_col <- split_defs[[split_name]]$split_col
  expected_levels <- split_defs[[split_name]]$expected_levels

  sample_df <- sample_df_from_manifest(manifest, split_col = split_col)
  level_vals <- levels(sample_df$group)
  if (!identical(level_vals, expected_levels)) {
    stop(
      "Unexpected levels for ", split_name, ": ",
      paste(level_vals, collapse = ", "),
      " (expected ", paste(expected_levels, collapse = ", "), ")"
    )
  }
  if (sum(sample_df$group == expected_levels[1]) != sum(sample_df$group == expected_levels[2])) {
    stop("Split is not balanced for ", split_name)
  }

  fit_rds <- file.path(celledger_cache_dir, paste0("fit_", celledger_stats_version, "_", split_name, ".rds"))
  split_stats_rds <- file.path(celledger_cache_dir, paste0("split_stats_", celledger_stats_version, "_", split_name, ".rds"))

  split_pvals <- NULL
  split_summary <- NULL

  if (file.exists(split_stats_rds) && !recompute_fit) {
    split_stats <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = split_name,
      step = "load_cached_split_stats",
      fn = function() readRDS(split_stats_rds)
    )
    if (is.list(split_stats) && is.data.frame(split_stats$split_pvals) && is.data.frame(split_stats$split_summary)) {
      split_pvals <- split_stats$split_pvals
      split_summary <- split_stats$split_summary
    }
  }

  if (is.null(split_pvals) || is.null(split_summary)) {
    res_volume <- NULL
    used_cached_fit <- FALSE

    if (read_fit_cache && file.exists(fit_rds) && !recompute_fit) {
      fit_obj <- tryCatch(
        run_timed_step(
          timings_csv = paths$timings_csv,
          run_id = run_id,
          engine = "celledger",
          split = split_name,
          step = "load_cached_fit",
          fn = function() readRDS(fit_rds)
        ),
        error = function(e) {
          message(
            "Cached fit load failed for split ", split_name, ": ", conditionMessage(e),
            " ; falling back to recompute."
          )
          NULL
        }
      )
      if (!is.null(fit_obj) && !is.null(fit_obj$res_volume)) {
        res_volume <- fit_obj$res_volume
        used_cached_fit <- TRUE
      }
      rm(fit_obj)
      gc(verbose = FALSE)
    }

    if (!used_cached_fit) {
      cellgraph_motifs_obj <- ensure_cellgraph_motifs()

      res_volume <- run_timed_step(
        timings_csv = paths$timings_csv,
        run_id = run_id,
        engine = "celledger",
        split = split_name,
        step = "motif_edger_volume",
        fn = function() {
          motif_edger(
            cellgraph = cellgraph_motifs_obj,
            sample_df = sample_df,
            design_formula = "~ group",
            verbose = TRUE
          )
        }
      )

      if (write_fit_cache) {
        run_timed_step(
          timings_csv = paths$timings_csv,
          run_id = run_id,
          engine = "celledger",
          split = split_name,
          step = "write_fit_cache",
          fn = function() {
            saveRDS(
              list(
                split_name = split_name,
                split_col = split_col,
                group_levels = expected_levels,
                offset_version = celledger_offset_version,
                res_volume = res_volume
              ),
              fit_rds
            )
            invisible(NULL)
          }
        )
      }
    }

    split_pvals <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = split_name,
      step = "collect_pvalues",
      fn = function() {
        edges_tbl <- top_edges(res_volume, coef = NULL, model = "full")
        motifs2_tbl <- top_motifs2(res_volume, coef = NULL, model = "full")

        out <- rbind(
          collect_pvals(edges_tbl, "volume_edges", split_name),
          collect_pvals(motifs2_tbl, "volume_motifs2", split_name)
        )

        out <- out[is.finite(out$p_value), , drop = FALSE]
        out
      }
    )

    split_summary <- run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = split_name,
      step = "uniformity_summary",
      fn = function() {
        if (!nrow(split_pvals)) {
          return(data.frame(
            split = character(),
            analysis = character(),
            n = numeric(),
            mean_p = numeric(),
            median_p = numeric(),
            frac_p_lt_0_05 = numeric(),
            ks_pvalue = numeric(),
            binom_pvalue = numeric(),
            stringsAsFactors = FALSE
          ))
        }
        by_analysis <- split(split_pvals$p_value, split_pvals$analysis)
        parts <- lapply(names(by_analysis), function(analysis_name) {
          out <- uniformity_stats(by_analysis[[analysis_name]])
          out$analysis <- analysis_name
          out
        })
        do.call(rbind, parts)
      }
    )

    run_timed_step(
      timings_csv = paths$timings_csv,
      run_id = run_id,
      engine = "celledger",
      split = split_name,
      step = "write_split_stats_cache",
      fn = function() {
        saveRDS(
          list(
            split_name = split_name,
            split_col = split_col,
            group_levels = expected_levels,
            split_pvals = split_pvals,
            split_summary = split_summary
          ),
          split_stats_rds
        )
        invisible(NULL)
      }
    )
    rm(res_volume)
    gc(verbose = FALSE)
  }

  if ("run_id" %in% names(split_pvals)) {
    split_pvals$run_id <- NULL
  }

  if (nrow(split_summary)) {
    split_summary$split <- split_name
    split_summary <- split_summary[, c("split", "analysis", "n", "mean_p", "median_p", "frac_p_lt_0_05", "ks_pvalue", "binom_pvalue")]
  }

  split_pvals$run_id <- run_id
  all_pvals[[split_name]] <- split_pvals
  all_uniformity[[split_name]] <- split_summary
}

pvals_out <- do.call(rbind, all_pvals)
uniformity_out <- do.call(rbind, all_uniformity)

manifest_used_csv <- file.path(celledger_results_dir, "sample_manifest_used.csv")
pvals_csv <- file.path(celledger_results_dir, "pvalues.csv")
uniformity_csv <- file.path(celledger_results_dir, "uniformity_summary.csv")

run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "celledger",
  split = "all",
  step = "write_results",
  fn = function() {
    utils::write.csv(manifest, manifest_used_csv, row.names = FALSE, quote = TRUE)
    utils::write.csv(pvals_out, pvals_csv, row.names = FALSE, quote = TRUE)
    utils::write.csv(uniformity_out, uniformity_csv, row.names = FALSE, quote = TRUE)
    invisible(NULL)
  }
)

message("CellEdgeR results written to: ", celledger_results_dir)
message("p-values rows: ", nrow(pvals_out))
message("uniformity rows: ", nrow(uniformity_out))
