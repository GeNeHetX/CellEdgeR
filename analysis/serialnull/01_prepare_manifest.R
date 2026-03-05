#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
})

source(file.path("analysis", "serialnull", "serialnull_utils.R"))

repo_root <- resolve_repo_root()
paths <- serialnull_paths(repo_root)

ensure_dir(paths$config_dir)
ensure_dir(paths$results_dir)

run_id <- Sys.getenv("SERIALNULL_RUN_ID", unset = format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))

manifest <- run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "shared",
  split = "all",
  step = "prepare_sample_manifest",
  fn = function() build_sample_manifest(paths$data_dir, n_expected = 56L)
)

run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "shared",
  split = "all",
  step = "write_sample_manifest",
  fn = function() {
    utils::write.csv(manifest, paths$manifest_csv, row.names = FALSE, quote = TRUE)
    invisible(NULL)
  }
)

odd_even_csv <- file.path(paths$config_dir, "sample_df_odd_even.csv")
half_csv <- file.path(paths$config_dir, "sample_df_first_half_second_half.csv")

run_timed_step(
  timings_csv = paths$timings_csv,
  run_id = run_id,
  engine = "shared",
  split = "all",
  step = "write_split_sample_dfs",
  fn = function() {
    odd_even_df <- data.frame(
      sample = manifest$sample,
      group = manifest$group_odd_even,
      stringsAsFactors = FALSE
    )
    half_df <- data.frame(
      sample = manifest$sample,
      group = manifest$group_first_half_second_half,
      stringsAsFactors = FALSE
    )
    utils::write.csv(odd_even_df, odd_even_csv, row.names = FALSE, quote = TRUE)
    utils::write.csv(half_df, half_csv, row.names = FALSE, quote = TRUE)
    invisible(NULL)
  }
)

message("Wrote sample manifest: ", paths$manifest_csv)
message("Wrote split sample_df files: ", odd_even_csv, " and ", half_csv)
message("Rows: ", nrow(manifest))
message("Odd/even counts: ", paste(table(manifest$group_odd_even), collapse = "/"))
message("Half split counts: ", paste(table(manifest$group_first_half_second_half), collapse = "/"))
