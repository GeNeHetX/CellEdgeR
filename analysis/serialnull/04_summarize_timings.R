#!/usr/bin/env Rscript

source(file.path("analysis", "serialnull", "serialnull_utils.R"))

repo_root <- resolve_repo_root()
paths <- serialnull_paths(repo_root)

if (!file.exists(paths$timings_csv)) {
  stop("Timing log not found: ", paths$timings_csv)
}

timings <- utils::read.csv(paths$timings_csv, stringsAsFactors = FALSE)
if (!nrow(timings)) {
  stop("Timing log is empty: ", paths$timings_csv)
}

run_id <- Sys.getenv("SERIALNULL_RUN_ID", unset = "")
if (!nzchar(run_id)) {
  run_id <- tail(timings$run_id, 1)
}

timings_run <- timings[timings$run_id == run_id, , drop = FALSE]
if (!nrow(timings_run)) {
  stop("No timing rows found for run_id: ", run_id)
}

timings_run$elapsed_sec <- as.numeric(timings_run$elapsed_sec)

summary_by_step <- aggregate(
  elapsed_sec ~ run_id + engine + split + step + status,
  data = timings_run,
  FUN = sum
)
summary_by_step <- summary_by_step[order(summary_by_step$engine, summary_by_step$split, summary_by_step$step), ]

totals <- aggregate(
  elapsed_sec ~ run_id + engine + split,
  data = timings_run[timings_run$status == "ok", , drop = FALSE],
  FUN = sum
)
totals <- totals[order(totals$engine, totals$split), ]

latest_step_csv <- file.path(paths$results_dir, "timings_summary_latest.csv")
latest_total_csv <- file.path(paths$results_dir, "timings_totals_latest.csv")
history_total_csv <- file.path(paths$results_dir, "timings_totals_history.csv")

utils::write.csv(summary_by_step, latest_step_csv, row.names = FALSE, quote = TRUE)
utils::write.csv(totals, latest_total_csv, row.names = FALSE, quote = TRUE)

if (!file.exists(history_total_csv)) {
  utils::write.csv(totals, history_total_csv, row.names = FALSE, quote = TRUE)
} else {
  utils::write.table(
    totals,
    file = history_total_csv,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE,
    append = TRUE
  )
}

message("Timing summary written: ", latest_step_csv)
message("Timing totals written: ", latest_total_csv)
message("Timing totals appended: ", history_total_csv)
