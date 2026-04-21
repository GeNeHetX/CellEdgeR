resolve_repo_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(start, dirname(start), dirname(dirname(start))))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && dir.exists(file.path(cand, "analysis"))) {
      return(cand)
    }
  }
  stop("Could not resolve repository root from: ", start)
}

serialnull_paths <- function(repo_root) {
  base_dir <- file.path(repo_root, "analysis", "serialnull")
  list(
    repo_root = repo_root,
    base_dir = base_dir,
    data_dir = file.path(repo_root, "analysis", "data", "raw", "serialNull"),
    config_dir = file.path(base_dir, "config"),
    cache_dir = file.path(base_dir, "cache"),
    results_dir = file.path(base_dir, "results"),
    manifest_csv = file.path(base_dir, "config", "sample_splits.csv"),
    timings_csv = file.path(base_dir, "results", "timings.csv")
  )
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

append_timing_row <- function(timings_csv, row) {
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  ensure_dir(dirname(timings_csv))
  if (!file.exists(timings_csv)) {
    utils::write.table(
      row,
      file = timings_csv,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = TRUE,
      append = FALSE
    )
  } else {
    utils::write.table(
      row,
      file = timings_csv,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      quote = TRUE,
      append = TRUE
    )
  }
}

run_timed_step <- function(timings_csv, run_id, engine, split, step, fn) {
  started_at <- Sys.time()
  status <- "ok"
  details <- ""
  value <- NULL

  value <- tryCatch(
    fn(),
    error = function(e) {
      status <<- "error"
      details <<- conditionMessage(e)
      NULL
    }
  )

  finished_at <- Sys.time()
  elapsed_sec <- as.numeric(difftime(finished_at, started_at, units = "secs"))

  append_timing_row(
    timings_csv = timings_csv,
    row = list(
      run_id = run_id,
      engine = engine,
      split = split,
      step = step,
      status = status,
      started_at = format(started_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      elapsed_sec = sprintf("%.6f", elapsed_sec),
      details = details
    )
  )

  if (!identical(status, "ok")) {
    stop("Step failed [", engine, " / ", split, " / ", step, "]: ", details)
  }

  value
}

decode_class_vector <- function(x) {
  if (is.factor(x)) return(as.character(x))

  if (is.character(x)) {
    if (length(dim(x)) > 1L) {
      out <- apply(x, 1L, paste0, collapse = "")
      return(trimws(out))
    }
    return(as.character(x))
  }

  if (is.list(x)) {
    return(as.character(unlist(x, use.names = FALSE)))
  }

  if (is.matrix(x) && (is.integer(x) || is.numeric(x))) {
    out <- apply(x, 1L, function(v) {
      bytes <- as.integer(v)
      bytes <- bytes[bytes > 0L]
      if (!length(bytes)) return("")
      rawToChar(as.raw(bytes))
    })
    return(out)
  }

  as.character(x)
}

read_serial_h5 <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package hdf5r is required to read raw serialnull .h5 slides.")
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  x <- as.numeric(h5[["x"]][])
  y <- as.numeric(h5[["y"]][])
  cls <- decode_class_vector(h5[["class"]][])

  if (length(x) != length(y) || length(x) != length(cls)) {
    stop(
      "Length mismatch in ", basename(path),
      ": x=", length(x), " y=", length(y), " class=", length(cls)
    )
  }

  data.frame(
    x = x,
    y = y,
    label = cls,
    stringsAsFactors = FALSE
  )
}

build_sample_manifest <- function(data_dir, n_expected = 56L) {
  h5_files <- sort(list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE))
  if (length(h5_files) != n_expected) {
    stop("Expected exactly ", n_expected, " .h5 files in ", data_dir, "; found ", length(h5_files), ".")
  }

  sample_names <- tools::file_path_sans_ext(basename(h5_files))
  if (anyDuplicated(sample_names)) {
    stop("Duplicate sample names after removing .h5 extension.")
  }

  n <- length(sample_names)
  if ((n %% 2L) != 0L) {
    stop("Number of samples must be even for odd/even split; found ", n, ".")
  }

  odd_even <- ifelse((seq_len(n) %% 2L) == 1L, "odd", "even")
  half_split <- ifelse(seq_len(n) <= (n / 2L), "first_half", "second_half")

  manifest <- data.frame(
    sample = sample_names,
    position = seq_len(n),
    group_odd_even = odd_even,
    group_first_half_second_half = half_split,
    stringsAsFactors = FALSE
  )

  if (sum(manifest$group_odd_even == "odd") != sum(manifest$group_odd_even == "even")) {
    stop("Odd/even split is unbalanced.")
  }
  if (sum(manifest$group_first_half_second_half == "first_half") !=
      sum(manifest$group_first_half_second_half == "second_half")) {
    stop("First-half/second-half split is unbalanced.")
  }

  manifest
}

sample_df_from_manifest <- function(manifest, split_col) {
  if (!split_col %in% names(manifest)) {
    stop("Split column not found in manifest: ", split_col)
  }

  group_vals <- as.character(manifest[[split_col]])
  levels_order <- unique(group_vals)

  data.frame(
    group = factor(group_vals, levels = levels_order),
    row.names = manifest$sample,
    stringsAsFactors = FALSE
  )
}

collect_pvals <- function(tbl, analysis, split_name) {
  if (is.null(tbl) || !nrow(tbl)) {
    return(data.frame(
      split = character(),
      analysis = character(),
      motif = character(),
      motif_type = character(),
      p_value = numeric(),
      fdr = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  fdr <- if ("FDR" %in% names(tbl)) {
    as.numeric(tbl$FDR)
  } else {
    stats::p.adjust(as.numeric(tbl$PValue), method = "BH")
  }
  motif_type <- if ("motif_type" %in% names(tbl)) as.character(tbl$motif_type) else NA_character_

  data.frame(
    split = split_name,
    analysis = analysis,
    motif = tbl$motif,
    motif_type = motif_type,
    p_value = as.numeric(tbl$PValue),
    fdr = fdr,
    stringsAsFactors = FALSE
  )
}

uniformity_stats <- function(p) {
  p <- p[is.finite(p)]
  n <- length(p)

  if (n < 2L) {
    return(data.frame(
      n = n,
      mean_p = NA_real_,
      median_p = NA_real_,
      frac_p_lt_0_05 = NA_real_,
      ks_pvalue = NA_real_,
      binom_pvalue = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ks <- stats::ks.test(p, "punif")
  bt <- stats::binom.test(sum(p < 0.05), n, p = 0.05)

  data.frame(
    n = n,
    mean_p = mean(p),
    median_p = stats::median(p),
    frac_p_lt_0_05 = mean(p < 0.05),
    ks_pvalue = ks$p.value,
    binom_pvalue = bt$p.value,
    stringsAsFactors = FALSE
  )
}
