# helper utilities -------------------------------------------------------------

standardize_sample_df <- function(sample_df, sample_name) {
  df <- as.data.frame(sample_df)
  if (ncol(df) < 3L) {
    stop("Sample '", sample_name, "' must have at least three columns (x, y, labels).")
  }
  coords <- df[, seq_len(2), drop = FALSE]
  numeric_ok <- vapply(coords, is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop("Sample '", sample_name, "' requires numeric coordinates in columns 1 and 2.")
  }
  label_column <- df[[3]]
  if (!is.atomic(label_column) || length(label_column) != nrow(df)) {
    stop("Sample '", sample_name, "' third column must be a vector of labels.")
  }
  data.frame(
    x = coords[[1]],
    y = coords[[2]],
    label = label_column,
    stringsAsFactors = FALSE
  )
}

validate_cells_by_sample <- function(cells_by_sample, max_labels = 100) {
  if (!is.list(cells_by_sample) || length(cells_by_sample) == 0) {
    stop("cells_by_sample must be a non-empty *named* list of data frames.")
  }
  samples <- names(cells_by_sample)
  if (is.null(samples) || anyNA(samples) || any(samples == "")) {
    stop("cells_by_sample must have non-empty, non-NA names for each sample.")
  }
  if (anyDuplicated(samples)) stop("cells_by_sample sample names must be unique.")
  cells_standard <- stats::setNames(
    lapply(samples, function(s) {
      df <- cells_by_sample[[s]]
      if (!is.data.frame(df)) stop("Sample '", s, "' is not a data.frame.")
      df_std <- standardize_sample_df(df, s)
      if (nrow(df_std) == 0) {
        warning("Sample '", s, "' has zero rows; it will be kept but contributes no motifs.")
      }
      if (any(!is.finite(df_std$x) | !is.finite(df_std$y))) {
        stop("Sample '", s, "' has non-finite coordinates.")
      }
      if (anyNA(df_std$label)) stop("Sample '", s, "' has NA labels.")
      df_std
    }),
    samples
  )
  all_labels <- unlist(lapply(cells_standard, function(d) as.character(d$label)), use.names = FALSE)
  if (length(unique(all_labels)) > max_labels) {
    stop("Too many unique labels (", length(unique(all_labels)), "); max allowed is ", max_labels, ".")
  }
  cells_standard
}

validate_graph_obj <- function(graph_obj) {
  if (!inherits(graph_obj, "cellEdgeR_graphs") && !inherits(graph_obj, "cellEdgeR_obj")) {
    stop("graph_obj must come from build_cell_graphs().")
  }
  if (is.null(graph_obj$samples) || !length(graph_obj$samples)) stop("graph_obj has no samples.")
  if (anyDuplicated(graph_obj$samples)) stop("graph_obj samples must be unique.")
  if (!is.list(graph_obj$per_sample)) stop("graph_obj$per_sample must be a list.")
  missing <- setdiff(graph_obj$samples, names(graph_obj$per_sample))
  if (length(missing)) stop("graph_obj missing per-sample entries: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

validate_motif_obj <- function(motif_obj, require_offsets = TRUE) {
  if (!is.list(motif_obj) || is.null(motif_obj$counts)) {
    stop("motif_obj must be the list returned by count_motifs_graphs().")
  }
  if (!inherits(motif_obj, "cellEdgeR_obj")) {
    stop("motif_obj must come from build_cell_graphs()/count_motifs_graphs().")
  }
  if (is.null(motif_obj$samples) || !length(motif_obj$samples)) {
    stop("motif_obj is missing sample identifiers.")
  }
  if (!is.list(motif_obj$counts)) stop("motif_obj$counts must be a list of sparse matrices.")
  if (is.null(motif_obj$exposure)) stop("motif_obj is missing exposure totals.")
  if (require_offsets && is.null(motif_obj$offsets)) {
    stop("motif_obj is missing offsets; rerun count_motifs_graphs() to compute them.")
  }
  invisible(TRUE)
}

empty_edge_result <- function() {
  list(edges = matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("from", "to"))), edge_len = numeric(0))
}

build_delaunay_edges <- function(df, verbose = FALSE) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("Package geometry is required. Install it with install.packages('geometry').")
  }
  xy <- as.matrix(df[, c("x", "y")])
  if (nrow(xy) == 0) return(empty_edge_result())
  keep_idx <- which(!duplicated(data.frame(x = xy[, 1], y = xy[, 2])))
  xy_unique <- xy[keep_idx, , drop = FALSE]
  if (NROW(xy_unique) < 2) return(empty_edge_result())
  tri <- tryCatch(
    geometry::delaunayn(xy_unique),
    error = function(e) stop("geometry::delaunayn failed: ", conditionMessage(e), call. = FALSE)
  )
  if (is.null(tri) || NROW(tri) == 0) return(empty_edge_result())
  tri <- as.matrix(tri)
  edges <- rbind(tri[, c(1, 2), drop = FALSE], tri[, c(2, 3), drop = FALSE], tri[, c(1, 3), drop = FALSE])
  edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  edges <- unique(edges)
  edges <- matrix(keep_idx[edges], ncol = 2, dimnames = list(NULL, c("from", "to")))
  storage.mode(edges) <- "integer"
  len <- sqrt(rowSums((xy[edges[, 1], , drop = FALSE] - xy[edges[, 2], , drop = FALSE])^2))
  list(edges = edges, edge_len = len)
}

format_label_triplets <- function(keys, lab_levels, prefix = NULL) {
  if (!length(keys)) return(character(0))
  keys <- gsub("\\|", "_", keys)
  ids_mat <- do.call(rbind, strsplit(keys, "_"))
  lab_mat <- matrix(lab_levels[as.integer(ids_mat)], ncol = ncol(ids_mat))
  base <- apply(lab_mat, 1, function(v) paste(v, collapse = "_"))
  if (!is.null(prefix)) paste0(prefix, "_", base) else base
}

prefix_key <- function(labels, prefix) paste0(prefix, "_", labels)

strip_prefix <- function(keys) sub("^[^_]+_", "", keys)

pair_key_vec <- function(u, v, prefix = "E") {
  idx <- u <= v
  out <- character(length(u))
  out[idx] <- paste(u[idx], v[idx], sep = "_")
  out[!idx] <- paste(v[!idx], u[!idx], sep = "_")
  paste0(prefix, "_", out)
}

split_pair_labels <- function(keys) {
  ab <- strsplit(strip_prefix(keys), "_")
  do.call(rbind, ab)
}

split_triplet_labels <- function(keys) {
  abc <- strsplit(strip_prefix(keys), "_")
  do.call(rbind, abc)
}

filter_edges_by_threshold <- function(edges, lengths, threshold) {
  if (is.null(edges) || length(edges) == 0) return(matrix(numeric(0), ncol = 2))
  if (!is.finite(threshold) || length(lengths) == 0) return(edges)
  keep <- lengths <= threshold
  if (!any(keep)) return(matrix(numeric(0), ncol = ncol(edges)))
  edges[keep, , drop = FALSE]
}

#' Build spatial graphs via Delaunay triangulation
#'
#' Compute the underlying Delaunay edges for each sample without pruning so that motif
#' counts can be recomputed with different thresholds or wedge options without rebuilding.
#'
#' @param cells_by_sample Named list of data frames where the first two columns are numeric coordinates
#'   and the third holds cell-type labels; column names are not required as long as the ordering is preserved.
#' @param n_cores Parallelism hint; values greater than 1 trigger `parallel::mclapply` on Unix-alikes, while 1 runs sequentially.
#' @param verbose Logical; print progress.
#' @return A cellgraph object (class `cellEdgeR_obj`) containing the sample names, global label set, and per-sample edges/labels. Motif slots (`counts`, `offsets`, etc.) are filled by [count_motifs_graphs()].
#' @export
build_cell_graphs <- function(
  cells_by_sample,
  n_cores = 1,
  verbose = TRUE,
  max_labels = 100
) {
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("Package geometry is required. Install it with install.packages('geometry').")
  }
  cells_standard <- validate_cells_by_sample(cells_by_sample, max_labels = max_labels)
  samples <- names(cells_standard)
  if (verbose) message("Samples: ", paste(samples, collapse = ", "))

  all_labels_chr <- unlist(
    lapply(cells_standard, function(d) as.character(d$label)),
    use.names = FALSE
  )
  lab_levels <- sort(unique(all_labels_chr))
  lab_to_id <- stats::setNames(seq_len(length(lab_levels)), lab_levels)

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    stop("n_cores must be a single positive number.")
  }
  run_parallel <- n_cores > 1 && .Platform$OS.type != "windows"
  if (verbose) {
    msg <- if (run_parallel) paste0(n_cores, " cores") else "sequential"
    message("Building Delaunay graphs (", msg, ")…")
  }

  per_sample <- if (run_parallel) {
    parallel::mclapply(samples, function(s, ...) {
      df <- cells_standard[[s]]
      labs_chr <- as.character(df$label)
      labs_id <- unname(lab_to_id[labs_chr])
      edges_info <- build_delaunay_edges(df, verbose = verbose)
      list(
        sample = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len,
        xy = as.matrix(df[, c("x", "y")])
      )
    }, mc.cores = n_cores)
  } else {
    lapply(samples, function(s, ...) {
      df <- cells_standard[[s]]
      labs_chr <- as.character(df$label)
      labs_id <- unname(lab_to_id[labs_chr])
      edges_info <- build_delaunay_edges(df, verbose = verbose)
      list(
        sample = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len,
        xy = as.matrix(df[, c("x", "y")])
      )
    })
  }
  names(per_sample) <- samples

  structure(
    list(
      samples = samples,
      label_levels = lab_levels,
      lab_to_id = lab_to_id,
      per_sample = per_sample,
      counts = NULL,
      exposure = NULL,
      offsets = NULL,
      norm_counts = NULL,
      meta = list(max_edge_len = NA_real_, include_wedges = FALSE, offset_pseudo = NA_real_, built_from = "cells")
    ),
    class = c("cellEdgeR_obj", "cellEdgeR_graphs")
  )
}

#' Count cell-type motifs on Delaunay graphs
#'
#' Build a Delaunay triangulation for each sample (when \code{graph_obj} is \code{NULL}), prune edges that exceed \code{max_edge_len},
#' and count cell-type singletons, unordered pairs (edges), unordered triplets (triangles), and optionally wedges.
#' Triangle counts are accelerated via the C++ helper exposed in \pkg{CellEdgeR}.
#'
#' @param cells_by_sample Named list of data frames where the first two columns are numeric coordinates and the third holds cell-type labels.
#'   Provide this argument when you do not already have a prebuilt graph object.
#' @param graph_obj Optional output of \code{build_cell_graphs()}; when supplied, the triangulation is reused and \code{n_cores} is ignored.
#' @param max_edge_len Numeric threshold; Delaunay edges longer than this are dropped. Set to \code{NA}, \code{NULL}, or \code{<= 0} to skip pruning.
#' @param n_cores Parallelism hint; values greater than 1 trigger \code{parallel::mclapply} (Unix-alike only) when building graphs from raw samples.
#' @param include_wedges Logical; if \code{TRUE}, returns open-triplet (wedge) counts alongside triangles.
#' @param verbose Logical; print progress.
#' @param offset_pseudo Small positive constant used inside offsets.
#' @param max_labels Maximum allowed number of unique labels (guards against malformed inputs).
#'
#' @return A cellgraph object (class `cellEdgeR_obj`) with the original graph info plus motif `counts`, `exposure`, `offsets`, `norm_counts`, and `relative_counts` (log2 residuals from intercept-only edgeR fits).
#' @examples
#' demo <- list(
#'   s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
#'   s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
#' )
#' graphs <- build_cell_graphs(demo, verbose = FALSE)
#' counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
#' counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedges = TRUE, verbose = FALSE)
#' str(counts_full$counts)
#' @rdname count_motifs_graphs
#' @export
count_motifs_graphs <- function(
  cells_by_sample = NULL,
  graph_obj = NULL,
  max_edge_len = NA_real_,
  n_cores = 1,
  include_wedges = FALSE,
  verbose = TRUE,
  offset_pseudo = 0.5,
  max_labels = 100
) {
  graphs <- NULL
  if (!is.null(graph_obj)) {
    validate_graph_obj(graph_obj)
    graphs <- graph_obj
  } else if (!is.null(cells_by_sample) && inherits(cells_by_sample, "cellEdgeR_graphs")) {
    validate_graph_obj(cells_by_sample)
    graphs <- cells_by_sample
    cells_by_sample <- NULL
  } else if (!is.null(cells_by_sample)) {
    graphs <- build_cell_graphs(cells_by_sample, n_cores = n_cores, verbose = verbose, max_labels = max_labels)
  } else {
    stop("Provide either cells_by_sample or graph_obj.")
  }

  counts_obj <- count_motifs_from_graphs(
    graphs = graphs,
    max_edge_len = max_edge_len,
    include_wedges = include_wedges,
    verbose = verbose
  )
  offsets <- compute_motif_offsets(counts_obj, offset_pseudo)
  norm_counts <- normalize_counts_simple(counts_obj$counts, offsets)

  graphs$counts <- counts_obj$counts
  graphs$exposure <- counts_obj$exposure
  graphs$offsets <- offsets
  graphs$norm_counts <- norm_counts
  graphs$relative_counts <- compute_relative_counts(graphs, offset_pseudo, verbose = verbose)
  graphs$meta$max_edge_len <- max_edge_len
  graphs$meta$include_wedges <- include_wedges
  graphs$meta$offset_pseudo <- offset_pseudo
  graphs
}

count_motifs_from_graphs <- function(graphs, max_edge_len, include_wedges, verbose) {
  if (!inherits(graphs, "cellEdgeR_graphs")) {
    stop("graphs must come from build_cell_graphs().")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package Matrix is required. Install it with install.packages('Matrix').")
  }
  samples <- graphs$samples
  lab_levels <- graphs$label_levels
  K <- length(lab_levels)

  prune_value <- if (length(max_edge_len)) max_edge_len[1] else NA_real_
  prune_threshold <- suppressWarnings(as.numeric(prune_value))
  prune_threshold <- if (is.null(prune_value) || is.na(prune_threshold) || prune_threshold <= 0) {
    Inf
  } else {
    prune_threshold
  }
  apply_prune <- is.finite(prune_threshold)
  if (verbose && !apply_prune) {
    message("Edge pruning disabled (max_edge_len is NA/NULL/<= 0).")
  }

  per_sample <- lapply(samples, function(s) {
    ps <- graphs$per_sample[[s]]
    edges <- filter_edges_by_threshold(ps$edges, ps$edge_len, prune_threshold)
    list(sample = s, n = ps$n, labels_id = ps$labels_id, labels_chr = ps$labels_chr, edges = edges)
  })
  names(per_sample) <- samples

  if (verbose) message("Counting size-1 (cells by label)…")
  sing_list <- lapply(per_sample, function(ps) {
    if (ps$n == 0) return(integer(0))
    tab <- tabulate(ps$labels_id, nbins = K)
    names(tab) <- prefix_key(lab_levels, "N")
    tab
  })
  Y1 <- do.call(cbind, sing_list)
  dimnames(Y1) <- list(prefix_key(lab_levels, "N"), samples)
  Y1 <- Matrix::Matrix(Y1, sparse = TRUE)

  if (verbose) message("Counting size-2 (edges by unordered label pair)…")
  pairs_by_sample <- vector("list", length(samples))
  names(pairs_by_sample) <- samples
  for (s in samples) {
    ps <- per_sample[[s]]
    e <- ps$edges
    if (length(e) == 0) {
      pairs_by_sample[[s]] <- integer(0)
      next
    }
    la <- ps$labels_chr[e[, 1]]
    lb <- ps$labels_chr[e[, 2]]
    key <- pair_key_vec(la, lb, prefix = "E")
    pairs_by_sample[[s]] <- sort(tapply(rep(1L, length(key)), key, sum))
  }
  pairs_all_keys <- unique(unlist(lapply(pairs_by_sample, names), use.names = FALSE))
  if (length(pairs_all_keys)) {
    total_pairs <- sum(vapply(pairs_by_sample, length, integer(1)))
    if (total_pairs == 0L) {
      Y2 <- Matrix::Matrix(0, nrow = length(pairs_all_keys), ncol = length(samples),
        sparse = TRUE, dimnames = list(pairs_all_keys, samples))
    } else {
      pair_rows <- integer(total_pairs)
      pair_cols <- integer(total_pairs)
      pair_vals <- integer(total_pairs)
      pos <- 0L
      key_to_idx <- stats::setNames(seq_along(pairs_all_keys), pairs_all_keys)
      for (col in seq_along(samples)) {
        tab <- pairs_by_sample[[samples[col]]]
        len <- length(tab)
        if (len == 0L) next
        idx_range <- seq_len(len) + pos
        pair_rows[idx_range] <- key_to_idx[names(tab)]
        pair_cols[idx_range] <- col
        pair_vals[idx_range] <- as.integer(tab)
        pos <- pos + len
      }
      Y2 <- Matrix::sparseMatrix(
        i = pair_rows, j = pair_cols, x = pair_vals,
        dims = c(length(pairs_all_keys), length(samples)),
        dimnames = list(pairs_all_keys, samples)
      )
    }
  } else {
    Y2 <- Matrix::Matrix(0, nrow = 0, ncol = length(samples), sparse = TRUE, dimnames = list(character(), samples))
  }

  if (verbose) {
    message("Counting size-3 (triangles by unordered label triplet) in C++…")
    if (include_wedges) message("Also collecting wedges (open triplets)…")
  }
  tris_counts_list <- vector("list", length(samples))
  names(tris_counts_list) <- samples
  tri_totals <- numeric(length(samples))
  names(tri_totals) <- samples
  wedge_counts_list <- vector("list", length(samples))
  names(wedge_counts_list) <- samples
  wedge_totals <- numeric(length(samples))
  names(wedge_totals) <- samples

  build_sparse_layer <- function(count_list, keys) {
    if (!length(keys)) {
      return(Matrix::Matrix(0, nrow = 0, ncol = length(samples), sparse = TRUE,
        dimnames = list(character(), samples)))
    }
    total <- sum(vapply(count_list, length, integer(1)))
    if (total == 0L) {
      return(Matrix::Matrix(0, nrow = length(keys), ncol = length(samples), sparse = TRUE,
        dimnames = list(keys, samples)))
    }
    rows <- integer(total)
    cols <- integer(total)
    vals <- integer(total)
    pos <- 0L
    key_to_idx <- stats::setNames(seq_along(keys), keys)
    for (col in seq_along(samples)) {
      tab <- count_list[[samples[col]]]
      len <- length(tab)
      if (len == 0L) next
      idx_range <- seq_len(len) + pos
      rows[idx_range] <- key_to_idx[names(tab)]
      cols[idx_range] <- col
      vals[idx_range] <- as.integer(tab)
      pos <- pos + len
    }
    Matrix::sparseMatrix(
      i = rows, j = cols, x = vals,
      dims = c(length(keys), length(samples)),
      dimnames = list(keys, samples)
    )
  }

  for (s in samples) {
    ps <- per_sample[[s]]
    if (nrow(ps$edges) == 0 || ps$n < 3) {
      tris_counts_list[[s]] <- integer(0)
      tri_totals[s] <- 0
      wedge_counts_list[[s]] <- integer(0)
      wedge_totals[s] <- 0
      next
    }
    out <- count_triangle_labels_cpp(
      n_nodes = ps$n,
      ei = ps$edges[, 1],
      ej = ps$edges[, 2],
      labels = ps$labels_id,
      label_ids = seq_len(K),
      count_wedges = include_wedges
    )
    keys <- as.character(out$keys)
    if (length(keys)) {
      keys_lab <- format_label_triplets(keys, lab_levels, prefix = "T")
      tris_counts_list[[s]] <- structure(as.integer(out$counts), names = keys_lab)
    } else {
      tris_counts_list[[s]] <- integer(0)
    }
    tri_totals[s] <- as.numeric(out$tri_total)
    if (include_wedges) {
      wedge_keys <- as.character(out$wedge_keys)
      if (length(wedge_keys)) {
        wedge_lab <- format_label_triplets(wedge_keys, lab_levels, prefix = "W")
        wedge_counts_list[[s]] <- structure(as.integer(out$wedge_counts), names = wedge_lab)
      } else {
        wedge_counts_list[[s]] <- integer(0)
      }
      wedge_totals[s] <- as.numeric(out$wedge_total)
    } else {
      wedge_counts_list[[s]] <- integer(0)
      wedge_totals[s] <- 0
    }
  }
  tris_all_keys <- unique(unlist(lapply(tris_counts_list, names), use.names = FALSE))
  Y3 <- build_sparse_layer(tris_counts_list, tris_all_keys)
  if (include_wedges) {
    wedge_all_keys <- unique(unlist(lapply(wedge_counts_list, names), use.names = FALSE))
    Yw <- build_sparse_layer(wedge_counts_list, wedge_all_keys)
  }

  n_edges <- vapply(per_sample, function(ps) nrow(ps$edges), 0L)
  names(n_edges) <- samples
  n_tris <- tri_totals

  if (verbose) {
    summary_msg <- paste0(
      "Counts ready: |labels|=", K,
      ", singles=", nrow(Y1),
      ", pairs=", nrow(Y2),
      ", triangles=", nrow(Y3)
    )
    if (include_wedges) summary_msg <- paste0(summary_msg, ", wedges=", nrow(Yw))
    message(summary_msg)
  }
  counts <- list(size1 = Y1, size2 = Y2, size3 = Y3)
  if (include_wedges) counts$wedges <- Yw
  exposure <- list(edges = n_edges, triangles = n_tris, cells = Matrix::colSums(Y1))
  if (include_wedges) exposure$wedges <- wedge_totals

  list(
    samples = samples,
    label_levels = lab_levels,
    counts = counts,
    exposure = exposure,
    meta = list(max_edge_len = max_edge_len, include_wedges = include_wedges)
  )
}

#' Build motif offsets for normalized counts and modeling
#'
#' @keywords internal
build_pair_offsets <- function(Y2, Y1, log_edges, log_cells, pseudo) {
  if (nrow(Y2) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y2)))
  }
  ab <- split_pair_labels(rownames(Y2))
  a <- ab[, 1]
  b <- ab[, 2]
  NA_ <- Matrix::Matrix(0, nrow = length(a), ncol = ncol(Y1), dimnames = list(rownames(Y2), colnames(Y1)))
  NB_ <- Matrix::Matrix(0, nrow = length(b), ncol = ncol(Y1), dimnames = list(rownames(Y2), colnames(Y1)))
  idx_a <- match(prefix_key(a, "N"), rownames(Y1))
  idx_b <- match(prefix_key(b, "N"), rownames(Y1))
  sel_a <- which(!is.na(idx_a))
  sel_b <- which(!is.na(idx_b))
  if (length(sel_a)) NA_[sel_a, ] <- as.matrix(Y1[idx_a[sel_a], , drop = FALSE])
  if (length(sel_b)) NB_[sel_b, ] <- as.matrix(Y1[idx_b[sel_b], , drop = FALSE])
  off <- log(pmax(NA_, pseudo)) + log(pmax(NB_, pseudo))
  off <- sweep(off, 2, log_edges[colnames(Y2)], FUN = "-")
  off <- sweep(off, 2, 2 * log_cells[colnames(Y2)], FUN = "-")
  dimnames(off) <- dimnames(Y2)
  off
}

#' @keywords internal
build_tri_offsets <- function(Y3, Y2, log_tris, pseudo) {
  if (nrow(Y3) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y3)))
  }
  abc <- split_triplet_labels(rownames(Y3))
  a <- abc[, 1]
  b <- abc[, 2]
  c <- abc[, 3]
  ABk <- pair_key_vec(a, b, prefix = "E")
  ACk <- pair_key_vec(a, c, prefix = "E")
  BCk <- pair_key_vec(b, c, prefix = "E")
  get_pair_mat <- function(keys) {
    idx <- match(keys, rownames(Y2))
    out <- matrix(0, nrow = length(keys), ncol = ncol(Y2), dimnames = list(keys, colnames(Y2)))
    sel <- which(!is.na(idx))
    if (length(sel)) out[sel, ] <- as.matrix(Y2[idx[sel], , drop = FALSE])
    out
  }
  AB <- get_pair_mat(ABk)
  AC <- get_pair_mat(ACk)
  BC <- get_pair_mat(BCk)
  off <- log(pmax(AB, pseudo)) + log(pmax(AC, pseudo)) + log(pmax(BC, pseudo))
  off <- sweep(off, 2, log_tris[colnames(Y3)], FUN = "-")
  dimnames(off) <- dimnames(Y3)
  off
}

build_wedge_offsets <- function(Yw, Y2, log_wedges, pseudo) {
  if (is.null(log_wedges) || nrow(Yw) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Yw)))
  }
  abc <- split_triplet_labels(rownames(Yw))
  a <- abc[, 1]
  b <- abc[, 2]
  c <- abc[, 3]
  ABk <- pair_key_vec(a, b, prefix = "E")
  ACk <- pair_key_vec(a, c, prefix = "E")
  BCk <- pair_key_vec(b, c, prefix = "E")
  get_pair_mat <- function(keys) {
    idx <- match(keys, rownames(Y2))
    out <- matrix(0, nrow = length(keys), ncol = ncol(Y2), dimnames = list(keys, colnames(Y2)))
    sel <- which(!is.na(idx))
    if (length(sel)) out[sel, ] <- as.matrix(Y2[idx[sel], , drop = FALSE])
    out
  }
  AB <- get_pair_mat(ABk)
  AC <- get_pair_mat(ACk)
  BC <- get_pair_mat(BCk)
  off <- log(pmax(AB, pseudo)) + log(pmax(AC, pseudo)) + log(pmax(BC, pseudo))
  off <- sweep(off, 2, log_wedges[colnames(Yw)], FUN = "-")
  dimnames(off) <- dimnames(Yw)
  off
}

compute_motif_offsets <- function(motif_obj, pseudo) {
  samples <- motif_obj$samples
  Y1 <- motif_obj$counts$size1
  Y2 <- motif_obj$counts$size2
  Y3 <- motif_obj$counts$size3
  Yw <- motif_obj$counts$wedges

  log_cells_vec <- log(pmax(as.numeric(motif_obj$exposure$cells), 1))
  names(log_cells_vec) <- samples
  log_cells <- matrix(
    log_cells_vec,
    nrow = max(1, nrow(Y1)),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(rownames(Y1), samples)
  )
  log_edges <- log(pmax(as.numeric(motif_obj$exposure$edges), 1))
  names(log_edges) <- samples
  log_tris <- log(pmax(as.numeric(motif_obj$exposure$triangles), 1))
  names(log_tris) <- samples
  log_wedges <- NULL
  if (!is.null(motif_obj$exposure$wedges)) {
    log_wedges <- log(pmax(as.numeric(motif_obj$exposure$wedges), 1))
    names(log_wedges) <- samples
  }

  offset2 <- build_pair_offsets(Y2, Y1, log_edges, log_cells_vec, pseudo)
  offset3 <- build_tri_offsets(Y3, Y2, log_tris, pseudo)
  offsetw <- NULL
  if (!is.null(Yw)) {
    offsetw <- build_wedge_offsets(Yw, Y2, log_wedges, pseudo)
  }
  list(
    size1 = Matrix::Matrix(log_cells, sparse = TRUE, dimnames = dimnames(Y1)),
    size2 = Matrix::Matrix(offset2, sparse = TRUE, dimnames = dimnames(Y2)),
    size3 = Matrix::Matrix(offset3, sparse = TRUE, dimnames = dimnames(Y3)),
    wedges = if (!is.null(offsetw)) Matrix::Matrix(offsetw, sparse = TRUE, dimnames = dimnames(Yw)) else NULL
  )
}

normalize_counts_simple <- function(counts, offsets) {
  norm_layer <- function(Y, off) {
    if (is.null(Y) || nrow(Y) == 0) {
      return(Matrix::Matrix(0, nrow = 0, ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    Matrix::Matrix(as.matrix(Y) / exp(as.matrix(off)), sparse = TRUE, dimnames = dimnames(Y))
  }
  list(
    size1 = norm_layer(counts$size1, offsets$size1),
    size2 = norm_layer(counts$size2, offsets$size2),
    size3 = norm_layer(counts$size3, offsets$size3),
    wedges = if (!is.null(counts$wedges) && !is.null(offsets$wedges)) norm_layer(counts$wedges, offsets$wedges) else NULL
  )
}

compute_relative_counts <- function(motif_obj, pseudo, verbose = FALSE, eps = 0.5) {
  samples <- motif_obj$samples
  offsets <- motif_obj$offsets
  counts <- motif_obj$counts
  design_null <- matrix(1, nrow = length(samples), ncol = 1)
  colnames(design_null) <- "Intercept"

  fit_layer <- function(Y, off) {
    if (is.null(Y) || nrow(Y) == 0) {
      return(Matrix::Matrix(0, nrow = 0, ncol = length(samples),
        sparse = TRUE, dimnames = list(NULL, samples)))
    }
    if (ncol(Y) < 2) {
      warning("Skipping model residuals for layer with <2 samples; returning zeros.")
      return(Matrix::Matrix(0, nrow = nrow(Y), ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    dge <- edgeR::DGEList(counts = Y)
    dge <- edgeR::calcNormFactors(dge)
    dge$offset <- as.matrix(off)
    dge <- tryCatch(edgeR::estimateDisp(dge, design = design_null), error = function(e) NULL)
    if (is.null(dge)) {
      warning("edgeR dispersion estimation failed for model residuals; returning zeros.")
      return(Matrix::Matrix(0, nrow = nrow(Y), ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    fit <- tryCatch(edgeR::glmQLFit(dge, design = design_null), error = function(e) NULL)
    if (is.null(fit)) {
      warning("edgeR fit failed for model residuals; returning zeros.")
      return(Matrix::Matrix(0, nrow = nrow(Y), ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    mu <- fit$fitted.values
    obs <- as.matrix(Y)
    Matrix::Matrix(log2(obs + eps) - log2(mu + eps), sparse = TRUE, dimnames = dimnames(Y))
  }

  list(
    size1 = fit_layer(counts$size1, offsets$size1),
    size2 = fit_layer(counts$size2, offsets$size2),
    size3 = fit_layer(counts$size3, offsets$size3),
    wedges = if (!is.null(counts$wedges) && !is.null(offsets$wedges)) fit_layer(counts$wedges, offsets$wedges) else NULL
  )
}


#' edgeR-normalized motif counts with optional null expectations
#'
#' Fits edgeR models per motif layer using supplied offsets and design, then returns normalized
#' counts or log-residual diagnostics relative to the fitted means. Useful when you want expectations
#' under a specific design (optionally removing a coefficient) rather than raw offsets alone.
#'
#' @param motif_obj Output of [count_motifs_graphs()].
#' @param pseudo Small positive constant added to counts in log/ratio computations.
#' @param return_log Logical; if `TRUE`, returns `log(obs + eps) - log(fitted + eps)`, otherwise `(obs + eps)/(fitted + eps)`.
#' @param eps Small stabilizer for ratios/logs.
#' @return List with normalized matrices (`size1`, `size2`, `size3`) 
#' @keywords internal
edgeRnorm_motif_counts  <- function(
  motif_obj,
  pseudo = 0.5,
  return_log = TRUE,
  eps = 0.5
) {
  .Deprecated("count_motifs_graphs()$norm_counts")
  offsets <- if (!is.null(motif_obj$offsets)) motif_obj$offsets else compute_motif_offsets(motif_obj, pseudo)
  base <- normalize_counts_simple(motif_obj$counts, offsets)
  if (!return_log) return(base)
  to_log <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0) return(mat)
    Matrix::Matrix(log2(as.matrix(mat) + eps), sparse = TRUE, dimnames = dimnames(mat))
  }
  list(
    size1 = to_log(base$size1),
    size2 = to_log(base$size2),
    size3 = to_log(base$size3),
    wedges = to_log(base$wedges)
  )
}

#' @rdname edgeRnorm_motif_counts
#' @keywords internal
normalize_motif_counts <- function(
  motif_obj,
  pseudo = 0.5,
  sample_df = NULL,
  design_formula = "~ 1",
  use_null = TRUE,
  return_log = FALSE,
  eps = 0.5,
  verbose = FALSE
) {
  .Deprecated("count_motifs_graphs()$norm_counts")
  if (is.null(motif_obj$norm_counts)) {
    offsets <- compute_motif_offsets(motif_obj, pseudo)
    return(normalize_counts_simple(motif_obj$counts, offsets))
  }
  motif_obj$norm_counts
}



#' Differential motif testing with edgeR and optional FDR control
#'
#' Fit edgeR QL GLMs for each motif size using offsets derived from expected counts,
#' apply a standard BH correction by default, and optionally run hierarchical DAG-based control
#' (the historical DAGGER workflow).
#'
#' @param motif_obj Output list from [count_motifs_graphs()].
#' @param sample_df Data frame with sample metadata; rownames must match sample IDs.
#' @param design_formula Formula string passed to `model.matrix`, e.g. `~ condition + batch`.
#' @param coef Target coefficient (name or column index) for the design matrix.
#' @param pseudo Small positive number to avoid zeros in log offsets.
#' @param alpha FDR threshold used when `fdr_method = "dagger"`.
#' @param verbose Logical; print progress.
#' @param fdr_method Character string selecting the multiple-testing correction strategy;
#'   `"bh"` (default) returns layerwise BH-adjusted p-values, `"dagger"` runs the
#'   DAGGER-style hierarchical correction.
#'
#' @return A list with `design`, `results` (data frames for each size),
#'   `fdr` metadata describing the method/threshold, and, when applicable, DAG nodes and edges.
#' @export
motif_edger <- function(
  motif_obj,
  sample_df,
  design_formula,
  coef = NULL,
  coef_variable = NULL,
  coef_level = NULL,
  pseudo = 0.5,
  alpha = 0.05,
  verbose = TRUE,
  fdr_method = c("bh", "dagger"),
  merge_triplets = FALSE
) {
  #"edgeR", "Matrix"
  
  validate_motif_obj(motif_obj, require_offsets = TRUE)
  samples <- motif_obj$samples
  if (is.null(rownames(sample_df))) stop("sample_df must have rownames = sample IDs")
  sample_df <- as.data.frame(sample_df[samples, , drop = FALSE])
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_df)
  resolve_coef <- function() {
    if (!is.null(coef)) {
      if (is.character(coef)) {
        idx <- match(coef, colnames(design))
        if (is.na(idx)) stop("coef name not found.")
        return(idx)
      }
      return(coef)
    }
    non_int <- which(colnames(design) != "(Intercept)")
    if (!is.null(coef_variable)) {
      target_cols <- grep(paste0("^", coef_variable), colnames(design))
      if (length(target_cols) == 0) stop("coef_variable not found in design.")
      if (!is.null(coef_level)) {
        target_name <- paste0(coef_variable, coef_level)
        idx <- match(target_name, colnames(design))
        if (!is.na(idx)) return(idx)
      }
      return(target_cols[1])
    }
    if (length(non_int) == 0) stop("No non-intercept term found; set `coef`.")
    non_int[1]
  }
  coef <- resolve_coef()

  fdr_method <- match.arg(fdr_method)

  Y1 <- motif_obj$counts$size1
  Y2 <- motif_obj$counts$size2
  Y3 <- motif_obj$counts$size3
  Yw <- motif_obj$counts$wedges

  offsets <- motif_obj$offsets
  if (is.null(offsets)) {
    pseudo_local <- if (!is.null(motif_obj$meta$offset_pseudo)) motif_obj$meta$offset_pseudo else pseudo
    offsets <- compute_motif_offsets(motif_obj, pseudo_local)
  }
  log_cells <- offsets$size1
  off2 <- offsets$size2
  off3 <- offsets$size3
  offw <- offsets$wedges

  # edgeR fits (per layer)
  fit_layer <- function(Y, off) {
    if (nrow(Y) == 0) return(NULL)
    dge <- edgeR::DGEList(counts = Y)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge$offset <- as.matrix(off) # keep offsets on the DGEList to avoid double-passing to edgeR
    design_rank <- qr(design)$rank
    resid_df <- ncol(Y) - design_rank
    if (resid_df <= 0) {
      if (verbose) message("Insufficient residual df (", resid_df, "); returning NA results for layer.")
      return(data.frame(
        motif = rownames(Y),
        logFC = rep(NA_real_, nrow(Y)),
        PValue = rep(NA_real_, nrow(Y)),
        FDR_BH = rep(NA_real_, nrow(Y)),
        row.names = NULL
      ))
    }

    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- tryCatch(
      edgeR::glmQLFit(dge, design = design),
      error = function(err) {
        if (grepl("NA dispersions not allowed", err$message, fixed = TRUE)) {
          if (verbose) message("edgeR GLM fit failed: ", err$message)
          return(NULL)
        }
        stop(err)
      }
    )
    if (is.null(fit)) {
      return(data.frame(
        motif = rownames(Y),
        logFC = rep(NA_real_, nrow(Y)),
        PValue = rep(NA_real_, nrow(Y)),
        FDR_BH = rep(NA_real_, nrow(Y)),
        row.names = NULL
      ))
    }
    tst <- edgeR::glmQLFTest(fit, coef = coef)
    tt <- edgeR::topTags(tst, n = Inf, sort.by = "none")$table
    data.frame(
      motif = rownames(Y),
      logFC = tt$logFC,
      PValue = tt$PValue,
      FDR_BH = stats::p.adjust(tt$PValue, method = "BH"),
      row.names = NULL
    )
  }

  merge_triplet_layer <- function(Y_tri, off_tri, Y_wedge, off_wedge) {
    counts_list <- list()
    offsets_list <- list()
    if (!is.null(Y_tri) && nrow(Y_tri)) {
      rownames(Y_tri) <- sub("^[^_]+_", "TW_", rownames(Y_tri))
      counts_list[[length(counts_list) + 1]] <- Y_tri
      offsets_list[[length(offsets_list) + 1]] <- off_tri
    }
    if (!is.null(Y_wedge) && nrow(Y_wedge)) {
      rownames(Y_wedge) <- sub("^[^_]+_", "TW_", rownames(Y_wedge))
      counts_list[[length(counts_list) + 1]] <- Y_wedge
      offsets_list[[length(offsets_list) + 1]] <- off_wedge
    }
    if (!length(counts_list)) return(list(counts = Matrix::Matrix(0, nrow = 0, ncol = ncol(Y2),
      sparse = TRUE, dimnames = list(character(), colnames(Y2))), offsets = Matrix::Matrix(0, nrow = 0, ncol = ncol(Y2))))
    list(
      counts = if (length(counts_list) == 1) counts_list[[1]] else rbind(counts_list[[1]], counts_list[[2]]),
      offsets = if (length(offsets_list) == 1) offsets_list[[1]] else rbind(offsets_list[[1]], offsets_list[[2]])
    )
  }

  if (merge_triplets) {
    merged <- merge_triplet_layer(Y3, off3, Yw, offw)
    Y3 <- merged$counts
    off3 <- merged$offsets
  }

  if (verbose) message("Fitting edgeR (QL)…")
  res1 <- fit_layer(Y1, log_cells)
  res2 <- fit_layer(Y2, off2)
  res3 <- fit_layer(Y3, off3)
  resw <- if (!merge_triplets && !is.null(Yw)) fit_layer(Yw, offw) else NULL

  dag_result <- NULL
  if (fdr_method == "dagger") {
    if (verbose) message("Applying DAG FDR…")
    size1 <- if (!is.null(res1)) res1$motif else character()
    size2 <- if (!is.null(res2)) res2$motif else character()
    size3 <- if (!is.null(res3)) res3$motif else character()
    node_order <- c(size1, size2, size3)
    node_id <- stats::setNames(seq_along(node_order), node_order)

    edges <- list()
    if (length(size1) && length(size2)) {
      ab <- split_pair_labels(size2)
      for (i in seq_along(size2)) {
        a <- prefix_key(ab[i, 1], "N")
        b <- prefix_key(ab[i, 2], "N")
        if (a %in% size1) edges[[length(edges) + 1]] <- c(node_id[[a]], node_id[[size2[i]]])
        if (b %in% size1) edges[[length(edges) + 1]] <- c(node_id[[b]], node_id[[size2[i]]])
      }
    }
    if (length(size2) && length(size3)) {
      abc <- split_triplet_labels(size3)
      for (i in seq_along(size3)) {
        A <- abc[i, 1]
        B <- abc[i, 2]
        C <- abc[i, 3]
        for (pp in c(pair_key_vec(A, B, "E"), pair_key_vec(A, C, "E"), pair_key_vec(B, C, "E"))) {
          if (pp %in% size2) edges[[length(edges) + 1]] <- c(node_id[[pp]], node_id[[size3[i]]])
        }
      }
    }
    E_dag <- if (length(edges)) do.call(rbind, edges) else matrix(numeric(0), ncol = 2)

    q_dag <- rep(NA_real_, length(node_order))
    q1 <- if (length(size1)) stats::p.adjust(res1$PValue, "BH") else numeric(0)
    keep1 <- if (length(q1)) (q1 <= alpha) else logical(0)
    q2 <- rep(NA_real_, length(size2))
    if (length(size2)) {
      ab <- split_pair_labels(size2)
      ok <- logical(length(size2))
      for (i in seq_along(size2)) {
        ok[i] <- (prefix_key(ab[i, 1], "N") %in% size1[keep1]) && (prefix_key(ab[i, 2], "N") %in% size1[keep1])
      }
      pv2 <- res2$PValue
      pv2[!ok] <- NA
      if (any(ok)) q2[ok] <- stats::p.adjust(pv2[ok], "BH")
    }
    q3 <- rep(NA_real_, length(size3))
    if (length(size3)) {
      abc <- split_triplet_labels(size3)
      ok3 <- logical(length(size3))
      keep2_set <- if (length(size2)) size2[which(q2 <= alpha)] else character(0)
      for (i in seq_along(size3)) {
        prs <- c(
          pair_key_vec(abc[i, 1], abc[i, 2], "E"),
          pair_key_vec(abc[i, 1], abc[i, 3], "E"),
          pair_key_vec(abc[i, 2], abc[i, 3], "E")
        )
        ok3[i] <- all(prs %in% keep2_set)
      }
      pv3 <- res3$PValue
      pv3[!ok3] <- NA
      if (any(ok3)) q3[ok3] <- stats::p.adjust(pv3[ok3], "BH")
    }
    if (length(size1)) q_dag[node_id[size1]] <- q1
    if (length(size2)) q_dag[node_id[size2]] <- q2
    if (length(size3)) q_dag[node_id[size3]] <- q3

    add_q <- function(res, layer_keys) {
      if (is.null(res) || !length(layer_keys)) return(res)
      res$FDR_DAG <- q_dag[node_id[layer_keys]]
      res
    }
    res1 <- add_q(res1, size1)
    res2 <- add_q(res2, size2)
    res3 <- add_q(res3, size3)

    dag_result <- list(nodes = node_order, edges = E_dag)
  } else if (verbose) {
    message("Using layerwise BH FDR; DAG correction skipped.")
  }

  fdr_meta <- list(
    method = fdr_method,
    alpha = if (fdr_method == "dagger") alpha else NA_real_
  )

list(
    design = design,
    results = {
      ord <- function(df) if (is.null(df)) NULL else df[order(df$PValue), ]
      out <- list(size1 = ord(res1), size2 = ord(res2), size3 = ord(res3))
      if (!is.null(resw)) out$wedges <- ord(resw)
      out
    },
    dag = dag_result,
    fdr = fdr_meta
  )
}

#' @rdname cellmotif_edger
#' @keywords internal
cellmotif_edger <- function(...) {
  .Deprecated("motif_edger")
  motif_edger(...)
}
