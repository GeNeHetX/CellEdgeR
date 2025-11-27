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

build_delaunay_edges <- function(df) {
  tri <- deldir::deldir(df$x, df$y)
  segs <- tri$delsgs
  if (is.null(segs) || NROW(segs) == 0) {
    return(list(edges = matrix(numeric(0), ncol = 2), edge_len = numeric(0)))
  }
  e <- unique(cbind(segs$ind1, segs$ind2))
  colnames(e) <- c("from", "to")
  xy <- as.matrix(df[, c("x", "y")])
  len <- sqrt(rowSums((xy[e[, 1], , drop = FALSE] - xy[e[, 2], , drop = FALSE])^2))
  list(edges = e, edge_len = len)
}

format_label_triplets <- function(keys, lab_levels) {
  if (!length(keys)) return(character(0))
  ids_mat <- do.call(rbind, strsplit(keys, "\\|"))
  lab_mat <- matrix(lab_levels[as.integer(ids_mat)], ncol = ncol(ids_mat))
  apply(lab_mat, 1, function(v) paste(v, collapse = "|"))
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
#' @return A list of class \code{cellEdgeR_graphs} containing the sample names, global label set, and per-sample edges/labels.
#' @export
build_cell_graphs <- function(
  cells_by_sample,
  n_cores = 1,
  verbose = TRUE
) {
  if (!requireNamespace("deldir", quietly = TRUE)) {
    stop("Package deldir is required. Install it with install.packages('deldir').")
  }
  samples <- names(cells_by_sample)
  if (is.null(samples) || length(samples) == 0) stop("cells_by_sample must be a *named* list.")
  if (verbose) message("Samples: ", paste(samples, collapse = ", "))

  if (verbose) message("Standardizing sample frames (columns 1/2 as coords, 3 as labels)…")
  cells_standard <- stats::setNames(
    lapply(samples, function(s) standardize_sample_df(cells_by_sample[[s]], s)),
    samples
  )

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
      edges_info <- build_delaunay_edges(df)
      list(
        sample = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len
      )
    }, mc.cores = n_cores)
  } else {
    lapply(samples, function(s, ...) {
      df <- cells_standard[[s]]
      labs_chr <- as.character(df$label)
      labs_id <- unname(lab_to_id[labs_chr])
      edges_info <- build_delaunay_edges(df)
      list(
        sample = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len
      )
    })
  }
  names(per_sample) <- samples

  structure(
    list(
      samples = samples,
      label_levels = lab_levels,
      lab_to_id = lab_to_id,
      per_sample = per_sample
    ),
    class = "cellEdgeR_graphs"
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
#'
#' @return A list with elements \code{samples}, \code{label_levels}, \code{counts} (sparse matrices for size1/2/3 plus optional wedges),
#'   \code{exposure} (edge, triangle, cell, and optional wedge totals), and \code{meta} (parameters).
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
  verbose = TRUE
) {
  graphs <- NULL
  if (!is.null(graph_obj)) {
    if (!inherits(graph_obj, "cellEdgeR_graphs")) {
      stop("graph_obj must come from build_cell_graphs().")
    }
    graphs <- graph_obj
  } else if (!is.null(cells_by_sample) && inherits(cells_by_sample, "cellEdgeR_graphs")) {
    graphs <- cells_by_sample
    cells_by_sample <- NULL
  } else if (!is.null(cells_by_sample)) {
    graphs <- build_cell_graphs(cells_by_sample, n_cores = n_cores, verbose = verbose)
  } else {
    stop("Provide either cells_by_sample or graph_obj.")
  }

  count_motifs_from_graphs(
    graphs = graphs,
    max_edge_len = max_edge_len,
    include_wedges = include_wedges,
    verbose = verbose
  )
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
    names(tab) <- lab_levels
    tab
  })
  Y1 <- do.call(cbind, sing_list)
  dimnames(Y1) <- list(lab_levels, samples)
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
    key <- ifelse(la <= lb, paste(la, lb, sep = "|"), paste(lb, la, sep = "|"))
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
      keys_lab <- format_label_triplets(keys, lab_levels)
      tris_counts_list[[s]] <- structure(as.integer(out$counts), names = keys_lab)
    } else {
      tris_counts_list[[s]] <- integer(0)
    }
    tri_totals[s] <- as.numeric(out$tri_total)
    if (include_wedges) {
      wedge_keys <- as.character(out$wedge_keys)
      if (length(wedge_keys)) {
        wedge_lab <- format_label_triplets(wedge_keys, lab_levels)
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
build_pair_offsets <- function(Y2, Y1, log_edges, pseudo) {
  if (nrow(Y2) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y2)))
  }
  ab <- do.call(rbind, strsplit(rownames(Y2), "\\|"))
  a <- ab[, 1]
  b <- ab[, 2]
  NA_ <- Matrix::Matrix(0, nrow = length(a), ncol = ncol(Y1), dimnames = list(rownames(Y2), colnames(Y1)))
  NB_ <- Matrix::Matrix(0, nrow = length(b), ncol = ncol(Y1), dimnames = list(rownames(Y2), colnames(Y1)))
  idx_a <- match(a, rownames(Y1))
  idx_b <- match(b, rownames(Y1))
  sel_a <- which(!is.na(idx_a))
  sel_b <- which(!is.na(idx_b))
  if (length(sel_a)) NA_[sel_a, ] <- as.matrix(Y1[idx_a[sel_a], , drop = FALSE])
  if (length(sel_b)) NB_[sel_b, ] <- as.matrix(Y1[idx_b[sel_b], , drop = FALSE])
  off <- log(pmax(NA_, pseudo)) + log(pmax(NB_, pseudo))
  off <- sweep(off, 2, log_edges[colnames(Y2)], FUN = "+")
  dimnames(off) <- dimnames(Y2)
  off
}

#' @keywords internal
build_tri_offsets <- function(Y3, Y2, log_tris, pseudo) {
  if (nrow(Y3) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y3)))
  }
  abc <- do.call(rbind, strsplit(rownames(Y3), "\\|"))
  a <- abc[, 1]
  b <- abc[, 2]
  c <- abc[, 3]
  pair_key_vec <- function(u, v) {
    idx <- u <= v
    out <- character(length(u))
    out[idx] <- paste(u[idx], v[idx], sep = "|")
    out[!idx] <- paste(v[!idx], u[!idx], sep = "|")
    out
  }
  ABk <- pair_key_vec(a, b)
  ACk <- pair_key_vec(a, c)
  BCk <- pair_key_vec(b, c)
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
  off <- sweep(off, 2, log_tris[colnames(Y3)], FUN = "+")
  dimnames(off) <- dimnames(Y3)
  off
}

compute_motif_offsets <- function(motif_obj, pseudo) {
  samples <- motif_obj$samples
  Y1 <- motif_obj$counts$size1
  Y2 <- motif_obj$counts$size2
  Y3 <- motif_obj$counts$size3

  log_cells <- matrix(
    log(pmax(as.numeric(motif_obj$exposure$cells), 1)),
    nrow = max(1, nrow(Y1)),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(rownames(Y1), samples)
  )
  log_edges <- log(pmax(as.numeric(motif_obj$exposure$edges), 1))
  names(log_edges) <- samples
  log_tris <- log(pmax(as.numeric(motif_obj$exposure$triangles), 1))
  names(log_tris) <- samples

  offset2 <- build_pair_offsets(Y2, Y1, log_edges, pseudo)
  offset3 <- build_tri_offsets(Y3, Y2, log_tris, pseudo)
  list(
    size1 = Matrix::Matrix(log_cells, sparse = TRUE, dimnames = dimnames(Y1)),
    size2 = Matrix::Matrix(offset2, sparse = TRUE, dimnames = dimnames(Y2)),
    size3 = Matrix::Matrix(offset3, sparse = TRUE, dimnames = dimnames(Y3))
  )
}

#' Normalize motif counts via the statistical offsets
#'
#' Divide each motif matrix entry by the per-sample offset used in the downstream edgeR fits.
#'
#' @param motif_obj Object returned by [count_motifs_graphs()].
#' @param pseudo Small positive constant that was also provided to [motif_edger()]; defaults to `0.5`.
#' @return List with `size1`, `size2`, and `size3` count matrices matching `motif_obj$counts`.
#' @export
normalize_motif_counts <- function(motif_obj, pseudo = 0.5) {
  requireNamespace("Matrix", quietly = TRUE) || stop("Package Matrix is required.")
  offsets <- compute_motif_offsets(motif_obj, pseudo)

  normalize_layer <- function(count_mat, log_offset) {
    if (nrow(count_mat) == 0) {
      return(Matrix::Matrix(0, nrow = 0, ncol = ncol(count_mat), sparse = TRUE,
        dimnames = dimnames(count_mat)))
    }
    norm <- as.matrix(count_mat)
    norm <- norm / exp(as.matrix(log_offset))
    Matrix::Matrix(norm, sparse = TRUE, dimnames = dimnames(count_mat))
  }

  list(
    size1 = normalize_layer(motif_obj$counts$size1, offsets$size1),
    size2 = normalize_layer(motif_obj$counts$size2, offsets$size2),
    size3 = normalize_layer(motif_obj$counts$size3, offsets$size3)
  )
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
cellmotif_edger <- function(
  motif_obj,
  sample_df,
  design_formula,
  coef = NULL,
  pseudo = 0.5,
  alpha = 0.05,
  verbose = TRUE,
  fdr_method = c("bh", "dagger")
) {
  #"edgeR", "Matrix"
  
  samples <- motif_obj$samples
  if (is.null(rownames(sample_df))) stop("sample_df must have rownames = sample IDs")
  sample_df <- as.data.frame(sample_df[samples, , drop = FALSE])
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_df)
  if (is.null(coef)) {
    coef <- which(colnames(design) != "(Intercept)")[1]
    if (length(coef) == 0) stop("No non-intercept term found; set `coef`.")
  } else if (is.character(coef)) {
    coef <- match(coef, colnames(design))
    if (is.na(coef)) stop("coef name not found.")
  }

  fdr_method <- match.arg(fdr_method)

  Y1 <- motif_obj$counts$size1
  Y2 <- motif_obj$counts$size2
  Y3 <- motif_obj$counts$size3

  offsets <- compute_motif_offsets(motif_obj, pseudo)
  log_cells <- offsets$size1
  off2 <- offsets$size2
  off3 <- offsets$size3

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

  if (verbose) message("Fitting edgeR (QL)…")
  res1 <- fit_layer(Y1, log_cells)
  res2 <- fit_layer(Y2, off2)
  res3 <- fit_layer(Y3, off3)

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
      ab <- do.call(rbind, strsplit(size2, "\\|"))
      for (i in seq_along(size2)) {
        a <- ab[i, 1]
        b <- ab[i, 2]
        if (a %in% size1) edges[[length(edges) + 1]] <- c(node_id[[a]], node_id[[size2[i]]])
        if (b %in% size1) edges[[length(edges) + 1]] <- c(node_id[[b]], node_id[[size2[i]]])
      }
    }
    if (length(size2) && length(size3)) {
      abc <- do.call(rbind, strsplit(size3, "\\|"))
      pair_key <- function(u, v) if (u <= v) paste(u, v, sep = "|") else paste(v, u, sep = "|")
      for (i in seq_along(size3)) {
        A <- abc[i, 1]
        B <- abc[i, 2]
        C <- abc[i, 3]
        for (pp in c(pair_key(A, B), pair_key(A, C), pair_key(B, C))) {
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
      ab <- do.call(rbind, strsplit(size2, "\\|"))
      ok <- logical(length(size2))
      for (i in seq_along(size2)) {
        ok[i] <- (ab[i, 1] %in% size1[keep1]) && (ab[i, 2] %in% size1[keep1])
      }
      pv2 <- res2$PValue
      pv2[!ok] <- NA
      if (any(ok)) q2[ok] <- stats::p.adjust(pv2[ok], "BH")
    }
    q3 <- rep(NA_real_, length(size3))
    if (length(size3)) {
      abc <- do.call(rbind, strsplit(size3, "\\|"))
      pair_key <- function(u, v) if (u <= v) paste(u, v, sep = "|") else paste(v, u, sep = "|")
      ok3 <- logical(length(size3))
      keep2_set <- if (length(size2)) size2[which(q2 <= alpha)] else character(0)
      for (i in seq_along(size3)) {
        prs <- c(
          pair_key(abc[i, 1], abc[i, 2]),
          pair_key(abc[i, 1], abc[i, 3]),
          pair_key(abc[i, 2], abc[i, 3])
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
    results = list(size1 = res1[order(res1$PValue),], size2 = res2[order(res2$PValue),], size3 = res3[order(res3$PValue),]),
    dag = dag_result,
    fdr = fdr_meta
  )
}

#' @rdname cellmotif_edger
#' @export
motif_edger <- cellmotif_edger
