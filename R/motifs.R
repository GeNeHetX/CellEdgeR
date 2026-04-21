# helper utilities -------------------------------------------------------------

CELL_EDGER_OFFSET_VERSION <- "volume_chung_lu_v2"

standardize_sample_df <- function(sample_df, sample_name) {
  df <- as.data.frame(sample_df)
  if (ncol(df) < 3L) {
    stop("Sample '", sample_name, "' must have at least three columns (x, y, labels).")
  }

  coord_cols <- if (all(c("x", "y") %in% names(df))) c("x", "y") else names(df)[seq_len(2)]
  coords <- df[, coord_cols, drop = FALSE]
  numeric_ok <- vapply(coords, is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop("Sample '", sample_name, "' requires numeric coordinates in columns 1 and 2 (or named x/y).")
  }

  label_column <- if ("label" %in% names(df)) df[["label"]] else df[[3]]
  if (!is.atomic(label_column) || length(label_column) != nrow(df)) {
    stop("Sample '", sample_name, "' label column must be a vector with one entry per row.")
  }

  out <- data.frame(
    x = coords[[1]],
    y = coords[[2]],
    label = label_column,
    stringsAsFactors = FALSE
  )

  boundary_cols <- c("boundary", "is_boundary", "on_boundary", "is_border", "on_border", "edge_cell", "boundary_flag")
  boundary_name <- boundary_cols[boundary_cols %in% names(df)][1]
  if (!is.na(boundary_name) && length(boundary_name)) {
    boundary_vec <- df[[boundary_name]]
    if (!is.atomic(boundary_vec) || length(boundary_vec) != nrow(df)) {
      stop("Sample '", sample_name, "' boundary column must be a vector with one entry per row.")
    }
    boundary_vec <- as.logical(boundary_vec)
    if (anyNA(boundary_vec)) {
      stop("Sample '", sample_name, "' boundary column contains NA values.")
    }
    out$boundary <- boundary_vec
  }
  out
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
  if (is.null(graph_obj$sample_name) || !length(graph_obj$sample_name)) stop("graph_obj has no sample names.")
  if (anyDuplicated(graph_obj$sample_name)) stop("graph_obj sample names must be unique.")
  if (is.null(graph_obj$per_sample_graph) || !length(graph_obj$per_sample_graph)) {
    stop("graph_obj is missing per_sample_graph; count_motifs_graphs() and plot_sample_graph() require per-sample graphs. ",
      "Rebuild with build_cell_graphs() or keep per_sample_graph when trimming to save memory.")
  }
  if (!is.list(graph_obj$per_sample_graph)) stop("graph_obj$per_sample_graph must be a list.")
  missing <- setdiff(graph_obj$sample_name, names(graph_obj$per_sample_graph))
  if (length(missing)) stop("graph_obj missing per-sample entries: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

available_offset_modes <- function(obj) {
  offs <- obj$offsets
  if (is.null(offs) || !is.list(offs) || !length(offs)) return(character())
  if (!is.null(names(offs)) && all(vapply(offs, function(x) is.list(x) && any(c("node", "edge", "triangle") %in% names(x)), logical(1)))) {
    return(names(offs))
  }
  if (all(c("node", "edge", "triangle") %in% names(offs))) {
    return("volume")
  }
  character()
}

validate_motif_obj <- function(motif_obj, require_offsets = TRUE) {
  if (!is.list(motif_obj) || is.null(motif_obj$raw_count)) {
    stop("Object must be the list returned by count_motifs_graphs().")
  }
  if (!inherits(motif_obj, "cellEdgeR_obj")) {
    stop("Object must come from build_cell_graphs()/count_motifs_graphs().")
  }
  if (is.null(motif_obj$sample_name) || !length(motif_obj$sample_name)) {
    stop("Object is missing sample names.")
  }
  if (!is.list(motif_obj$raw_count)) stop("raw_count must be a list of sparse matrices.")
  if (is.null(motif_obj$exposure)) stop("Object is missing exposure totals.")
  if (require_offsets && length(available_offset_modes(motif_obj)) == 0) {
    stop("Object is missing offsets; rerun count_motifs_graphs() to compute them.")
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

validate_edge_indices <- function(edges, n, sample_name) {
  if (is.null(edges) || length(edges) == 0) return(edges)
  if (!is.matrix(edges) || ncol(edges) != 2L) {
    stop("Sample '", sample_name, "' edges must be a 2-column matrix.")
  }
  storage.mode(edges) <- "integer"
  if (anyNA(edges)) {
    stop("Sample '", sample_name, "' edges contain NA indices.")
  }
  if (any(edges < 1L | edges > n)) {
    stop("Sample '", sample_name, "' edges have indices outside 1..n.")
  }
  edges
}

coerce_boundary_mask <- function(boundary, n, sample_name, source = "boundary") {
  if (is.null(boundary)) return(NULL)
  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("Sample '", sample_name, "' has invalid node count for ", source, ".")
  }
  if (is.logical(boundary)) {
    if (length(boundary) != n) {
      stop("Sample '", sample_name, "' ", source, " mask must be length ", n, ".")
    }
    if (anyNA(boundary)) stop("Sample '", sample_name, "' ", source, " mask contains NA values.")
    return(boundary)
  }
  if (is.numeric(boundary)) {
    if (length(boundary) == n && all(boundary %in% c(0, 1))) {
      boundary <- as.logical(boundary)
      if (anyNA(boundary)) stop("Sample '", sample_name, "' ", source, " mask contains NA values.")
      return(boundary)
    }
    if (any(boundary %% 1 != 0)) {
      stop("Sample '", sample_name, "' ", source, " indices must be integers.")
    }
    idx <- as.integer(boundary)
    if (anyNA(idx) || any(idx < 1L | idx > n)) {
      stop("Sample '", sample_name, "' ", source, " indices are out of range.")
    }
    mask <- rep(FALSE, n)
    if (length(idx)) mask[idx] <- TRUE
    return(mask)
  }
  stop("Sample '", sample_name, "' ", source, " must be a logical mask or integer indices.")
}

normalize_erosion_cells <- function(erosion_cells, samples) {
  if (is.null(erosion_cells)) return(NULL)
  if (!is.list(erosion_cells)) stop("erosion_cells must be NULL or a list.")
  if (is.null(names(erosion_cells)) || all(names(erosion_cells) == "")) {
    if (length(erosion_cells) == length(samples)) {
      names(erosion_cells) <- samples
    } else if (length(erosion_cells) == 1L) {
      erosion_cells <- rep(erosion_cells, length(samples))
      names(erosion_cells) <- samples
    } else {
      stop("erosion_cells must be named or length 1 or length = number of samples.")
    }
  }
  erosion_cells
}

#' Build spatial graphs via Delaunay triangulation
#'
#' Compute the underlying Delaunay edges for each sample without pruning so that motif
#' counts can be recomputed with different thresholds or wedge options without rebuilding.
#'
#' @param cells_by_sample Named list of data frames where the first two columns are numeric coordinates
#'   and the third holds cell-type labels (or named `x`, `y`, `label`). An optional logical boundary column
#'   (e.g., `boundary`/`is_boundary`) can be supplied to enable erosion during motif counting.
#' @param n_cores Parallelism hint; values greater than 1 trigger `parallel::mclapply` on Unix-alikes, while 1 runs sequentially.
#' @param verbose Logical; print progress.
#' @return A cellgraph object (class `cellEdgeR_obj`) containing:
#' \describe{
#' \item{sample_name}{Sample names in the order used throughout.}
#' \item{label_levels/lab_to_id}{Canonical label set and integer map.}
#' \item{per_sample_graph}{Per-sample list with node count, labels (chr/int), edges, edge lengths, and coordinates.}
#' \item{raw_count/exposure/offsets/norm_counts/relative_counts}{Empty placeholders until [count_motifs_graphs()] is called.}
#' \item{edger}{Empty placeholder until [motif_edger()] is called.}
#' \item{parameters}{Build parameters (placeholders for max_edge_len, offset_pseudo, built_from).}
#' }
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
      boundary <- NULL
      if ("boundary" %in% names(df)) {
        boundary <- as.logical(df$boundary)
        if (length(boundary) != nrow(df)) {
          stop("Sample '", s, "' boundary mask must match number of rows.")
        }
        if (anyNA(boundary)) stop("Sample '", s, "' boundary mask contains NA values.")
      }
      list(
        sample_name = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len,
        xy = as.matrix(df[, c("x", "y")]),
        boundary = boundary
      )
    }, mc.cores = n_cores)
  } else {
    lapply(samples, function(s, ...) {
      df <- cells_standard[[s]]
      labs_chr <- as.character(df$label)
      labs_id <- unname(lab_to_id[labs_chr])
      edges_info <- build_delaunay_edges(df, verbose = verbose)
      boundary <- NULL
      if ("boundary" %in% names(df)) {
        boundary <- as.logical(df$boundary)
        if (length(boundary) != nrow(df)) {
          stop("Sample '", s, "' boundary mask must match number of rows.")
        }
        if (anyNA(boundary)) stop("Sample '", s, "' boundary mask contains NA values.")
      }
      list(
        sample_name = s,
        n = nrow(df),
        labels_id = labs_id,
        labels_chr = labs_chr,
        edges = edges_info$edges,
        edge_len = edges_info$edge_len,
        xy = as.matrix(df[, c("x", "y")]),
        boundary = boundary
      )
    })
  }
  names(per_sample) <- samples

  structure(
    list(
      sample_name = samples,
      label_levels = lab_levels,
      lab_to_id = lab_to_id,
      per_sample_graph = per_sample,
      raw_count = NULL,
      exposure = NULL,
      offsets = NULL,
      norm_counts = NULL,
      relative_counts = NULL,
      offset_part_id = NULL,
      offset_part_values = NULL,
      edger = NULL,
      parameters = list(max_edge_len = NA_real_, include_wedge = FALSE, offset_pseudo = NA_real_, built_from = "cells")
    ),
    class = c("cellEdgeR_obj", "cellEdgeR_graphs")
  )
}

#' Count cell-type motifs on Delaunay graphs
#'
#' Reuse a prebuilt Delaunay triangulation, prune edges that exceed \code{max_edge_len},
#' and count cell-type singletons, unordered pairs (edges), unordered triplets (triangles), and optionally wedge.
#' Triangle counts are accelerated via the C++ helper exposed in \pkg{CellEdgeR}.
#'
#' @param graph_obj Output of \code{build_cell_graphs()}; the triangulation is reused and only pruning/wedge parameters are applied here.
#' @param max_edge_len Numeric threshold; Delaunay edges longer than this are dropped. Set to \code{NA}, \code{NULL}, or \code{<= 0} to skip pruning.
#' @param include_wedge Logical; if \code{TRUE}, returns open-triplet (wedge) counts alongside triangles.
#' @param verbose Logical; print progress.
#' @param offset_pseudo Small positive constant used inside offsets.
#' @param n_cores Parallelism hint for motif counting; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.
#' @param erosion Logical; if \code{TRUE} (default), exclude motifs that touch boundary cells. Boundary
#'   cells can be provided as a logical column (e.g., `boundary`/`is_boundary`) in the input data or via
#'   \code{erosion_cells}.
#' @param erosion_cells Optional list of boundary masks/indices by sample name. Each entry can be a logical
#'   vector (length = number of cells) or integer indices to exclude.
#'
#' @return A cellgraph object (class `cellEdgeR_obj`) with the original graph info plus:
#' \describe{
#' \item{raw_count}{Sparse matrices for `node`, `edge`, `triangle`, and optional `wedge` motifs (rows=motifs, cols=samples).}
#' \item{exposure}{Totals used in offsets: `cells`, `edges`, `triangles`, `volumes` per label×sample,
#'   `center_pairs` (sum of degree*(degree-1) per label×sample), plus `wedge` and `triples` (centered open+closed triples)
#'   when requested.}
#' \item{offsets}{Volume offsets containing log-expected matrices per layer.
#'   Node offsets are retained for normalization/submotif inspection; tested layers use structural volume offsets.}
#' \item{norm_counts}{Normalized counts for the volume offset (counts / exp(offset)).}
#' \item{relative_counts}{edgeR intercept-only log2 residuals for the volume offset.}
#' \item{offset_part_id}{Components contributing to each motif's volume offset.}
#' \item{offset_part_values}{Numeric values referenced by `offset_part_id` (e.g., volumes, 2m, edge posteriors).}
#' \item{edger}{edgeR fits/tests once [motif_edger()] is run.}
#' \item{parameters}{Run parameters (max_edge_len, include_wedge, offset_pseudo, available offset_modes, layer names,
#'   node_tmm_offsets flag).}
#' }
#' @examples
#' demo <- list(
#'   s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
#'   s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
#' )
#' graphs <- build_cell_graphs(demo, verbose = FALSE)
#' counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
#' counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
#' str(counts_full$raw_count)
#' @rdname count_motifs_graphs
#' @export
count_motifs_graphs <- function(
  graph_obj,
  max_edge_len = NA,
  include_wedge = FALSE,
  verbose = TRUE,
  offset_pseudo = 1,
  n_cores = 1,
  erosion = TRUE,
  erosion_cells = NULL
) {
  validate_graph_obj(graph_obj)
  graphs <- graph_obj

  counts_obj <- count_motifs_from_graphs(
    graphs = graphs,
    max_edge_len = max_edge_len,
    include_wedge = include_wedge,
    verbose = verbose,
    n_cores = n_cores,
    erosion = erosion,
    erosion_cells = erosion_cells
  )
  offsets_volume <- compute_offsets_volume(counts_obj, offset_pseudo)
  offset_results <- list(volume = offsets_volume)
  offset_results <- lapply(offset_results, function(res) {
    list(
      offsets = res$offsets,
      parts = res$parts,
      part_values = res$part_values,
      norm_counts = normalize_counts_simple(counts_obj$raw_count, res$offsets),
      relative_counts = compute_relative_counts(counts_obj$raw_count, res$offsets, samples = graphs$sample_name, verbose = verbose)
    )
  })

  graphs$raw_count <- counts_obj$raw_count
  graphs$exposure <- counts_obj$exposure
  graphs$offsets <- lapply(offset_results, `[[`, "offsets")
  graphs$offset_part_id <- lapply(offset_results, `[[`, "parts")
  graphs$offset_part_values <- lapply(offset_results, `[[`, "part_values")
  graphs$norm_counts <- lapply(offset_results, `[[`, "norm_counts")
  graphs$relative_counts <- lapply(offset_results, `[[`, "relative_counts")
  graphs$edger <- NULL
  graphs$parameters$max_edge_len <- max_edge_len
  graphs$parameters$include_wedge <- include_wedge
  graphs$parameters$offset_pseudo <- offset_pseudo
  graphs$parameters$offset_modes <- names(offset_results)
  graphs$parameters$offset_version <- CELL_EDGER_OFFSET_VERSION
  graphs$parameters$node_tmm_offsets <- isTRUE(offsets_volume$node_tmm_applied)
  graphs$parameters$layer_names <- list(node = "node", edge = "edge", triangle = "triangle", wedge = "wedge")
  graphs$parameters$erosion <- isTRUE(erosion)
  graphs$parameters$erosion_cells_provided <- !is.null(erosion_cells)
  graphs
}

#' Merge two motif objects
#'
#' Combine two `cellEdgeR_obj` outputs from [count_motifs_graphs()] into a single
#' object with combined samples. Offsets and normalized counts are recomputed
#' for the merged data.
#'
#' @param motif_obj_a,motif_obj_b `cellEdgeR_obj` objects returned by [count_motifs_graphs()].
#' @param verbose Logical; print progress while recomputing offsets.
#' @return A merged `cellEdgeR_obj` containing all samples; offsets, normalized
#'   counts, and relative counts are recomputed for the combined data. Per-sample
#'   graphs are kept only when present in both inputs.
#' @export
merge_motif_objs <- function(motif_obj_a, motif_obj_b, verbose = TRUE) {
  validate_motif_obj(motif_obj_a, require_offsets = FALSE)
  validate_motif_obj(motif_obj_b, require_offsets = FALSE)

  check_samples <- function(samples, obj_name) {
    if (is.null(samples) || !length(samples)) stop(obj_name, " has no sample names.")
    if (anyNA(samples) || any(samples == "")) stop(obj_name, " sample names must be non-empty, non-NA.")
    if (anyDuplicated(samples)) stop(obj_name, " sample names must be unique.")
    invisible(TRUE)
  }
  align_cols <- function(mat, samples, layer, obj_name) {
    if (is.null(mat)) stop(obj_name, " raw_count$", layer, " is missing.")
    if (is.null(colnames(mat))) stop(obj_name, " raw_count$", layer, " is missing column names.")
    if (anyDuplicated(colnames(mat))) {
      stop(obj_name, " raw_count$", layer, " has duplicated column names.")
    }
    if (!setequal(colnames(mat), samples)) {
      stop(obj_name, " raw_count$", layer, " columns do not match sample names.")
    }
    if (!identical(colnames(mat), samples)) mat <- mat[, samples, drop = FALSE]
    mat
  }
  align_exposure_matrix <- function(mat, samples, label, obj_name) {
    if (is.null(mat)) stop(obj_name, " exposure$", label, " is missing.")
    if (is.null(colnames(mat))) stop(obj_name, " exposure$", label, " is missing column names.")
    if (anyDuplicated(colnames(mat))) {
      stop(obj_name, " exposure$", label, " has duplicated column names.")
    }
    if (!setequal(colnames(mat), samples)) {
      stop(obj_name, " exposure$", label, " columns do not match sample names.")
    }
    if (!identical(colnames(mat), samples)) mat <- mat[, samples, drop = FALSE]
    if (is.null(rownames(mat))) stop(obj_name, " exposure$", label, " is missing row names.")
    mat
  }
  align_node_rows <- function(mat, labels, obj_name) {
    expected <- prefix_key(labels, "N")
    if (is.null(rownames(mat))) stop(obj_name, " raw_count$node is missing row names.")
    if (anyNA(rownames(mat)) || any(rownames(mat) == "")) {
      stop(obj_name, " raw_count$node has empty row names.")
    }
    if (anyDuplicated(rownames(mat))) stop(obj_name, " raw_count$node has duplicated row names.")
    if (!setequal(rownames(mat), expected)) {
      stop(obj_name, " raw_count$node row names must match label_levels.")
    }
    if (!identical(rownames(mat), expected)) mat <- mat[expected, , drop = FALSE]
    mat
  }
  check_layer_keys <- function(mat, labels, layer, obj_name) {
    keys <- rownames(mat)
    if (is.null(keys)) stop(obj_name, " raw_count$", layer, " is missing row names.")
    if (!length(keys)) return(invisible(TRUE))
    if (anyNA(keys) || any(keys == "")) stop(obj_name, " raw_count$", layer, " has empty motif keys.")
    if (anyDuplicated(keys)) stop(obj_name, " raw_count$", layer, " has duplicated motif keys.")
    prefix <- switch(layer, edge = "E_", triangle = "T_", wedge = "W_")
    if (!all(startsWith(keys, prefix))) {
      stop(obj_name, " raw_count$", layer, " uses unexpected motif prefixes.")
    }
    parts <- strsplit(strip_prefix(keys), "_", fixed = TRUE)
    expect_len <- if (layer == "edge") 2L else 3L
    lens <- vapply(parts, length, integer(1))
    if (any(lens != expect_len)) {
      stop(obj_name, " raw_count$", layer, " has malformed motif keys.")
    }
    used_labels <- unique(unlist(parts))
    missing <- setdiff(used_labels, labels)
    if (length(missing)) {
      stop(obj_name, " raw_count$", layer, " references labels not in label_levels: ",
        paste(missing, collapse = ", "))
    }
    invisible(TRUE)
  }
  align_vector <- function(vec, samples, label, obj_name) {
    if (is.null(vec)) stop(obj_name, " exposure$", label, " is missing.")
    if (is.null(names(vec))) stop(obj_name, " exposure$", label, " must be named by sample names.")
    if (anyDuplicated(names(vec))) stop(obj_name, " exposure$", label, " has duplicated sample names.")
    if (!setequal(names(vec), samples)) {
      stop(obj_name, " exposure$", label, " names do not match sample names.")
    }
    vec[samples]
  }
  align_volume <- function(mat, labels, samples, obj_name) {
    if (is.null(mat)) stop(obj_name, " exposure$volumes is missing.")
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop(obj_name, " exposure$volumes must have row and column names.")
    }
    if (anyDuplicated(rownames(mat))) stop(obj_name, " exposure$volumes has duplicated row names.")
    if (anyDuplicated(colnames(mat))) stop(obj_name, " exposure$volumes has duplicated column names.")
    if (!setequal(rownames(mat), labels)) {
      stop(obj_name, " exposure$volumes rows do not match label_levels.")
    }
    if (!setequal(colnames(mat), samples)) {
      stop(obj_name, " exposure$volumes columns do not match sample names.")
    }
    mat[labels, samples, drop = FALSE]
  }
  align_label_matrix <- function(mat, labels, samples, label, obj_name) {
    if (is.null(mat)) stop(obj_name, " exposure$", label, " is missing.")
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop(obj_name, " exposure$", label, " must have row and column names.")
    }
    if (anyDuplicated(rownames(mat))) stop(obj_name, " exposure$", label, " has duplicated row names.")
    if (anyDuplicated(colnames(mat))) stop(obj_name, " exposure$", label, " has duplicated column names.")
    if (!setequal(rownames(mat), labels)) {
      stop(obj_name, " exposure$", label, " rows do not match label_levels.")
    }
    if (!setequal(colnames(mat), samples)) {
      stop(obj_name, " exposure$", label, " columns do not match sample names.")
    }
    mat[labels, samples, drop = FALSE]
  }
  expand_sparse_rows <- function(mat, target_rows, layer, obj_name) {
    if (is.null(target_rows) || !length(target_rows)) {
      return(Matrix::Matrix(0, nrow = 0L, ncol = ncol(mat), sparse = TRUE,
        dimnames = list(character(), colnames(mat))))
    }
    row_idx <- match(rownames(mat), target_rows)
    if (anyNA(row_idx)) {
      stop(obj_name, " raw_count$", layer, " has motifs not present in the merged set.")
    }
    if (nrow(mat) == 0L) {
      return(Matrix::Matrix(0, nrow = length(target_rows), ncol = ncol(mat), sparse = TRUE,
        dimnames = list(target_rows, colnames(mat))))
    }
    sm <- Matrix::summary(mat)
    if (!nrow(sm)) {
      return(Matrix::Matrix(0, nrow = length(target_rows), ncol = ncol(mat), sparse = TRUE,
        dimnames = list(target_rows, colnames(mat))))
    }
    Matrix::sparseMatrix(
      i = row_idx[sm$i],
      j = sm$j,
      x = sm$x,
      dims = c(length(target_rows), ncol(mat)),
      dimnames = list(target_rows, colnames(mat))
    )
  }
  expand_dense_rows <- function(mat, target_rows, obj_name, label = "volumes") {
    if (is.null(rownames(mat))) stop(obj_name, " exposure$", label, " is missing row names.")
    out <- matrix(0, nrow = length(target_rows), ncol = ncol(mat),
      dimnames = list(target_rows, colnames(mat)))
    if (!nrow(mat)) return(out)
    idx <- match(rownames(mat), target_rows)
    if (anyNA(idx)) stop(obj_name, " exposure$volumes has unknown labels.")
    out[idx, ] <- mat
    out
  }
  remap_graphs <- function(per_sample_graph, old_labels, new_map, obj_name) {
    lapply(per_sample_graph, function(ps) {
      labels_chr <- ps$labels_chr
      if (is.null(labels_chr) && !is.null(ps$labels_id)) {
        labels_chr <- old_labels[ps$labels_id]
        ps$labels_chr <- labels_chr
      }
      if (is.null(labels_chr)) {
        stop(obj_name, " per_sample_graph entry is missing labels.")
      }
      new_ids <- unname(new_map[as.character(labels_chr)])
      if (anyNA(new_ids)) stop(obj_name, " per_sample_graph has labels not in merged label_levels.")
      ps$labels_id <- as.integer(new_ids)
      ps
    })
  }
  pick_param <- function(name, default = NULL) {
    val_a <- motif_obj_a$parameters[[name]]
    val_b <- motif_obj_b$parameters[[name]]
    if (!is.null(val_a) && !is.null(val_b)) {
      if (!isTRUE(all.equal(val_a, val_b))) {
        stop("motif_obj_a and motif_obj_b have different ", name, " values.")
      }
      return(val_a)
    }
    if (!is.null(val_a)) return(val_a)
    if (!is.null(val_b)) return(val_b)
    default
  }

  samples_a <- as.character(motif_obj_a$sample_name)
  samples_b <- as.character(motif_obj_b$sample_name)
  check_samples(samples_a, "motif_obj_a")
  check_samples(samples_b, "motif_obj_b")
  dup_samples <- intersect(samples_a, samples_b)
  if (length(dup_samples)) {
    stop("Sample names overlap between motif_obj_a and motif_obj_b: ", paste(dup_samples, collapse = ", "))
  }

  labels_a <- as.character(motif_obj_a$label_levels)
  labels_b <- as.character(motif_obj_b$label_levels)
  if (!length(labels_a) || !length(labels_b)) stop("Both motif objects must have label_levels.")
  if (anyNA(labels_a) || any(labels_a == "")) stop("motif_obj_a label_levels must be non-empty, non-NA.")
  if (anyNA(labels_b) || any(labels_b == "")) stop("motif_obj_b label_levels must be non-empty, non-NA.")
  if (anyDuplicated(labels_a)) stop("motif_obj_a label_levels must be unique.")
  if (anyDuplicated(labels_b)) stop("motif_obj_b label_levels must be unique.")
  overlap <- intersect(labels_a, labels_b)
  if (!length(overlap)) {
    stop("motif_obj_a and motif_obj_b have no overlapping labels; cannot merge.")
  }

  required_layers <- c("node", "edge", "triangle")
  missing_a <- setdiff(required_layers, names(motif_obj_a$raw_count))
  missing_b <- setdiff(required_layers, names(motif_obj_b$raw_count))
  if (length(missing_a)) stop("motif_obj_a is missing raw_count layers: ", paste(missing_a, collapse = ", "))
  if (length(missing_b)) stop("motif_obj_b is missing raw_count layers: ", paste(missing_b, collapse = ", "))

  has_wedge_a <- "wedge" %in% names(motif_obj_a$raw_count) && !is.null(motif_obj_a$raw_count$wedge)
  has_wedge_b <- "wedge" %in% names(motif_obj_b$raw_count) && !is.null(motif_obj_b$raw_count$wedge)
  if (xor(has_wedge_a, has_wedge_b)) {
    stop("Cannot merge motif objects when only one includes wedge counts.")
  }
  include_wedge <- has_wedge_a

  param_include_a <- motif_obj_a$parameters$include_wedge
  param_include_b <- motif_obj_b$parameters$include_wedge
  if (!is.null(param_include_a) && !is.na(param_include_a) && isTRUE(param_include_a) != include_wedge) {
    stop("motif_obj_a include_wedge parameter does not match raw_count layers.")
  }
  if (!is.null(param_include_b) && !is.na(param_include_b) && isTRUE(param_include_b) != include_wedge) {
    stop("motif_obj_b include_wedge parameter does not match raw_count layers.")
  }

  max_edge_len <- pick_param("max_edge_len", default = NA_real_)
  offset_pseudo <- pick_param("offset_pseudo", default = 1)
  built_from <- pick_param("built_from", default = "cells")
  erosion <- pick_param("erosion", default = TRUE)
  erosion_cells_provided <- pick_param("erosion_cells_provided", default = FALSE)

  raw_a <- motif_obj_a$raw_count
  raw_b <- motif_obj_b$raw_count
  raw_a$node <- align_node_rows(align_cols(raw_a$node, samples_a, "node", "motif_obj_a"), labels_a, "motif_obj_a")
  raw_b$node <- align_node_rows(align_cols(raw_b$node, samples_b, "node", "motif_obj_b"), labels_b, "motif_obj_b")
  raw_a$edge <- align_cols(raw_a$edge, samples_a, "edge", "motif_obj_a")
  raw_b$edge <- align_cols(raw_b$edge, samples_b, "edge", "motif_obj_b")
  raw_a$triangle <- align_cols(raw_a$triangle, samples_a, "triangle", "motif_obj_a")
  raw_b$triangle <- align_cols(raw_b$triangle, samples_b, "triangle", "motif_obj_b")
  check_layer_keys(raw_a$edge, labels_a, "edge", "motif_obj_a")
  check_layer_keys(raw_b$edge, labels_b, "edge", "motif_obj_b")
  check_layer_keys(raw_a$triangle, labels_a, "triangle", "motif_obj_a")
  check_layer_keys(raw_b$triangle, labels_b, "triangle", "motif_obj_b")
  if (include_wedge) {
    raw_a$wedge <- align_cols(raw_a$wedge, samples_a, "wedge", "motif_obj_a")
    raw_b$wedge <- align_cols(raw_b$wedge, samples_b, "wedge", "motif_obj_b")
    check_layer_keys(raw_a$wedge, labels_a, "wedge", "motif_obj_a")
    check_layer_keys(raw_b$wedge, labels_b, "wedge", "motif_obj_b")
  }

  exp_a <- motif_obj_a$exposure
  exp_b <- motif_obj_b$exposure
  exp_a$edges <- align_vector(exp_a$edges, samples_a, "edges", "motif_obj_a")
  exp_b$edges <- align_vector(exp_b$edges, samples_b, "edges", "motif_obj_b")
  exp_a$triangles <- align_vector(exp_a$triangles, samples_a, "triangles", "motif_obj_a")
  exp_b$triangles <- align_vector(exp_b$triangles, samples_b, "triangles", "motif_obj_b")
  exp_a$cells <- align_vector(exp_a$cells, samples_a, "cells", "motif_obj_a")
  exp_b$cells <- align_vector(exp_b$cells, samples_b, "cells", "motif_obj_b")
  exp_a$volumes <- align_volume(exp_a$volumes, labels_a, samples_a, "motif_obj_a")
  exp_b$volumes <- align_volume(exp_b$volumes, labels_b, samples_b, "motif_obj_b")
  exp_a$center_pairs <- if (!is.null(exp_a$center_pairs)) {
    align_label_matrix(exp_a$center_pairs, labels_a, samples_a, "center_pairs", "motif_obj_a")
  } else {
    exp_a$volumes
  }
  exp_b$center_pairs <- if (!is.null(exp_b$center_pairs)) {
    align_label_matrix(exp_b$center_pairs, labels_b, samples_b, "center_pairs", "motif_obj_b")
  } else {
    exp_b$volumes
  }
  if (include_wedge) {
    exp_a$wedge <- align_vector(exp_a$wedge, samples_a, "wedge", "motif_obj_a")
    exp_b$wedge <- align_vector(exp_b$wedge, samples_b, "wedge", "motif_obj_b")
  }
  has_triples_a <- !is.null(exp_a$triples)
  has_triples_b <- !is.null(exp_b$triples)
  if (xor(has_triples_a, has_triples_b)) {
    stop("Cannot merge motif objects when only one includes triple exposures.")
  }
  if (has_triples_a) {
    exp_a$triples <- align_exposure_matrix(exp_a$triples, samples_a, "triples", "motif_obj_a")
    exp_b$triples <- align_exposure_matrix(exp_b$triples, samples_b, "triples", "motif_obj_b")
  }

  label_levels <- sort(unique(c(labels_a, labels_b)))
  lab_to_id <- stats::setNames(seq_len(length(label_levels)), label_levels)
  sample_name <- c(samples_a, samples_b)

  node_rows <- prefix_key(label_levels, "N")
  edge_rows <- unique(c(rownames(raw_a$edge), rownames(raw_b$edge)))
  tri_rows <- unique(c(rownames(raw_a$triangle), rownames(raw_b$triangle)))
  wedge_rows <- if (include_wedge) unique(c(rownames(raw_a$wedge), rownames(raw_b$wedge))) else character(0)

  Y1 <- cbind(
    expand_sparse_rows(raw_a$node, node_rows, "node", "motif_obj_a"),
    expand_sparse_rows(raw_b$node, node_rows, "node", "motif_obj_b")
  )
  Y2 <- cbind(
    expand_sparse_rows(raw_a$edge, edge_rows, "edge", "motif_obj_a"),
    expand_sparse_rows(raw_b$edge, edge_rows, "edge", "motif_obj_b")
  )
  Y3 <- cbind(
    expand_sparse_rows(raw_a$triangle, tri_rows, "triangle", "motif_obj_a"),
    expand_sparse_rows(raw_b$triangle, tri_rows, "triangle", "motif_obj_b")
  )
  raw_count <- list(node = Y1, edge = Y2, triangle = Y3)
  if (include_wedge) {
    Yw <- cbind(
      expand_sparse_rows(raw_a$wedge, wedge_rows, "wedge", "motif_obj_a"),
      expand_sparse_rows(raw_b$wedge, wedge_rows, "wedge", "motif_obj_b")
    )
    raw_count$wedge <- Yw
  }

  volumes <- cbind(
    expand_dense_rows(exp_a$volumes, label_levels, "motif_obj_a"),
    expand_dense_rows(exp_b$volumes, label_levels, "motif_obj_b")
  )
  center_pairs <- cbind(
    expand_dense_rows(exp_a$center_pairs, label_levels, "motif_obj_a", label = "center_pairs"),
    expand_dense_rows(exp_b$center_pairs, label_levels, "motif_obj_b", label = "center_pairs")
  )
  exposure <- list(
    edges = c(exp_a$edges, exp_b$edges),
    triangles = c(exp_a$triangles, exp_b$triangles),
    cells = c(exp_a$cells, exp_b$cells),
    volumes = volumes,
    center_pairs = center_pairs
  )
  if (include_wedge) exposure$wedge <- c(exp_a$wedge, exp_b$wedge)
  if (has_triples_a) {
    triple_rows <- unique(c(rownames(exp_a$triples), rownames(exp_b$triples)))
    triples <- cbind(
      expand_sparse_rows(exp_a$triples, triple_rows, "triples", "motif_obj_a"),
      expand_sparse_rows(exp_b$triples, triple_rows, "triples", "motif_obj_b")
    )
    exposure$triples <- triples
  }

  per_sample_graph <- NULL
  if (is.list(motif_obj_a$per_sample_graph) && length(motif_obj_a$per_sample_graph) &&
    is.list(motif_obj_b$per_sample_graph) && length(motif_obj_b$per_sample_graph)) {
    if (!setequal(names(motif_obj_a$per_sample_graph), samples_a) ||
      !setequal(names(motif_obj_b$per_sample_graph), samples_b)) {
      if (verbose) message("Dropping per_sample_graph (entries missing or unnamed).")
    } else {
      psa <- remap_graphs(motif_obj_a$per_sample_graph[samples_a], labels_a, lab_to_id, "motif_obj_a")
      psb <- remap_graphs(motif_obj_b$per_sample_graph[samples_b], labels_b, lab_to_id, "motif_obj_b")
      per_sample_graph <- c(psa, psb)
      per_sample_graph <- per_sample_graph[sample_name]
    }
  }

  merged <- list(
    sample_name = sample_name,
    label_levels = label_levels,
    lab_to_id = lab_to_id,
    per_sample_graph = per_sample_graph,
    raw_count = raw_count,
    exposure = exposure,
    offsets = NULL,
    norm_counts = NULL,
    relative_counts = NULL,
    offset_part_id = NULL,
    offset_part_values = NULL,
    parameters = list(
      max_edge_len = max_edge_len,
      include_wedge = include_wedge,
      offset_pseudo = offset_pseudo,
      built_from = built_from,
      erosion = erosion,
      erosion_cells_provided = erosion_cells_provided
    )
  )

  if (verbose) message("Recomputing offsets for merged object...")
  offsets_volume <- compute_offsets_volume(merged, offset_pseudo)
  offset_results <- list(volume = offsets_volume)
  offset_results <- lapply(offset_results, function(res) {
    list(
      offsets = res$offsets,
      parts = res$parts,
      part_values = res$part_values,
      norm_counts = normalize_counts_simple(raw_count, res$offsets),
      relative_counts = compute_relative_counts(raw_count, res$offsets, samples = sample_name, verbose = verbose)
    )
  })
  merged$offsets <- lapply(offset_results, `[[`, "offsets")
  merged$offset_part_id <- lapply(offset_results, `[[`, "parts")
  merged$offset_part_values <- lapply(offset_results, `[[`, "part_values")
  merged$norm_counts <- lapply(offset_results, `[[`, "norm_counts")
  merged$relative_counts <- lapply(offset_results, `[[`, "relative_counts")
  merged$parameters$offset_modes <- names(offset_results)
  merged$parameters$offset_version <- CELL_EDGER_OFFSET_VERSION
  merged$parameters$node_tmm_offsets <- isTRUE(offsets_volume$node_tmm_applied)
  merged$parameters$layer_names <- list(node = "node", edge = "edge", triangle = "triangle", wedge = "wedge")
  merged$parameters <- merged$parameters[c(
    "max_edge_len", "include_wedge", "offset_pseudo", "built_from",
    "offset_modes", "offset_version", "node_tmm_offsets", "layer_names",
    "erosion", "erosion_cells_provided"
  )]
  merged$edger <- NULL

  structure(merged, class = c("cellEdgeR_obj", "cellEdgeR_graphs"))
}

count_motifs_from_graphs <- function(graphs, max_edge_len, include_wedge, verbose, n_cores = 1,
                                     erosion = TRUE, erosion_cells = NULL) {
  if (!inherits(graphs, "cellEdgeR_graphs")) {
    stop("graphs must come from build_cell_graphs().")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package Matrix is required. Install it with install.packages('Matrix').")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    stop("n_cores must be a single positive number.")
  }
  run_parallel <- n_cores > 1 && .Platform$OS.type != "windows"
  samples <- graphs$sample_name
  lab_levels <- graphs$label_levels
  K <- length(lab_levels)
  erosion_cells <- normalize_erosion_cells(erosion_cells, samples)
  has_boundary <- any(vapply(samples, function(s) !is.null(graphs$per_sample_graph[[s]]$boundary), logical(1)))
  if (isTRUE(erosion) && isTRUE(verbose) && !has_boundary && is.null(erosion_cells)) {
    message("Erosion enabled but no boundary masks provided; counting all cells.")
  }

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
    ps <- graphs$per_sample_graph[[s]]
    edges <- filter_edges_by_threshold(ps$edges, ps$edge_len, prune_threshold)
    boundary_mask <- NULL
    if (isTRUE(erosion)) {
      boundary_mask <- coerce_boundary_mask(ps$boundary, ps$n, s, "boundary")
      if (!is.null(erosion_cells) && s %in% names(erosion_cells)) {
        extra_mask <- coerce_boundary_mask(erosion_cells[[s]], ps$n, s, "erosion_cells")
        if (!is.null(extra_mask)) {
          boundary_mask <- if (is.null(boundary_mask)) extra_mask else (boundary_mask | extra_mask)
        }
      }
    }
    if (!is.null(boundary_mask) && any(boundary_mask)) {
      keep_idx <- which(!boundary_mask)
      if (!length(keep_idx)) {
        return(list(sample_name = s, n = 0L, labels_id = integer(0), labels_chr = character(0),
          edges = matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("from", "to")))))
      }
      if (nrow(edges)) {
        keep_edges <- !boundary_mask[edges[, 1]] & !boundary_mask[edges[, 2]]
        edges <- edges[keep_edges, , drop = FALSE]
      }
      idx_map <- integer(ps$n)
      idx_map[keep_idx] <- seq_along(keep_idx)
      if (nrow(edges)) {
        edges <- matrix(idx_map[edges], ncol = 2, dimnames = list(NULL, c("from", "to")))
      } else {
        edges <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("from", "to")))
      }
      edges <- validate_edge_indices(edges, length(keep_idx), s)
      list(
        sample_name = s,
        n = length(keep_idx),
        labels_id = ps$labels_id[keep_idx],
        labels_chr = ps$labels_chr[keep_idx],
        edges = edges
      )
    } else {
      edges <- validate_edge_indices(edges, ps$n, s)
      list(sample_name = s, n = ps$n, labels_id = ps$labels_id, labels_chr = ps$labels_chr, edges = edges)
    }
  })
  names(per_sample) <- samples

  lapply_fun <- if (run_parallel) function(...) parallel::mclapply(..., mc.cores = n_cores) else lapply

  if (verbose) message("Counting node motifs (cells by label)…")
  sing_list <- lapply_fun(per_sample, function(ps) {
    if (ps$n == 0) return(integer(0))
    tab <- tabulate(ps$labels_id, nbins = K)
    names(tab) <- prefix_key(lab_levels, "N")
    tab
  })
  Y1 <- do.call(cbind, sing_list)
  dimnames(Y1) <- list(prefix_key(lab_levels, "N"), samples)
  Y1 <- Matrix::Matrix(Y1, sparse = TRUE)

  if (verbose) message("Counting edge motifs (unordered label pairs)…")
  pairs_by_sample <- lapply_fun(samples, function(s) {
    ps <- per_sample[[s]]
    e <- ps$edges
    if (length(e) == 0) {
      return(integer(0))
    }
    la <- ps$labels_chr[e[, 1]]
    lb <- ps$labels_chr[e[, 2]]
    key <- pair_key_vec(la, lb, prefix = "E")
    sort(tapply(rep(1L, length(key)), key, sum))
  })
  names(pairs_by_sample) <- samples
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
    message("Counting triangle motifs (unordered label triplets) in C++…")
    if (include_wedge) message("Also collecting wedge (open triplets)…")
  }
  tris_counts_list <- vector("list", length(samples))
  names(tris_counts_list) <- samples
  tri_totals <- numeric(length(samples))
  names(tri_totals) <- samples
  wedge_counts_list <- vector("list", length(samples))
  names(wedge_counts_list) <- samples
  wedge_totals <- numeric(length(samples))
  names(wedge_totals) <- samples
  triplet_counts_list <- vector("list", length(samples))
  names(triplet_counts_list) <- samples

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

  tri_wedge_fun <- function(s) {
    ps <- per_sample[[s]]
    if (nrow(ps$edges) == 0 || ps$n < 3) {
      return(list(tri = integer(0), tri_total = 0, wedge = integer(0), wedge_total = 0,
        triplet = integer(0)))
    }
    out <- count_triangle_labels_cpp(
      n_nodes = ps$n,
      ei = ps$edges[, 1],
      ej = ps$edges[, 2],
      labels = ps$labels_id,
      label_ids = seq_len(K),
      count_wedges = include_wedge
    )
    keys <- as.character(out$keys)
    if (length(keys)) {
      keys_lab <- format_label_triplets(keys, lab_levels, prefix = "T")
      tris <- structure(as.integer(out$counts), names = keys_lab)
    } else {
      tris <- integer(0)
    }
    tri_total <- as.numeric(out$tri_total)
    if (include_wedge) {
      wedge_keys <- as.character(out$wedge_keys)
      if (length(wedge_keys)) {
        wedge_lab <- format_label_triplets(wedge_keys, lab_levels, prefix = "W")
        wedge <- structure(as.integer(out$wedge_counts), names = wedge_lab)
      } else {
        wedge <- integer(0)
      }
      wedge_total <- as.numeric(out$wedge_total)
      triplet_keys <- as.character(out$triplet_keys)
      if (length(triplet_keys)) {
        triplet_lab <- format_label_triplets(triplet_keys, lab_levels, prefix = "T")
        triplet <- structure(as.integer(out$triplet_counts), names = triplet_lab)
      } else {
        triplet <- integer(0)
      }
    } else {
      wedge <- integer(0)
      wedge_total <- 0
      triplet <- integer(0)
    }
    list(tri = tris, tri_total = tri_total, wedge = wedge, wedge_total = wedge_total,
      triplet = triplet)
  }
  tri_res_list <- lapply_fun(samples, tri_wedge_fun)
  for (i in seq_along(samples)) {
    s <- samples[i]
    tris_counts_list[[s]] <- tri_res_list[[i]]$tri
    tri_totals[s] <- tri_res_list[[i]]$tri_total
    wedge_counts_list[[s]] <- tri_res_list[[i]]$wedge
    wedge_totals[s] <- tri_res_list[[i]]$wedge_total
    triplet_counts_list[[s]] <- tri_res_list[[i]]$triplet
  }
  tris_all_keys <- unique(unlist(lapply(tris_counts_list, names), use.names = FALSE))
  Y3 <- build_sparse_layer(tris_counts_list, tris_all_keys)
  if (include_wedge) {
    wedge_all_keys <- unique(unlist(lapply(wedge_counts_list, names), use.names = FALSE))
    Yw <- build_sparse_layer(wedge_counts_list, wedge_all_keys)
    triplet_all_keys <- unique(unlist(lapply(triplet_counts_list, names), use.names = FALSE))
    Yt <- build_sparse_layer(triplet_counts_list, triplet_all_keys)
  }

  n_edges <- vapply(per_sample, function(ps) nrow(ps$edges), 0L)
  names(n_edges) <- samples
  n_tris <- tri_totals
  volumes <- matrix(0, nrow = length(lab_levels), ncol = length(samples),
    dimnames = list(lab_levels, samples))
  center_pairs <- matrix(0, nrow = length(lab_levels), ncol = length(samples),
    dimnames = list(lab_levels, samples))
  for (s in samples) {
    ps <- per_sample[[s]]
    if (nrow(ps$edges)) {
      deg <- integer(ps$n)
      tab <- table(ps$edges)
      deg[as.integer(names(tab))] <- as.integer(tab)
      vols <- tapply(deg, ps$labels_chr, sum)
      volumes[names(vols), s] <- vols
      deg_pairs <- deg * pmax(deg - 1L, 0L)
      cp <- tapply(deg_pairs, ps$labels_chr, sum)
      center_pairs[names(cp), s] <- cp
    }
  }

  if (verbose) {
    summary_msg <- paste0(
      "Counts ready: |labels|=", K,
      ", singles=", nrow(Y1),
      ", pairs=", nrow(Y2),
      ", triangles=", nrow(Y3)
    )
    if (include_wedge) summary_msg <- paste0(summary_msg, ", wedge=", nrow(Yw))
    message(summary_msg)
  }
  counts <- list(node = Y1, edge = Y2, triangle = Y3)
  if (include_wedge) counts$wedge <- Yw
  exposure <- list(edges = n_edges, triangles = n_tris, cells = Matrix::colSums(Y1))
  exposure$volumes <- volumes
  exposure$center_pairs <- center_pairs
  if (include_wedge) {
    exposure$wedge <- wedge_totals
    exposure$triples <- Yt
  }

  list(
    sample_name = samples,
    label_levels = lab_levels,
    raw_count = counts,
    exposure = exposure,
    parameters = list(max_edge_len = max_edge_len, include_wedge = include_wedge)
  )
}

#' Build motif offsets for normalized counts and modeling
#'
#' @keywords internal
build_pair_offsets <- function(Y2, log_vols, log_2m, pseudo) {
  if (nrow(Y2) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y2)))
  }
  ab <- split_pair_labels(rownames(Y2))
  a <- ab[, 1]
  b <- ab[, 2]
  vol_a <- matrix(pseudo, nrow = length(a), ncol = ncol(log_vols), dimnames = list(rownames(Y2), colnames(log_vols)))
  vol_b <- matrix(pseudo, nrow = length(b), ncol = ncol(log_vols), dimnames = list(rownames(Y2), colnames(log_vols)))
  idx_a <- match(a, rownames(log_vols))
  idx_b <- match(b, rownames(log_vols))
  sel_a <- which(!is.na(idx_a))
  sel_b <- which(!is.na(idx_b))
  if (length(sel_a)) vol_a[sel_a, ] <- exp(log_vols[idx_a[sel_a], , drop = FALSE])
  if (length(sel_b)) vol_b[sel_b, ] <- exp(log_vols[idx_b[sel_b], , drop = FALSE])
  off <- log(pmax(vol_a, pseudo)) + log(pmax(vol_b, pseudo))
  off <- sweep(off, 2, log_2m[colnames(Y2)], FUN = "-")
  same_label <- a == b
  if (any(same_label)) {
    off[same_label, ] <- off[same_label, ] - log(2)
  }
  dimnames(off) <- dimnames(Y2)
  off
}

#' @keywords internal
build_tri_offsets <- function(Y3, log_vols, log_2m, pseudo) {
  if (nrow(Y3) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Y3)))
  }
  abc <- split_triplet_labels(rownames(Y3))
  a <- abc[, 1]
  b <- abc[, 2]
  c <- abc[, 3]
  vol_a <- matrix(pseudo, nrow = length(a), ncol = ncol(log_vols), dimnames = list(rownames(Y3), colnames(log_vols)))
  vol_b <- matrix(pseudo, nrow = length(b), ncol = ncol(log_vols), dimnames = list(rownames(Y3), colnames(log_vols)))
  vol_c <- matrix(pseudo, nrow = length(c), ncol = ncol(log_vols), dimnames = list(rownames(Y3), colnames(log_vols)))
  idx_a <- match(a, rownames(log_vols))
  idx_b <- match(b, rownames(log_vols))
  idx_c <- match(c, rownames(log_vols))
  sel_a <- which(!is.na(idx_a))
  sel_b <- which(!is.na(idx_b))
  sel_c <- which(!is.na(idx_c))
  if (length(sel_a)) vol_a[sel_a, ] <- exp(log_vols[idx_a[sel_a], , drop = FALSE])
  if (length(sel_b)) vol_b[sel_b, ] <- exp(log_vols[idx_b[sel_b], , drop = FALSE])
  if (length(sel_c)) vol_c[sel_c, ] <- exp(log_vols[idx_c[sel_c], , drop = FALSE])
  off <- log(pmax(vol_a, pseudo)) + log(pmax(vol_b, pseudo)) + log(pmax(vol_c, pseudo))
  off <- sweep(off, 2, 2 * log_2m[colnames(Y3)], FUN = "-")
  multiplicity <- vapply(seq_along(a), function(i) {
    labs <- c(a[i], b[i], c[i])
    tab <- table(labs)
    factorial(3L) / prod(factorial(as.integer(tab)))
  }, numeric(1))
  off <- off + log(multiplicity)
  dimnames(off) <- dimnames(Y3)
  off
}

build_wedge_offsets <- function(Yw, log_vols, log_2m, center_pairs, pseudo) {
  # 1. Check for empty input
  if (nrow(Yw) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(Yw)))
  }

  # 2. Split labels (A, B, C)
  abc <- split_triplet_labels(rownames(Yw))
  a <- abc[, 1]
  b <- abc[, 2]
  c <- abc[, 3]

  # 3. Map center-pair exposure for the center label, and volumes for leaves
  center_mat <- matrix(pseudo, nrow = length(a), ncol = ncol(log_vols), dimnames = list(rownames(Yw), colnames(log_vols)))
  vol_b <- matrix(pseudo, nrow = length(b), ncol = ncol(log_vols), dimnames = list(rownames(Yw), colnames(log_vols)))
  vol_c <- matrix(pseudo, nrow = length(c), ncol = ncol(log_vols), dimnames = list(rownames(Yw), colnames(log_vols)))

  idx_center <- match(a, rownames(center_pairs))
  idx_b <- match(b, rownames(log_vols))
  idx_c <- match(c, rownames(log_vols))

  sel_center <- which(!is.na(idx_center))
  sel_b <- which(!is.na(idx_b))
  sel_c <- which(!is.na(idx_c))

  if (length(sel_center)) center_mat[sel_center, ] <- center_pairs[idx_center[sel_center], , drop = FALSE]
  if (length(sel_b)) vol_b[sel_b, ] <- exp(log_vols[idx_b[sel_b], , drop = FALSE])
  if (length(sel_c)) vol_c[sel_c, ] <- exp(log_vols[idx_c[sel_c], , drop = FALSE])

  # 4. Compute Offset
  # Expected centered triples for A-B-C:
  # - If B != C: center_pairs(A) * VolB * VolC / (2m)^2
  # - If B == C: 0.5 * center_pairs(A) * VolB^2 / (2m)^2
  off <- log(pmax(center_mat, pseudo)) + log(pmax(vol_b, pseudo)) + log(pmax(vol_c, pseudo))
  off <- sweep(off, 2, 2 * log_2m[colnames(Yw)], FUN = "-")
  same_leaf <- b == c
  if (any(same_leaf)) {
    off[same_leaf, ] <- off[same_leaf, ] - log(2)
  }

  dimnames(off) <- dimnames(Yw)
  off
}

compute_offsets_volume <- function(motif_obj, pseudo) {
  samples <- motif_obj$sample_name
  Y1 <- motif_obj$raw_count$node
  Y2 <- motif_obj$raw_count$edge
  Y3 <- motif_obj$raw_count$triangle
  Yw <- motif_obj$raw_count$wedge
  log_2m <- log(pmax(2 * as.numeric(motif_obj$exposure$edges), 1))
  names(log_2m) <- samples
  log_vols <- log(pmax(motif_obj$exposure$volumes, pseudo))
  center_pairs <- motif_obj$exposure$center_pairs
  if (is.null(center_pairs)) {
    center_pairs <- motif_obj$exposure$volumes
  }
  offset2 <- build_pair_offsets(Y2, log_vols, log_2m, pseudo)
  offset3 <- build_tri_offsets(Y3, log_vols, log_2m, pseudo)
  offsetw <- if (!is.null(Yw)) build_wedge_offsets(Yw, log_vols, log_2m, center_pairs, pseudo) else NULL
  parts2 <- if (nrow(Y2)) {
    ab <- split_pair_labels(rownames(Y2))
    setNames(lapply(seq_len(nrow(Y2)), function(i) c(paste0("vol_", ab[i, 1]), paste0("vol_", ab[i, 2]), "2m")), rownames(Y2))
  } else list()
  parts3 <- if (nrow(Y3)) {
    abc <- split_triplet_labels(rownames(Y3))
    setNames(lapply(seq_len(nrow(Y3)), function(i) c(paste0("vol_", abc[i, 1]), paste0("vol_", abc[i, 2]), paste0("vol_", abc[i, 3]), "2m")), rownames(Y3))
  } else list()
  partsw <- if (!is.null(Yw) && nrow(Yw)) {
    abc <- split_triplet_labels(rownames(Yw))
    setNames(lapply(seq_len(nrow(Yw)), function(i) {
      parts <- c(paste0("center_pairs_", abc[i, 1]), paste0("vol_", abc[i, 2]), paste0("vol_", abc[i, 3]), "2m")
      if (abc[i, 2] == abc[i, 3]) parts <- c(parts, "half")
      parts
    }), rownames(Yw))
  } else list()
  log_cells <- matrix(
    log(pmax(as.numeric(motif_obj$exposure$cells), 1)),
    nrow = max(1, nrow(Y1)),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(rownames(Y1), samples)
  )
  node_tmm_applied <- FALSE
  if (!is.null(Y1) && nrow(Y1)) {
    y_nodes <- as.matrix(Y1)
    if (!identical(colnames(y_nodes), samples)) y_nodes <- y_nodes[, samples, drop = FALSE]
    dge_nodes <- edgeR::calcNormFactors(edgeR::DGEList(counts = y_nodes), method = "TMM")
    norm_factors <- dge_nodes$samples$norm.factors
    if (!is.null(norm_factors)) {
      if (!is.null(names(norm_factors))) norm_factors <- norm_factors[samples]
      if (length(norm_factors) == length(samples) &&
        all(is.finite(norm_factors)) && all(norm_factors > 0)) {
        log_cells <- sweep(log_cells, 2, log(norm_factors), "+")
        node_tmm_applied <- TRUE
      }
    }
  }
  offsets <- list(
    node = Matrix::Matrix(log_cells, sparse = TRUE, dimnames = dimnames(Y1)),
    edge = Matrix::Matrix(offset2, sparse = TRUE, dimnames = dimnames(Y2)),
    triangle = Matrix::Matrix(offset3, sparse = TRUE, dimnames = dimnames(Y3)),
    wedge = if (!is.null(offsetw)) Matrix::Matrix(offsetw, sparse = TRUE, dimnames = dimnames(Yw)) else NULL
  )
  list(
    offsets = offsets,
    parts = list(edge = parts2, triangle = parts3, wedge = partsw),
    part_values = list(
      volumes = motif_obj$exposure$volumes,
      center_pairs = center_pairs,
      two_m = 2 * as.numeric(motif_obj$exposure$edges),
      half = 0.5
    ),
    node_tmm_applied = node_tmm_applied
  )
}

fit_edge_posterior <- function(Y2, offset_edge, design) {
  dge <- edgeR::DGEList(counts = Y2)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge$offset <- as.matrix(offset_edge)
  fit <- tryCatch({
    edgeR::estimateGLMRobustDisp(dge, design = design)
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  fit <- tryCatch(edgeR::glmQLFit(fit, design = design), error = function(e) NULL)
  fit
}

build_hier_offsets_from_mu <- function(motif_obj, pseudo, edge_mu) {
  Y1 <- motif_obj$raw_count$node
  Y2 <- motif_obj$raw_count$edge
  Y3 <- motif_obj$raw_count$triangle
  Yw <- motif_obj$raw_count$wedge
  samples <- motif_obj$sample_name
  log_2m <- log(pmax(2 * as.numeric(motif_obj$exposure$edges), 1))
  names(log_2m) <- samples
  log_vols <- log(pmax(motif_obj$exposure$volumes, pseudo))
  center_pairs <- motif_obj$exposure$center_pairs
  if (is.null(center_pairs)) {
    center_pairs <- motif_obj$exposure$volumes
  }
  offset2 <- build_pair_offsets(Y2, log_vols, log_2m, pseudo) # edges stay volume-based
  tri_from_mu <- function(Y, edge_mat) {
    if (is.null(Y) || nrow(Y) == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = ncol(edge_mat), dimnames = list(character(), colnames(edge_mat))))
    }
    abc <- split_triplet_labels(rownames(Y))
    get_mu <- function(keys) {
      idx <- match(keys, rownames(edge_mat))
      out <- matrix(pseudo, nrow = length(keys), ncol = ncol(edge_mat), dimnames = list(keys, colnames(edge_mat)))
      sel <- which(!is.na(idx))
      if (length(sel)) out[sel, ] <- as.matrix(edge_mat[idx[sel], , drop = FALSE])
      out
    }
    ABk <- pair_key_vec(abc[, 1], abc[, 2], prefix = "E")
    ACk <- pair_key_vec(abc[, 1], abc[, 3], prefix = "E")
    BCk <- pair_key_vec(abc[, 2], abc[, 3], prefix = "E")
    AB <- get_mu(ABk); AC <- get_mu(ACk); BC <- get_mu(BCk)
    off <- log(pmax(AB, pseudo)) + log(pmax(AC, pseudo)) + log(pmax(BC, pseudo))
    off <- sweep(off, 2, 2 * log_2m[colnames(Y)], FUN = "-")
    dimnames(off) <- dimnames(Y)
    off
  }
  offset3 <- tri_from_mu(Y3, edge_mu)
  offsetw <- if (!is.null(Yw)) build_wedge_offsets(Yw, log_vols, log_2m, center_pairs, pseudo) else NULL
  parts2 <- if (nrow(Y2)) {
    ab <- split_pair_labels(rownames(Y2))
    setNames(lapply(seq_len(nrow(Y2)), function(i) c(paste0("vol_", ab[i, 1]), paste0("vol_", ab[i, 2]), "2m")), rownames(Y2))
  } else list()
  parts3 <- if (nrow(Y3)) {
    abc <- split_triplet_labels(rownames(Y3))
    setNames(lapply(seq_len(nrow(Y3)), function(i) c(pair_key_vec(abc[i, 1], abc[i, 2], "E"), pair_key_vec(abc[i, 1], abc[i, 3], "E"), pair_key_vec(abc[i, 2], abc[i, 3], "E"), "2m")), rownames(Y3))
  } else list()
  partsw <- if (!is.null(Yw) && nrow(Yw)) {
    abc <- split_triplet_labels(rownames(Yw))
    setNames(lapply(seq_len(nrow(Yw)), function(i) {
      parts <- c(paste0("center_pairs_", abc[i, 1]), paste0("vol_", abc[i, 2]), paste0("vol_", abc[i, 3]), "2m")
      if (abc[i, 2] == abc[i, 3]) parts <- c(parts, "half")
      parts
    }), rownames(Yw))
  } else list()
  log_cells <- matrix(
    log(pmax(as.numeric(motif_obj$exposure$cells), 1)),
    nrow = max(1, nrow(Y1)),
    ncol = length(samples),
    byrow = TRUE,
    dimnames = list(rownames(Y1), samples)
  )
  if (!is.null(Y1) && nrow(Y1)) {
    y_nodes <- as.matrix(Y1)
    if (!identical(colnames(y_nodes), samples)) y_nodes <- y_nodes[, samples, drop = FALSE]
    dge_nodes <- edgeR::calcNormFactors(edgeR::DGEList(counts = y_nodes), method = "TMM")
    norm_factors <- dge_nodes$samples$norm.factors
    if (!is.null(norm_factors)) {
      if (!is.null(names(norm_factors))) norm_factors <- norm_factors[samples]
      if (length(norm_factors) == length(samples) &&
        all(is.finite(norm_factors)) && all(norm_factors > 0)) {
        log_cells <- sweep(log_cells, 2, log(norm_factors), "+")
      }
    }
  }
  offsets <- list(
    node = Matrix::Matrix(log_cells, sparse = TRUE, dimnames = dimnames(Y1)),
    edge = Matrix::Matrix(offset2, sparse = TRUE, dimnames = dimnames(Y2)),
    triangle = Matrix::Matrix(offset3, sparse = TRUE, dimnames = dimnames(Y3)),
    wedge = if (!is.null(offsetw)) Matrix::Matrix(offsetw, sparse = TRUE, dimnames = dimnames(Yw)) else NULL
  )
  list(
    offsets = offsets,
    parts = list(edge = parts2, triangle = parts3, wedge = partsw),
    part_values = list(
      two_m = 2 * as.numeric(motif_obj$exposure$edges),
      edge_mu = edge_mu,
      center_pairs = center_pairs,
      half = 0.5
    )
  )
}

compute_offsets_hierarchical <- function(motif_obj, pseudo, mode = c("null", "full"), design = NULL) {
  mode <- match.arg(mode)
  # Fit edge layer to get posterior means
  volume_offsets <- compute_offsets_volume(motif_obj, pseudo)$offsets
  design_null <- matrix(1, nrow = length(motif_obj$sample_name), ncol = 1)
  colnames(design_null) <- "Intercept"
  des <- if (mode == "null") design_null else design
  fit <- fit_edge_posterior(motif_obj$raw_count$edge, volume_offsets$edge, des)
  if (is.null(fit)) {
    if (mode == "null") warning("edgeR fit failed for edge-derived offsets (null); falling back to counts.")
    else warning("edgeR fit failed for edge-derived offsets (full); falling back to counts.")
    edge_mu <- as.matrix(motif_obj$raw_count$edge)
  } else {
    edge_mu <- fit$fitted.values
  }
  build_hier_offsets_from_mu(motif_obj, pseudo, edge_mu)
}

normalize_counts_simple <- function(counts, offsets) {
  norm_layer <- function(Y, off) {
    if (is.null(Y) || nrow(Y) == 0) {
      return(Matrix::Matrix(0, nrow = 0, ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    Matrix::Matrix(as.matrix(Y) / exp(as.matrix(off)), sparse = TRUE, dimnames = dimnames(Y))
  }
  list(
    node = norm_layer(counts$node, offsets$node),
    edge = norm_layer(counts$edge, offsets$edge),
    triangle = norm_layer(counts$triangle, offsets$triangle),
    wedge = if (!is.null(counts$wedge) && !is.null(offsets$wedge)) norm_layer(counts$wedge, offsets$wedge) else NULL
  )
}

compute_relative_counts <- function(counts, offsets, samples, verbose = FALSE, eps = 1) {
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
    dge$samples$norm.factors <- 1
    dge$offset <- as.matrix(off)
    dge <- suppressWarnings(tryCatch(edgeR::estimateGLMRobustDisp(dge, design = design_null), error = function(e) NULL))
    if (is.null(dge)) {
      warning("edgeR dispersion estimation failed for model residuals; returning zeros.")
      return(Matrix::Matrix(0, nrow = nrow(Y), ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    fit <- suppressWarnings(tryCatch(edgeR::glmQLFit(dge, design = design_null), error = function(e) NULL))
    if (is.null(fit)) {
      warning("edgeR fit failed for model residuals; returning zeros.")
      return(Matrix::Matrix(0, nrow = nrow(Y), ncol = ncol(Y), sparse = TRUE, dimnames = dimnames(Y)))
    }
    mu <- fit$fitted.values
    obs <- as.matrix(Y)
    Matrix::Matrix(log2(obs + eps) - log2(mu + eps), sparse = TRUE, dimnames = dimnames(Y))
  }

  list(
    node = fit_layer(counts$node, offsets$node),
    edge = fit_layer(counts$edge, offsets$edge),
    triangle = fit_layer(counts$triangle, offsets$triangle),
    wedge = if (!is.null(counts$wedge) && !is.null(offsets$wedge)) fit_layer(counts$wedge, offsets$wedge) else NULL
  )
}

#' Retrieve normalized counts as a single data frame
#'
#' Pulls all normalized motifs from a cellgraph (`cellEdgeR_obj`) into one data frame by row-binding layers.
#'
#' @param cellgraph A `cellEdgeR_obj` returned by [count_motifs_graphs()].
#' @param offset_mode Which offset mode to pull normalized counts from; defaults to the first available when multiple are stored.
#' @param log2 Logical; if `TRUE` (default), apply `log2(count + pseudocount)` to the normalized counts.
#' @param pseudocount Small positive value added before log transformation (default 1 to match offset construction).
#' @return A matrix/data frame with samples as rownames and motifs as column names.
#' @export
get_norm_counts <- function(cellgraph, offset_mode = NULL, log2 = TRUE, pseudocount = 1) {
  validate_motif_obj(cellgraph, require_offsets = FALSE)
  norm_container <- cellgraph$norm_counts
  if (is.null(norm_container) || !is.list(norm_container)) {
    stop("cellgraph is missing norm_counts; run count_motifs_graphs() first.")
  }
  # Handle legacy single-layer vs multi-offset storage
  if (all(c("node", "edge", "triangle") %in% names(norm_container))) {
    mats <- norm_container
  } else {
    avail <- names(norm_container)
    if (is.null(offset_mode)) offset_mode <- avail[1]
    if (!offset_mode %in% avail) stop("offset_mode not found in norm_counts; available: ", paste(avail, collapse = ", "))
    mats <- norm_container[[offset_mode]]
  }
  mats <- Filter(function(m) !is.null(m) && nrow(m) > 0, mats)
  if (!length(mats)) return(data.frame())
  transform_mat <- function(m) {
    mat <- as.matrix(m)
    if (log2) mat <- base::log2(mat + pseudocount)
    mat
  }
  mats <- lapply(mats, transform_mat)
  combined <- do.call(rbind, mats) # motifs in rows, samples in columns
  out <- t(combined)               # samples in rows, motifs in columns
  data.frame(out, check.names = FALSE)
}

#' Retrieve motif values (raw, normalized)
#'
#' Pull raw counts or normalized counts for a motif (or its lower-order submotifs)
#' into a samples-by-motifs data frame.
#'
#' @param cellgraph A `cellEdgeR_obj` returned by [count_motifs_graphs()].
#' @param motif_key Character vector of motif keys (e.g., `E_A_B`, `T_A_B_C`). If `NULL`, return all motifs.
#' @param include_submotifs Logical; if `TRUE`, include lower-order submotifs implied by the motif key
#'   (edges for triangles/wedges, nodes for edges).
#' @param value Which value to return: `raw` or `norm`.
#' @return A data frame with sample names as rownames and one column per motif.
#'   When `include_submotifs = TRUE`, requested motifs appear first, followed by unique submotifs.
#' @details
#' Normalized values are computed directly from the volume offsets so they match the offsets used by
#' [motif_edger()] (including node TMM adjustments when present).
#' @export
get_motif_values <- function(
  cellgraph,
  motif_key = NULL,
  include_submotifs = FALSE,
  value = c("raw", "norm")
) {
  validate_motif_obj(cellgraph, require_offsets = FALSE)
  samples <- cellgraph$sample_name
  value <- match.arg(value)

  raw_container <- cellgraph$raw_count
  if (is.null(raw_container) || !is.list(raw_container)) {
    stop("cellgraph is missing raw_count; run count_motifs_graphs() first.")
  }
  motif_type_map <- character()
  all_motifs <- character()
  for (layer in names(raw_container)) {
    mat <- raw_container[[layer]]
    if (is.null(mat)) next
    keys <- rownames(mat)
    if (is.null(keys)) stop("raw_count$", layer, " is missing row names.")
    if (!length(keys)) next
    all_motifs <- c(all_motifs, keys)
    motif_type_map <- c(motif_type_map, stats::setNames(rep(layer, length(keys)), keys))
  }
  all_motifs <- unique(all_motifs)
  if (!length(all_motifs)) {
    out <- data.frame(row.names = cellgraph$sample_name)
    return(out)
  }

  if (!is.null(motif_key)) {
    motif_key <- as.character(motif_key)
    motif_key <- motif_key[!duplicated(motif_key)]
    missing <- setdiff(motif_key, all_motifs)
    if (length(missing)) {
      stop("motif_key not found in raw_count: ", paste(missing, collapse = ", "))
    }
  }
  if (include_submotifs && is.null(motif_key)) {
    stop("include_submotifs requires motif_key.")
  }

  needs_offsets <- value == "norm"
  if (needs_offsets) {
    modes <- available_offset_modes(cellgraph)
    if (!length(modes)) stop("cellgraph is missing offsets; run count_motifs_graphs() first.")
    if (!"volume" %in% modes) {
      stop("cellgraph is missing volume offsets; run count_motifs_graphs() first.")
    }
  }

  if (is.null(motif_key)) {
    ordered_motifs <- all_motifs
  } else if (!include_submotifs) {
    ordered_motifs <- motif_key
  } else {
    unique_keep <- function(x) x[!duplicated(x)]
    submotifs_for_key <- function(key) {
      layer <- motif_type_map[[key]]
      if (is.null(layer) || is.na(layer)) return(character())
      labels <- strsplit(strip_prefix(key), "_", fixed = TRUE)[[1]]
      if (layer %in% c("triangle", "wedge")) {
        if (length(labels) != 3L) return(character())
        a <- labels[1]
        b <- labels[2]
        c <- labels[3]
        c(
          pair_key_vec(a, b, prefix = "E"),
          pair_key_vec(a, c, prefix = "E"),
          pair_key_vec(b, c, prefix = "E")
        )
      } else if (layer == "edge") {
        if (length(labels) != 2L) return(character())
        prefix_key(labels, "N")
      } else {
        character()
      }
    }
    submotifs <- unique_keep(unlist(lapply(motif_key, submotifs_for_key), use.names = FALSE))
    submotifs <- intersect(submotifs, all_motifs)
    ordered_motifs <- c(motif_key, setdiff(submotifs, motif_key))
  }

  if (!length(ordered_motifs)) {
    out <- data.frame(row.names = cellgraph$sample_name)
    return(out)
  }

  select_mode_container <- function(container, label) {
    if (is.null(container) || !is.list(container)) {
      stop("cellgraph is missing ", label, "; run count_motifs_graphs() first.")
    }
    if (all(c("node", "edge", "triangle") %in% names(container))) {
      return(list(mats = container, mode = "volume"))
    }
    if (!"volume" %in% names(container)) {
      stop(label, " does not contain volume offsets.")
    }
    list(mats = container[["volume"]], mode = "volume")
  }

  tmm_baked <- isTRUE(cellgraph$parameters$node_tmm_offsets)

  matrix_from_layers <- function(container, motifs, samples, fill = NA_real_) {
    mat <- matrix(fill, nrow = length(motifs), ncol = length(samples),
      dimnames = list(motifs, samples))
    for (layer in names(container)) {
      src <- container[[layer]]
      if (is.null(src) || !nrow(src)) next
      if (is.null(rownames(src))) stop("Container layer is missing row names.")
      if (!identical(colnames(src), samples)) {
        if (setequal(colnames(src), samples)) {
          src <- src[, samples, drop = FALSE]
        } else {
          stop("Container layer columns do not match sample names.")
        }
      }
      hit <- intersect(rownames(src), motifs)
      if (!length(hit)) next
      sub <- as.matrix(src[hit, samples, drop = FALSE])
      mat[hit, ] <- sub
    }
    mat
  }

  get_offsets_for_mode <- function() {
    sel <- select_mode_container(cellgraph$offsets, "offsets")
    offsets_layers <- sel$mats
    if (tmm_baked) return(offsets_layers)
    if (is.null(offsets_layers$node) || is.null(raw_container$node) || nrow(raw_container$node) == 0) {
      return(offsets_layers)
    }
    node_counts <- raw_container$node
    if (is.null(colnames(node_counts))) stop("raw_count$node is missing column names.")
    if (!identical(colnames(node_counts), samples)) {
      if (setequal(colnames(node_counts), samples)) {
        node_counts <- node_counts[, samples, drop = FALSE]
      } else {
        stop("raw_count$node columns do not match sample names.")
      }
    }
    node_offsets <- offsets_layers$node
    if (is.null(colnames(node_offsets))) stop("offsets$node is missing column names.")
    if (!identical(colnames(node_offsets), samples)) {
      if (setequal(colnames(node_offsets), samples)) {
        node_offsets <- node_offsets[, samples, drop = FALSE]
      } else {
        stop("offsets$node columns do not match sample names.")
      }
    }
    dge_nodes <- edgeR::calcNormFactors(edgeR::DGEList(counts = as.matrix(node_counts)), method = "TMM")
    norm_factors <- dge_nodes$samples$norm.factors
    if (!is.null(names(norm_factors))) norm_factors <- norm_factors[samples]
    if (length(norm_factors) == length(samples) &&
      all(is.finite(norm_factors)) && all(norm_factors > 0)) {
      node_mat <- as.matrix(node_offsets)
      node_mat <- sweep(node_mat, 2, log(norm_factors), "+")
      offsets_layers$node <- Matrix::Matrix(node_mat, sparse = TRUE, dimnames = dimnames(node_mat))
    }
    offsets_layers
  }

  value_mat <- switch(
    value,
    raw = matrix_from_layers(raw_container, ordered_motifs, samples, fill = 0),
    norm = {
      offsets_layers <- get_offsets_for_mode()
      offset_mat <- matrix_from_layers(offsets_layers, ordered_motifs, samples, fill = NA_real_)
      raw_mat <- matrix_from_layers(raw_container, ordered_motifs, samples, fill = 0)
      raw_mat / exp(offset_mat)
    }
  )

  out <- as.data.frame(t(value_mat), check.names = FALSE)
  rownames(out) <- samples
  out
}

normalize_triplet_key <- function(keys, prefix) {
  if (!length(keys)) return(character(0))
  parts <- strsplit(strip_prefix(keys), "_", fixed = TRUE)
  sorted <- vapply(parts, function(x) paste(sort(x), collapse = "_"), character(1))
  paste0(prefix, "_", sorted)
}

aggregate_rows_by_group <- function(mat, groups, n_groups, group_names, samples) {
  if (is.null(mat) || !nrow(mat) || !length(groups)) {
    return(Matrix::Matrix(0, nrow = n_groups, ncol = length(samples),
      sparse = TRUE, dimnames = list(group_names, samples)))
  }
  keep <- !is.na(groups)
  if (!any(keep)) {
    return(Matrix::Matrix(0, nrow = n_groups, ncol = length(samples),
      sparse = TRUE, dimnames = list(group_names, samples)))
  }
  agg <- Matrix::sparseMatrix(
    i = groups[keep],
    j = which(keep),
    x = 1,
    dims = c(n_groups, length(groups))
  )
  out <- agg %*% mat
  rownames(out) <- group_names
  out
}

build_triplet_map <- function(tri_keys, wedge_keys, prefix = "TP") {
  tri_norm <- normalize_triplet_key(tri_keys, prefix)
  wedge_norm <- normalize_triplet_key(wedge_keys, prefix)
  triplet_keys <- sort(unique(c(tri_norm, wedge_norm)))
  if (!length(triplet_keys)) {
    return(list(
      keys = character(),
      tri_group = integer(0),
      wedge_group = integer(0),
      n_components = integer(0),
      rep_source = character(0),
      rep_index = integer(0)
    ))
  }
  tri_group <- if (length(tri_norm)) match(tri_norm, triplet_keys) else integer(0)
  wedge_group <- if (length(wedge_norm)) match(wedge_norm, triplet_keys) else integer(0)

  n_components <- integer(length(triplet_keys))
  if (length(tri_group)) n_components <- n_components + tabulate(tri_group, nbins = length(triplet_keys))
  if (length(wedge_group)) n_components <- n_components + tabulate(wedge_group, nbins = length(triplet_keys))

  rep_tri <- rep(NA_integer_, length(triplet_keys))
  if (length(tri_group)) {
    idx <- tapply(seq_along(tri_group), tri_group, function(v) v[1])
    rep_tri[as.integer(names(idx))] <- as.integer(idx)
  }
  rep_wedge <- rep(NA_integer_, length(triplet_keys))
  if (length(wedge_group)) {
    idx <- tapply(seq_along(wedge_group), wedge_group, function(v) v[1])
    rep_wedge[as.integer(names(idx))] <- as.integer(idx)
  }

  rep_source <- ifelse(!is.na(rep_tri), "triangle", "wedge")
  rep_index <- ifelse(!is.na(rep_tri), rep_tri, rep_wedge)
  if (anyNA(rep_index)) {
    stop("Triplet merge failed to identify representative rows.")
  }

  list(
    keys = triplet_keys,
    tri_group = tri_group,
    wedge_group = wedge_group,
    n_components = n_components,
    rep_source = rep_source,
    rep_index = rep_index
  )
}

collapse_triplet_counts <- function(tri_mat, wedge_mat, map, samples) {
  if (!length(map$keys)) {
    return(Matrix::Matrix(0, nrow = 0, ncol = length(samples),
      sparse = TRUE, dimnames = list(character(), samples)))
  }
  tri_sum <- aggregate_rows_by_group(tri_mat, map$tri_group, length(map$keys), map$keys, samples)
  wedge_sum <- aggregate_rows_by_group(wedge_mat, map$wedge_group, length(map$keys), map$keys, samples)
  tri_sum + wedge_sum
}

collapse_triplet_offsets <- function(tri_off, wedge_off, map, samples) {
  if (!length(map$keys)) {
    return(matrix(numeric(0), nrow = 0, ncol = length(samples),
      dimnames = list(character(), samples)))
  }
  tri_mu <- if (!is.null(tri_off) && nrow(tri_off)) {
    aggregate_rows_by_group(exp(as.matrix(tri_off)), map$tri_group, length(map$keys), map$keys, samples)
  } else {
    Matrix::Matrix(0, nrow = length(map$keys), ncol = length(samples),
      sparse = TRUE, dimnames = list(map$keys, samples))
  }
  wedge_mu <- if (!is.null(wedge_off) && nrow(wedge_off)) {
    aggregate_rows_by_group(exp(as.matrix(wedge_off)), map$wedge_group, length(map$keys), map$keys, samples)
  } else {
    Matrix::Matrix(0, nrow = length(map$keys), ncol = length(samples),
      sparse = TRUE, dimnames = list(map$keys, samples))
  }
  mu <- tri_mu + wedge_mu
  log_mu <- log(as.matrix(mu))
  dimnames(log_mu) <- list(map$keys, samples)
  log_mu
}

sum_wedges_to_triangles <- function(wedge_mat, tri_keys, samples) {
  if (is.null(wedge_mat) || !nrow(wedge_mat) || !length(tri_keys)) {
    return(Matrix::Matrix(0, nrow = length(tri_keys), ncol = length(samples),
      sparse = TRUE, dimnames = list(tri_keys, samples)))
  }
  wedge_keys <- normalize_triplet_key(rownames(wedge_mat), prefix = "T")
  groups <- match(wedge_keys, tri_keys)
  keep <- !is.na(groups)
  if (!any(keep)) {
    return(Matrix::Matrix(0, nrow = length(tri_keys), ncol = length(samples),
      sparse = TRUE, dimnames = list(tri_keys, samples)))
  }
  agg <- Matrix::sparseMatrix(
    i = groups[keep],
    j = which(keep),
    x = 1,
    dims = c(length(tri_keys), nrow(wedge_mat))
  )
  out <- agg %*% wedge_mat
  rownames(out) <- tri_keys
  out
}


#' Differential motif testing with edgeR
#'
#' Fit a single volume-offset edgeR QL model across edge, wedge, and triangle
#' motifs. Node counts are used internally to build graph volumes but are not
#' tested by `motif_edger()`. Wedges and triangles are kept as separate motifs;
#' there is no merged triplet test. Low-count motifs are filtered with
#' [edgeR::filterByExpr()] before model fitting, so motifs carried by too few
#' samples are not tested.
#'
#' @param cellgraph Output list from [count_motifs_graphs()].
#' @param sample_df Data frame with sample metadata; rownames must match `cellgraph$sample_name`.
#' @param design_formula Formula or formula string passed to `model.matrix`, e.g. `~ condition + batch`.
#' @param verbose Logical; print progress.
#'
#' @return The input `cellgraph` augmented with `edger`, containing one stored
#'   strategy named `volume`. Use [top_edges()] for edge results and
#'   [top_motifs2()] for wedge/triangle results with edge-submotif statistics.
#'   Filter details are stored in `cellgraph$edger$filter`.
#' @export
motif_edger <- function(
  cellgraph,
  sample_df,
  design_formula,
  verbose = TRUE
) {
  validate_motif_obj(cellgraph, require_offsets = TRUE)
  if (!"Matrix" %in% loadedNamespaces()) {
    requireNamespace("Matrix")
  }
  samples <- cellgraph$sample_name
  if (is.null(rownames(sample_df))) stop("sample_df must have rownames = sample names")
  missing_rows <- setdiff(samples, rownames(sample_df))
  if (length(missing_rows)) {
    stop("sample_df is missing rows for: ", paste(missing_rows, collapse = ", "))
  }
  sample_df <- as.data.frame(sample_df[samples, , drop = FALSE])
  sample_df[] <- lapply(sample_df, function(col) if (is.character(col)) factor(col) else col)

  formula_obj <- stats::as.formula(design_formula)
  used_vars <- all.vars(formula_obj)
  if ("." %in% used_vars) used_vars <- names(sample_df)
  if (length(used_vars)) {
    single_level <- used_vars[vapply(sample_df[, used_vars, drop = FALSE], function(col) {
      if (is.factor(col) || is.character(col)) length(unique(as.character(col))) < 2 else FALSE
    }, logical(1))]
    if (length(single_level)) {
      stop("Design variable(s) have fewer than 2 levels after matching samples: ",
        paste(single_level, collapse = ", "),
        ". Provide at least two groups or use '~ 1' for intercept-only.")
    }
  }

  design <- stats::model.matrix(formula_obj, data = sample_df)
  if (ncol(design) == 0) stop("Design formula produced an empty model matrix.")
  if (is.null(colnames(design))) colnames(design) <- paste0("coef", seq_len(ncol(design)))
  design_null <- stats::model.matrix(~ 1, data = sample_df)
  if (is.null(colnames(design_null))) colnames(design_null) <- "(Intercept)"

  counts <- cellgraph$raw_count
  align_counts <- function(mat, layer) {
    if (is.null(mat) || nrow(mat) == 0) return(NULL)
    if (is.null(colnames(mat))) stop("raw_count$", layer, " is missing column names.")
    if (!identical(colnames(mat), samples)) {
      if (!setequal(colnames(mat), samples)) stop("raw_count$", layer, " columns do not match sample names.")
      mat <- mat[, samples, drop = FALSE]
    }
    mat
  }

  counts_layers <- list(
    edge = align_counts(counts$edge, "edge"),
    triangle = align_counts(counts$triangle, "triangle"),
    wedge = align_counts(counts$wedge, "wedge")
  )

  get_volume_offsets <- function() {
    offs <- cellgraph$offsets
    if (is.null(offs) || !is.list(offs) || !length(offs)) {
      stop("Volume offsets are missing; rerun count_motifs_graphs().")
    }
    if (all(c("node", "edge", "triangle") %in% names(offs))) return(offs)
    if (!"volume" %in% names(offs)) stop("Volume offsets are missing; rerun count_motifs_graphs().")
    offs$volume
  }
  offsets_volume <- get_volume_offsets()
  align_offsets <- function(counts_mat, offset_mat, layer) {
    if (is.null(counts_mat) || nrow(counts_mat) == 0) return(NULL)
    if (is.null(offset_mat)) stop("Offsets missing for layer: ", layer)
    if (is.null(colnames(offset_mat))) stop("Offsets for layer ", layer, " are missing column names.")
    if (!identical(colnames(offset_mat), samples)) {
      if (!setequal(colnames(offset_mat), samples)) stop("Offsets for layer ", layer, " columns do not match sample names.")
      offset_mat <- offset_mat[, samples, drop = FALSE]
    }
    if (!identical(rownames(offset_mat), rownames(counts_mat))) {
      if (!setequal(rownames(offset_mat), rownames(counts_mat))) {
        stop("Offsets row names do not match counts for layer: ", layer)
      }
      offset_mat <- offset_mat[rownames(counts_mat), , drop = FALSE]
    }
    offset_mat
  }
  offset_layers <- list(
    edge = align_offsets(counts_layers$edge, offsets_volume$edge, "edge"),
    triangle = align_offsets(counts_layers$triangle, offsets_volume$triangle, "triangle"),
    wedge = align_offsets(counts_layers$wedge, offsets_volume$wedge, "wedge")
  )

  model_layers <- c("edge", "triangle", "wedge")
  mats <- Filter(Negate(is.null), counts_layers[model_layers])
  off_mats <- Filter(Negate(is.null), offset_layers[model_layers])
  if (length(mats) != length(off_mats)) {
    stop("Counts and offsets are not aligned across modeled motif layers.")
  }
  Y_all <- if (length(mats)) do.call(rbind, mats) else NULL
  offsets_all <- if (length(off_mats)) do.call(rbind, off_mats) else NULL
  motif_type <- unlist(lapply(names(mats), function(layer) rep(layer, nrow(mats[[layer]]))), use.names = FALSE)
  motif_info <- data.frame(
    motif = if (is.null(Y_all)) character() else rownames(Y_all),
    motif_type = motif_type,
    stringsAsFactors = FALSE
  )
  filter_info <- list(
    method = "edgeR::filterByExpr",
    min_count = 10,
    min_total_count = 15,
    n_total = if (is.null(Y_all)) 0L else nrow(Y_all),
    n_tested = if (is.null(Y_all)) 0L else nrow(Y_all),
    n_filtered = 0L,
    filtered_motifs = character()
  )

  fit_model <- function(dge_base, design_mat, label) {
    dge_fit <- tryCatch(edgeR::estimateGLMRobustDisp(dge_base, design = design_mat), error = function(err) NULL)
    if (is.null(dge_fit)) {
      dge_fit <- tryCatch(edgeR::estimateDisp(dge_base, design = design_mat, robust = TRUE), error = function(err) NULL)
    }
    if (is.null(dge_fit)) {
      if (verbose) message("edgeR dispersion estimation failed for ", label, " model; tests will be empty.")
      return(list(dge = dge_base, fit = NULL, tests = list()))
    }
    fit <- tryCatch(edgeR::glmQLFit(dge_fit, design = design_mat, robust = TRUE), error = function(err) NULL)
    if (is.null(fit)) {
      if (verbose) message("edgeR GLM fit failed for ", label, " model; tests will be empty.")
      return(list(dge = dge_fit, fit = NULL, tests = list()))
    }
    tests <- lapply(seq_len(ncol(design_mat)), function(i) {
      tryCatch(edgeR::glmQLFTest(fit, coef = i), error = function(err) NULL)
    })
    names(tests) <- colnames(design_mat)
    list(dge = dge_fit, fit = fit, tests = tests)
  }

  if (is.null(Y_all) || !nrow(Y_all)) {
    empty_strategy <- list(
      type = "edgeR",
      offset_mode = "volume",
      model = "volume_chung_lu",
      design = design,
      design_formula = design_formula,
      full = list(design = design, design_formula = design_formula, dge = NULL, fit = NULL, tests = list()),
      null = list(design = design_null, design_formula = "~ 1", dge = NULL, fit = NULL, tests = list())
    )
    cellgraph$edger <- list(
      strategies = list(volume = empty_strategy),
      motif_info = motif_info,
      sample_df = sample_df,
      filter = filter_info,
      result_sets = c("edges", "motifs2"),
      statistics = "volume_chung_lu",
      offset_version = cellgraph$parameters$offset_version
    )
    return(cellgraph)
  }

  dge_base <- edgeR::DGEList(counts = Y_all)
  dge_base$samples$norm.factors <- 1
  keep <- tryCatch(
    edgeR::filterByExpr(dge_base, design = design),
    error = function(err) {
      warning("edgeR::filterByExpr failed; testing all motifs. Reason: ", conditionMessage(err))
      rep(TRUE, nrow(Y_all))
    }
  )
  if (!is.logical(keep) || length(keep) != nrow(Y_all)) {
    warning("edgeR::filterByExpr returned an invalid filter; testing all motifs.")
    keep <- rep(TRUE, nrow(Y_all))
  }
  keep[is.na(keep)] <- FALSE
  names(keep) <- rownames(Y_all)
  filter_info$n_tested <- sum(keep)
  filter_info$n_filtered <- sum(!keep)
  filter_info$filtered_motifs <- rownames(Y_all)[!keep]

  if (verbose) {
    message(
      "Filtering low-count motifs with edgeR::filterByExpr: kept ",
      filter_info$n_tested, " / ", filter_info$n_total, "."
    )
  }

  if (!any(keep)) {
    empty_strategy <- list(
      type = "edgeR",
      offset_mode = "volume",
      model = "volume_chung_lu",
      design = design,
      design_formula = design_formula,
      full = list(design = design, design_formula = design_formula, dge = NULL, fit = NULL, tests = list()),
      null = list(design = design_null, design_formula = "~ 1", dge = NULL, fit = NULL, tests = list())
    )
    cellgraph$edger <- list(
      strategies = list(volume = empty_strategy),
      motif_info = motif_info[keep, , drop = FALSE],
      sample_df = sample_df,
      filter = filter_info,
      result_sets = c("edges", "motifs2"),
      statistics = "volume_chung_lu",
      offset_version = cellgraph$parameters$offset_version
    )
    return(cellgraph)
  }

  Y_all <- Y_all[keep, , drop = FALSE]
  offsets_all <- offsets_all[keep, , drop = FALSE]
  motif_info <- motif_info[keep, , drop = FALSE]

  if (verbose) message("Fitting edgeR (QL) for volume offsets on edges, wedges, and triangles...")
  dge_base <- edgeR::DGEList(counts = Y_all)
  dge_base$samples$norm.factors <- 1
  dge_base$offset <- as.matrix(offsets_all)

  full_res <- fit_model(dge_base, design, "volume/full")
  null_res <- if (identical(design, design_null)) full_res else fit_model(dge_base, design_null, "volume/null")

  cellgraph$edger <- list(
    strategies = list(
      volume = list(
        type = "edgeR",
        offset_mode = "volume",
        model = "volume_chung_lu",
        design = design,
        design_formula = design_formula,
        full = full_res,
        null = null_res
      )
    ),
    motif_info = motif_info,
    sample_df = sample_df,
    filter = filter_info,
    result_sets = c("edges", "motifs2"),
    statistics = "volume_chung_lu",
    offset_version = cellgraph$parameters$offset_version
  )
  cellgraph
}
infer_motif_type <- function(motifs) {
  out <- rep(NA_character_, length(motifs))
  out[startsWith(motifs, "N_")] <- "node"
  out[startsWith(motifs, "E_")] <- "edge"
  out[startsWith(motifs, "T_")] <- "triangle"
  out[startsWith(motifs, "W_")] <- "wedge"
  out[startsWith(motifs, "TP_")] <- "triplet"
  out
}

resolve_coef_name <- function(design_mat, coef, coef_names) {
  if (is.null(coef_names)) coef_names <- colnames(design_mat)
  if (is.null(coef_names)) coef_names <- paste0("coef", seq_len(ncol(design_mat)))
  if (is.null(coef)) {
    non_int <- which(coef_names != "(Intercept)")
    if (length(non_int)) return(coef_names[non_int[1]])
    return(coef_names[1])
  }
  if (is.character(coef)) {
    if (!coef %in% coef_names) stop("coef name not found in model design.")
    return(coef)
  }
  if (!is.numeric(coef) || length(coef) != 1 || is.na(coef)) {
    stop("coef must be a single name or numeric index.")
  }
  if (coef < 1 || coef > length(coef_names)) stop("coef index is out of range.")
  coef_names[coef]
}

get_strategy_table <- function(cellgraph, strategy_name, coef = NULL, model = c("full", "null")) {
  validate_motif_obj(cellgraph, require_offsets = FALSE)
  edger <- cellgraph$edger
  if (is.null(edger) || !is.list(edger)) {
    stop("cellgraph$edger is missing; run motif_edger() first.")
  }
  if (is.null(edger$strategies) || !is.list(edger$strategies)) {
    stop("cellgraph$edger$strategies is missing; rerun motif_edger().")
  }
  if (!strategy_name %in% names(edger$strategies)) {
    stop("Strategy not found: ", strategy_name, ". Available: ", paste(names(edger$strategies), collapse = ", "))
  }
  model <- match.arg(model)

  motif_info <- edger$motif_info
  if (is.null(motif_info) || !is.data.frame(motif_info)) {
    motif_info <- data.frame(motif = character(), motif_type = character(), stringsAsFactors = FALSE)
  }
  if (!all(c("motif", "motif_type") %in% names(motif_info))) {
    motif_info <- data.frame(motif = motif_info$motif, motif_type = NA_character_, stringsAsFactors = FALSE)
  }

  strat <- edger$strategies[[strategy_name]]
  if (is.null(strat) || !is.list(strat)) {
    stop("Strategy not found: ", strategy_name)
  }
  if (identical(strat$type, "edgeR")) {
    model_res <- strat[[model]]
    if (is.null(model_res) || !is.list(model_res)) {
      stop("Strategy ", strategy_name, " has no ", model, " model results.")
    }
    tests <- model_res$tests
    design_mat <- strat$design
    if (is.null(design_mat) || !ncol(design_mat)) stop("Design matrix missing for strategy: ", strategy_name)
    coef_name <- resolve_coef_name(design_mat, coef, colnames(design_mat))
    if (is.null(tests) || !length(tests) || is.null(tests[[coef_name]])) {
      warning("No edgeR tests stored for ", strategy_name, " coef: ", coef_name, ". Returning NA results.")
      out <- data.frame(
        motif = motif_info$motif,
        logFC = rep(NA_real_, nrow(motif_info)),
        PValue = rep(NA_real_, nrow(motif_info)),
        stringsAsFactors = FALSE
      )
    } else {
      tt <- edgeR::topTags(tests[[coef_name]], n = Inf, sort.by = "none")$table
      out <- data.frame(
        motif = rownames(tt),
        logFC = tt$logFC,
        PValue = tt$PValue,
        stringsAsFactors = FALSE
      )
    }
  } else if (identical(strat$type, "submotif_adj")) {
    coef_names <- strat$coef_names
    design_mat <- strat$design
    coef_name <- resolve_coef_name(design_mat, coef, coef_names)
    if (is.null(strat$logFC) || is.null(strat$PValue)) {
      warning("Submotif-adjusted results missing for strategy: ", strategy_name, ". Returning NA results.")
      out <- data.frame(
        motif = motif_info$motif,
        logFC = rep(NA_real_, nrow(motif_info)),
        PValue = rep(NA_real_, nrow(motif_info)),
        stringsAsFactors = FALSE
      )
    } else if (!coef_name %in% colnames(strat$logFC)) {
      stop("coef name not found in submotif-adjusted results.")
    } else {
      out <- data.frame(
        motif = rownames(strat$logFC),
        logFC = strat$logFC[, coef_name],
        PValue = strat$PValue[, coef_name],
        stringsAsFactors = FALSE
      )
    }
  } else {
    stop("Unknown strategy type for: ", strategy_name)
  }
  out$motif_type <- motif_info$motif_type[match(out$motif, motif_info$motif)]
  missing_type <- is.na(out$motif_type)
  if (any(missing_type)) {
    out$motif_type[missing_type] <- infer_motif_type(out$motif[missing_type])
  }
  out
}

#' Top motif results
#'
#' Convenience helpers for retrieving ranked motif results from [motif_edger()].
#' [top_edges()] returns edge motifs only. [top_motifs2()] returns wedge and
#' triangle motifs, with the corresponding edge-submotif statistics added as
#' columns. [motif_results()] returns both tables as a named list.
#'
#' @param cellgraph A `cellEdgeR_obj` with `edger` results.
#' @param coef Coefficient name or index; defaults to the first non-intercept coefficient.
#' @param model Which stored model to use: `full` or `null`.
#' @param n Number of motifs to return; defaults to all.
#' @param fdr_method Multiple testing correction method for `p.adjust` (default `BH`).
#' @param strategy Deprecated compatibility argument for [top_triplets()]; ignored.
#' @return A data frame with columns `motif`, `motif_type`, `logFC`,
#'   `PValue`, `FDR`, and `model_used`. `top_motifs2()` also includes
#'   `edge12`, `edge13`, `edge23` and their edge-level `logFC`, `PValue`,
#'   and `FDR` values.
#' @name top_edges
#' @aliases top_motifs2 top_triplets motif_results
#' @export
top_edges <- function(
  cellgraph,
  coef = NULL,
  model = c("full", "null"),
  n = Inf,
  fdr_method = "BH"
) {
  model <- match.arg(model)
  tbl <- get_strategy_table(cellgraph, strategy_name = "volume", coef = coef, model = model)
  tbl <- tbl[tbl$motif_type %in% "edge", , drop = FALSE]
  if (!nrow(tbl)) {
    return(data.frame(
      motif = character(),
      motif_type = character(),
      logFC = numeric(),
      PValue = numeric(),
      FDR = numeric(),
      model_used = character(),
      stringsAsFactors = FALSE
    ))
  }
  tbl$FDR <- stats::p.adjust(tbl$PValue, method = fdr_method)
  tbl$model_used <- "volume"
  tbl <- tbl[, c("motif", "motif_type", "logFC", "PValue", "FDR", "model_used")]
  tbl <- tbl[order(tbl$PValue, na.last = TRUE), , drop = FALSE]
  if (!is.infinite(n)) {
    n <- as.integer(n[1])
    if (is.na(n) || n < 0) stop("n must be non-negative.")
    tbl <- utils::head(tbl, n)
  }
  tbl
}

#' @rdname top_edges
#' @export
top_motifs2 <- function(
  cellgraph,
  coef = NULL,
  model = c("full", "null"),
  n = Inf,
  fdr_method = "BH"
) {
  model <- match.arg(model)
  tbl <- get_strategy_table(cellgraph, strategy_name = "volume", coef = coef, model = model)
  tbl <- tbl[tbl$motif_type %in% c("triangle", "wedge"), , drop = FALSE]
  if (!nrow(tbl)) {
    return(data.frame(
      motif = character(),
      motif_type = character(),
      logFC = numeric(),
      PValue = numeric(),
      FDR = numeric(),
      model_used = character(),
      edge12 = character(),
      edge13 = character(),
      edge23 = character(),
      edge12_logFC = numeric(),
      edge13_logFC = numeric(),
      edge23_logFC = numeric(),
      edge12_PValue = numeric(),
      edge13_PValue = numeric(),
      edge23_PValue = numeric(),
      edge12_FDR = numeric(),
      edge13_FDR = numeric(),
      edge23_FDR = numeric(),
      submotif_min_FDR = numeric(),
      submotif_max_abs_logFC = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  tbl$FDR <- stats::p.adjust(tbl$PValue, method = fdr_method)
  tbl$model_used <- "volume"
  edge_tbl <- top_edges(cellgraph, coef = coef, model = model, n = Inf, fdr_method = fdr_method)
  tbl <- add_motif2_submotif_stats(tbl, edge_tbl)
  tbl <- tbl[order(tbl$PValue, na.last = TRUE), , drop = FALSE]
  if (!is.infinite(n)) {
    n <- as.integer(n[1])
    if (is.na(n) || n < 0) stop("n must be non-negative.")
    tbl <- utils::head(tbl, n)
  }
  tbl
}

motif2_edge_keys <- function(motif, motif_type) {
  labels <- strsplit(strip_prefix(motif), "_", fixed = TRUE)[[1]]
  if (!motif_type %in% c("triangle", "wedge") || length(labels) != 3L) {
    return(c(edge12 = NA_character_, edge13 = NA_character_, edge23 = NA_character_))
  }
  c(
    edge12 = pair_key_vec(labels[1], labels[2], prefix = "E"),
    edge13 = pair_key_vec(labels[1], labels[3], prefix = "E"),
    edge23 = pair_key_vec(labels[2], labels[3], prefix = "E")
  )
}

add_motif2_submotif_stats <- function(tbl, edge_tbl) {
  edge_keys <- t(vapply(seq_len(nrow(tbl)), function(i) {
    motif2_edge_keys(tbl$motif[i], tbl$motif_type[i])
  }, character(3)))
  edge_keys <- as.data.frame(edge_keys, stringsAsFactors = FALSE)
  names(edge_keys) <- c("edge12", "edge13", "edge23")
  out <- cbind(tbl[, c("motif", "motif_type", "logFC", "PValue", "FDR", "model_used"), drop = FALSE], edge_keys)

  add_edge_stat <- function(stat_col, out_prefix) {
    for (edge_col in c("edge12", "edge13", "edge23")) {
      values <- rep(NA_real_, nrow(out))
      if (!is.null(edge_tbl) && nrow(edge_tbl) && stat_col %in% names(edge_tbl)) {
        idx <- match(out[[edge_col]], edge_tbl$motif)
        ok <- !is.na(idx)
        values[ok] <- as.numeric(edge_tbl[[stat_col]][idx[ok]])
      }
      out[[paste0(edge_col, "_", out_prefix)]] <<- values
    }
  }
  add_edge_stat("logFC", "logFC")
  add_edge_stat("PValue", "PValue")
  add_edge_stat("FDR", "FDR")

  fdr_mat <- as.matrix(out[, c("edge12_FDR", "edge13_FDR", "edge23_FDR"), drop = FALSE])
  logfc_mat <- as.matrix(out[, c("edge12_logFC", "edge13_logFC", "edge23_logFC"), drop = FALSE])
  out$submotif_min_FDR <- apply(fdr_mat, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  out$submotif_max_abs_logFC <- apply(logfc_mat, 1, function(x) if (all(is.na(x))) NA_real_ else max(abs(x), na.rm = TRUE))
  out
}

#' @rdname top_edges
#' @export
top_triplets <- function(
  cellgraph,
  coef = NULL,
  model = c("full", "null"),
  n = Inf,
  fdr_method = "BH",
  strategy = NULL
) {
  if (!is.null(strategy)) {
    warning("top_triplets(strategy=...) is deprecated; returning top_motifs2() volume results and ignoring strategy.")
  }
  top_motifs2(cellgraph = cellgraph, coef = coef, model = model, n = n, fdr_method = fdr_method)
}

#' @rdname top_edges
#' @export
motif_results <- function(
  cellgraph,
  coef = NULL,
  model = c("full", "null"),
  n = Inf,
  fdr_method = "BH"
) {
  model <- match.arg(model)
  list(
    edges = top_edges(cellgraph, coef = coef, model = model, n = n, fdr_method = fdr_method),
    motifs2 = top_motifs2(cellgraph, coef = coef, model = model, n = n, fdr_method = fdr_method)
  )
}

#' Count the number of possible tested motif labels
#'
#' Compute the combinatorial number of possible edge, wedge, and triangle label
#' motifs implied by a `cellEdgeR_obj`. Counts reflect label combinations, not
#' graph isomorphism classes. Node motifs and merged triplets are not part of the
#' tested motif space.
#'
#' @param cellgraph A `cellEdgeR_obj` (typically after [motif_edger()]).
#' @param include_wedge Logical; whether to include wedge motifs in the count.
#'   Defaults to `cellgraph$parameters$include_wedge` when available.
#' @return A list with `labels`, `include_wedge`, `counts` (data frame), and
#'   `total` (sum of possible tested motifs).
#' @export
motif_space_size <- function(cellgraph, include_wedge = NULL) {
  if (!inherits(cellgraph, "cellEdgeR_obj")) {
    stop("cellgraph must be a cellEdgeR_obj.")
  }
  K <- length(cellgraph$label_levels)
  if (K < 1) stop("cellgraph has no label levels.")
  if (is.null(include_wedge)) {
    include_wedge <- isTRUE(cellgraph$parameters$include_wedge)
  }
  edge_n <- K * (K + 1L) / 2L
  tri_n <- choose(K + 2L, 3L)
  wedge_n <- if (include_wedge) K * (K * (K + 1L) / 2L) else 0L
  counts <- data.frame(
    layer = c("edge", "triangle", if (include_wedge) "wedge"),
    n_possible = c(edge_n, tri_n, if (include_wedge) wedge_n),
    stringsAsFactors = FALSE
  )
  list(
    labels = K,
    include_wedge = include_wedge,
    counts = counts,
    total = sum(counts$n_possible)
  )
}
