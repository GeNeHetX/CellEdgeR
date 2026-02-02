#' Boxplot of normalized motif counts (internal helper)
#'
#' Visualize normalized motif counts (e.g., from [get_norm_counts()]) across samples for a given motif.
#'
#' @param norm_counts List of normalized counts (e.g., `cellgraph$norm_counts` from [count_motifs_graphs()]); when multiple offset modes are present, select one via `offset_mode`.
#' @param motif_key Character motif identifier with prefix (e.g., `"E_A_B"`).
#' @param layer Which layer to plot; when `NULL`, inferred from the motif prefix (`N_`, `E_`, `T_`, or `W_`).
#' @param sample_df Optional data frame of sample metadata; rownames must match the motif columns.
#' @param group_var Optional column name in `sample_df` to use for grouping/coloring the boxplot.
#' @return A `ggplot` object.
#' @keywords internal
plot_motif_box <- function(norm_counts, motif_key, layer = NULL,
                           sample_df = NULL, group_var = NULL,
                           offset_mode = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  if (!is.list(norm_counts)) {
    stop("norm_counts must be the list returned by count_motifs_graphs()$norm_counts or similar.")
  }
  # Pick offset mode when multiple are present
  if (!all(c("node", "edge", "triangle") %in% names(norm_counts))) {
    avail <- names(norm_counts)
    if (is.null(offset_mode)) offset_mode <- avail[1]
    if (!offset_mode %in% avail) stop("offset_mode not found in norm_counts; available: ", paste(avail, collapse = ", "))
    norm_counts <- norm_counts[[offset_mode]]
  }
  infer_layer <- function(key) {
    prefix <- sub("_.*", "", key)
    switch(prefix,
      N = "node",
      E = "edge",
      T = "triangle",
      W = "wedge",
      stop("Cannot infer motif layer from key: ", key)
    )
  }
  if (is.null(layer)) layer <- infer_layer(motif_key)
  if (length(layer) != 1 || !layer %in% names(norm_counts)) stop("Layer ", layer, " not found in norm_counts.")
  mat <- norm_counts[[layer]]
  if (!motif_key %in% rownames(mat)) stop("Motif key not found in layer ", layer, ".")
  vals <- as.numeric(mat[motif_key, ])
  df <- data.frame(
    sample = colnames(mat),
    value = vals,
    stringsAsFactors = FALSE
  )
  if (!is.null(sample_df)) {
    if (is.null(rownames(sample_df))) stop("sample_df must have rownames matching samples.")
    sample_df <- sample_df[df$sample, , drop = FALSE]
    df <- cbind(df, sample_df)
    if (!is.null(group_var)) {
      if (!group_var %in% names(sample_df)) stop("group_var not found in sample_df.")
    }
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = if (!is.null(group_var)) .data[[group_var]] else "all",
                                        y = value,
                                        fill = if (!is.null(group_var)) .data[[group_var]] else NULL)) +
    ggplot2::geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = 21, outlier.fill = "white") +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.6, size = 2, color = "#1b4f72") +
    ggplot2::ylab(paste0("Normalized count: ", motif_key)) +
    ggplot2::xlab(if (!is.null(group_var)) group_var else "") +
    ggplot2::ggtitle(paste0("Normalized ", layer, " motif: ", motif_key)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  p
}

resolve_layer_container <- function(container, offset_mode, layer, label) {
  if (is.null(container) || !is.list(container) || !length(container)) {
    stop(label, " is missing; run count_motifs_graphs() first.")
  }
  if (all(c("node", "edge", "triangle") %in% names(container))) {
    if (is.null(offset_mode)) offset_mode <- "volume"
    if (!identical(offset_mode, "volume")) {
      stop("offset_mode not found in ", label, "; available: volume")
    }
    layer_list <- container
  } else {
    if (is.null(offset_mode)) offset_mode <- names(container)[1]
    if (!offset_mode %in% names(container)) {
      stop("offset_mode not found in ", label, "; available: ", paste(names(container), collapse = ", "))
    }
    layer_list <- container[[offset_mode]]
  }
  if (is.null(layer_list) || !is.list(layer_list) || is.null(layer_list[[layer]])) {
    stop(label, " does not include layer: ", layer)
  }
  list(layer = layer_list[[layer]], offset_mode = offset_mode)
}

align_layer_matrix <- function(mat, ref_rows, ref_cols, label) {
  if (is.null(mat)) stop(label, " is missing.")
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(label, " must have row and column names.")
  }
  if (!setequal(colnames(mat), ref_cols)) {
    stop(label, " columns do not match sample names.")
  }
  if (!identical(colnames(mat), ref_cols)) {
    mat <- mat[, ref_cols, drop = FALSE]
  }
  if (!setequal(rownames(mat), ref_rows)) {
    stop(label, " row names do not match motif keys.")
  }
  if (!identical(rownames(mat), ref_rows)) {
    mat <- mat[ref_rows, , drop = FALSE]
  }
  mat
}

#' Slope test for geometric scaling
#'
#' Compare observed triangle counts to the structural expectation from sub-motif offsets
#' using a log-log regression.
#'
#' @param cellgraph Output of [count_motifs_graphs()] (or [motif_edger()]) with stored offsets.
#' @param offset_mode Offset set to use for the expectation; defaults to \code{"hier_null"}.
#' @param layer Motif layer to plot; defaults to \code{"triangle"}.
#' @param motif_key Optional character vector of motif keys to include; defaults to all motifs.
#' @param log_base Logarithm base for the regression; defaults to natural log.
#' @param pseudocount Positive value added before log transform; defaults to \code{cellgraph$parameters$offset_pseudo} or 1.
#' @return A \code{ggplot} object with the fitted slope annotation.
#' @export
plot_motif_slope_test <- function(cellgraph,
                                  offset_mode = "hier_null",
                                  layer = "triangle",
                                  motif_key = NULL,
                                  log_base = exp(1),
                                  pseudocount = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  validate_motif_obj(cellgraph, require_offsets = TRUE)
  counts <- cellgraph$raw_count[[layer]]
  if (is.null(counts) || nrow(counts) == 0) stop("raw_count layer is missing or empty: ", layer)
  offsets <- resolve_layer_container(cellgraph$offsets, offset_mode, layer, "offsets")
  offs_layer <- align_layer_matrix(offsets$layer, rownames(counts), colnames(counts), "offsets layer")
  if (!is.null(motif_key)) {
    motif_key <- as.character(motif_key)
    keep <- rownames(counts) %in% motif_key
    if (!any(keep)) stop("motif_key not found in ", layer, " layer.")
    counts <- counts[keep, , drop = FALSE]
    offs_layer <- offs_layer[rownames(counts), , drop = FALSE]
  }
  pseudocount <- if (is.null(pseudocount)) {
    pc <- cellgraph$parameters$offset_pseudo
    if (is.numeric(pc) && length(pc) == 1 && is.finite(pc) && pc > 0) pc else 1
  } else {
    as.numeric(pseudocount)[1]
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) stop("pseudocount must be a positive number.")
  if (!is.finite(log_base) || log_base <= 0 || log_base == 1) stop("log_base must be a positive number not equal to 1.")

  obs <- as.matrix(counts)
  log_obs <- log(as.numeric(obs) + pseudocount, base = log_base)
  log_exp <- as.numeric(offs_layer) / log(log_base)
  df <- data.frame(
    motif = rep(rownames(obs), times = ncol(obs)),
    sample_name = rep(colnames(obs), each = nrow(obs)),
    log_obs = log_obs,
    log_exp = log_exp,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$log_obs) & is.finite(df$log_exp), , drop = FALSE]
  if (nrow(df) < 2) stop("Not enough finite points to fit the slope test.")

  fit <- stats::lm(log_obs ~ log_exp, data = df)
  slope <- unname(stats::coef(fit)[2])
  r2 <- summary(fit)$r.squared
  label_txt <- sprintf("slope = %.3f\\nR^2 = %.3f", slope, r2)
  title_txt <- paste0("Slope test (", offsets$offset_mode, " offsets, ", layer, ")")
  xlab_txt <- paste0("log", if (isTRUE(all.equal(log_base, exp(1)))) "" else paste0("[base ", base::sprintf("%g", log_base), "]"),
    " expected counts")
  ylab_txt <- paste0("log", if (isTRUE(all.equal(log_base, exp(1)))) "" else paste0("[base ", base::sprintf("%g", log_base), "]"),
    " observed counts")

  ggplot2::ggplot(df, ggplot2::aes(x = log_exp, y = log_obs)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.3, color = "#1b4f72") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#e67e22", linewidth = 0.8) +
    ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.2,
      label = label_txt, size = 3.2) +
    ggplot2::labs(
      title = title_txt,
      x = xlab_txt,
      y = ylab_txt
    ) +
    ggplot2::theme_minimal()
}

#' Artifact check: residuals versus edge density
#'
#' Plot per-sample log residuals from edge-derived offsets against total edge density.
#'
#' @param cellgraph Output of [count_motifs_graphs()] (or [motif_edger()]) with stored relative counts.
#' @param offset_mode Offset set to use for residuals; defaults to \code{"hier_null"}.
#' @param layer Motif layer to plot; defaults to \code{"triangle"}.
#' @param motif_key Optional character vector of motif keys to include; defaults to all motifs.
#' @param log_base Logarithm base for the x-axis; defaults to 2 to match residuals.
#' @param pseudocount Positive value added before log transform; defaults to \code{cellgraph$parameters$offset_pseudo} or 1.
#' @param smooth Logical; add a linear trend line.
#' @return A \code{ggplot} object.
#' @export
plot_motif_artifact_check <- function(cellgraph,
                                      offset_mode = "hier_null",
                                      layer = "triangle",
                                      motif_key = NULL,
                                      log_base = 2,
                                      pseudocount = NULL,
                                      smooth = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  validate_motif_obj(cellgraph, require_offsets = TRUE)
  rel_counts <- resolve_layer_container(cellgraph$relative_counts, offset_mode, layer, "relative_counts")
  rel_layer <- rel_counts$layer
  samples <- cellgraph$sample_name
  rel_layer <- align_layer_matrix(rel_layer, rownames(rel_layer), samples, "relative_counts layer")
  if (!nrow(rel_layer) || !ncol(rel_layer)) {
    stop("relative_counts layer is empty for ", layer, ". ",
      "Try a different layer or rerun count_motifs_graphs() with settings that yield motifs.")
  }
  if (!is.null(motif_key)) {
    motif_key <- as.character(motif_key)
    keep <- rownames(rel_layer) %in% motif_key
    if (!any(keep)) stop("motif_key not found in ", layer, " layer.")
    rel_layer <- rel_layer[keep, , drop = FALSE]
  }
  if (!nrow(rel_layer)) {
    stop("No motifs available for layer ", layer, " after applying motif_key.")
  }
  pseudocount <- if (is.null(pseudocount)) {
    pc <- cellgraph$parameters$offset_pseudo
    if (is.numeric(pc) && length(pc) == 1 && is.finite(pc) && pc > 0) pc else 1
  } else {
    as.numeric(pseudocount)[1]
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) stop("pseudocount must be a positive number.")
  if (!is.finite(log_base) || log_base <= 0 || log_base == 1) stop("log_base must be a positive number not equal to 1.")

  edges <- cellgraph$exposure$edges
  if (is.null(edges) || !length(edges)) stop("cellgraph$exposure$edges is missing.")
  if (is.null(names(edges))) {
    if (length(edges) != length(samples)) {
      stop("cellgraph$exposure$edges must be named by sample names.")
    }
    names(edges) <- samples
  } else if (!setequal(names(edges), samples)) {
    stop("cellgraph$exposure$edges names do not match sample names.")
  }
  edges <- edges[samples]
  edge_log <- stats::setNames(log(as.numeric(edges) + pseudocount, base = log_base), names(edges))

  res_mat <- as.matrix(rel_layer)
  res_vec <- as.numeric(res_mat)
  if (!isTRUE(all.equal(log_base, 2))) {
    res_vec <- res_vec / log2(log_base)
  }
  df <- data.frame(
    motif = rep(rownames(res_mat), times = ncol(res_mat)),
    sample_name = rep(colnames(res_mat), each = nrow(res_mat)),
    log_edges = rep(edge_log[colnames(res_mat)], each = nrow(res_mat)),
    residual = res_vec,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$log_edges) & is.finite(df$residual), , drop = FALSE]
  if (!nrow(df)) {
    stop("No finite values available for the artifact check. ",
      "Check that residuals and edge totals are finite for the selected layer.")
  }

  title_txt <- paste0("Artifact check (", rel_counts$offset_mode, " residuals, ", layer, ")")
  xlab_txt <- paste0("log", if (isTRUE(all.equal(log_base, exp(1)))) "" else paste0("[base ", base::sprintf("%g", log_base), "]"),
    " total edges")
  ylab_txt <- if (isTRUE(all.equal(log_base, 2))) {
    "Residual (log2)"
  } else {
    paste0("Residual (log base ", base::sprintf("%g", log_base), ")")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log_edges, y = residual)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.3, color = "#1b4f72") +
    ggplot2::labs(
      title = title_txt,
      x = xlab_txt,
      y = ylab_txt
    ) +
    ggplot2::theme_minimal()
  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#e67e22", linewidth = 0.8)
  }
  p
}

#' Plot a sample graph with highlighted motifs
#'
#' Draw the Delaunay edges and cell coordinates for a single sample, optionally highlighting
#' a motif by labels (and edges connecting those labels).
#'
#' @param graph_obj Output of [build_cell_graphs()] with stored coordinates.
#' @param sample_name Sample name to plot.
#' @param max_edge_len Optional numeric threshold to prune long edges for display; set `Inf` to keep all.
#' @param highlight_labels Optional character vector of labels to emphasize.
#' @param motif_key Optional motif identifier (e.g., `"E_A_B"`); when provided, the labels in the key
#'   are highlighted and edges connecting those labels are accentuated.
#' @param motif_layer Layer for the motif key; when `NULL`, inferred from the motif prefix.
#' @param cells_by_sample Optional named list of raw sample data frames; used only when \code{graph_obj}
#'   lacks stored coordinates (backward compatibility with older objects).
#' @param motif_node_size Size for nodes participating in the highlighted motif.
#' @param dim_node_nonmotif Factor to shrink nodes that are not in the highlighted motif.
#' @param alpha_node_nonmotif Alpha for nodes that are not in the highlighted motif.
#' @param alpha_edge_nonmotif Alpha for edges that are not in the highlighted motif.
#' @return A `ggplot` object.
#' @export
plot_sample_graph <- function(graph_obj, sample_name, max_edge_len = Inf, highlight_labels = NULL,
                              motif_key = NULL, motif_layer = NULL,
                              cells_by_sample = NULL,
                              motif_node_size = 3,
                              dim_node_nonmotif = 0.6,
                              alpha_node_nonmotif = 0.4,
                              alpha_edge_nonmotif = 0.3) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  validate_graph_obj(graph_obj)
  if (!sample_name %in% graph_obj$sample_name) stop("sample_name not found in graph_obj.")
  ps <- graph_obj$per_sample_graph[[sample_name]]
  if (is.null(ps)) stop("graph_obj has no per-sample graph for sample_name: ", sample_name)
  if (is.null(ps$xy)) {
    if (!is.null(cells_by_sample) && sample_name %in% names(cells_by_sample)) {
      df <- standardize_sample_df(cells_by_sample[[sample_name]], sample_name)
      xy_tmp <- as.matrix(df[, c("x", "y")])
      labs_chr_tmp <- as.character(df$label)
      if (!identical(labs_chr_tmp, ps$labels_chr)) {
        stop("Labels from cells_by_sample do not match graph_obj labels; cannot plot.")
      }
      ps$xy <- xy_tmp
    } else {
      stop("graph_obj does not store coordinates; rebuild graphs with the current version of build_cell_graphs() or provide cells_by_sample.")
    }
  }
  xy <- ps$xy
  motif_labels <- NULL
  motif_pairs <- NULL
  nodes_in_motif <- rep(FALSE, nrow(ps$xy))
  if (!is.null(motif_key)) {
    infer_layer <- function(key) {
      prefix <- sub("_.*", "", key)
      switch(prefix,
        N = "node",
        E = "edge",
        T = "triangle",
      W = "wedge",
        TW = "triangle",
        stop("Cannot infer motif layer from key: ", key)
      )
    }
    if (is.null(motif_layer)) motif_layer <- infer_layer(motif_key)
    motif_labels <- strsplit(sub("^[^_]+_", "", motif_key), "_")[[1]]
    if (motif_layer == "edge") {
      motif_pairs <- paste(sort(motif_labels[seq_len(min(2, length(motif_labels)))]), collapse = "_")
    }
    highlight_labels <- unique(c(highlight_labels, motif_labels))
  }
  if (is.null(motif_layer)) motif_layer <- "triangle"

  nodes <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    label = ps$labels_chr,
    highlight = if (is.null(highlight_labels)) FALSE else ps$labels_chr %in% highlight_labels,
    stringsAsFactors = FALSE
  )
  edges <- ps$edges
  if (length(edges)) {
    cutoff <- if (is.finite(max_edge_len)) {
      max_edge_len
    } else if (!is.null(graph_obj$parameters$max_edge_len) && is.finite(graph_obj$parameters$max_edge_len)) {
      graph_obj$parameters$max_edge_len
    } else {
      Inf
    }
    if (is.finite(cutoff)) {
      keep <- ps$edge_len <= cutoff
      edges <- edges[keep, , drop = FALSE]
    }
  }

  highlight_edges <- logical(nrow(edges))
  if (length(edges) && !is.null(motif_labels)) {
    la <- ps$labels_chr[edges[, 1]]
    lb <- ps$labels_chr[edges[, 2]]
    if (motif_layer == "edge" && !is.null(motif_pairs)) {
      pairs <- ifelse(la <= lb, paste(la, lb, sep = "_"), paste(lb, la, sep = "_"))
      highlight_edges <- pairs == motif_pairs
    } else {
      highlight_edges <- la %in% motif_labels & lb %in% motif_labels
    }
    if (any(highlight_edges)) {
      nodes_in_motif[unique(c(edges[highlight_edges, 1], edges[highlight_edges, 2]))] <- TRUE
    }
  } else if (!is.null(motif_labels) && motif_layer == "node") {
    nodes_in_motif <- ps$labels_chr %in% motif_labels
  }

  edge_df <- if (length(edges)) {
    data.frame(
      x1 = xy[edges[, 1], 1],
      y1 = xy[edges[, 1], 2],
      x2 = xy[edges[, 2], 1],
      y2 = xy[edges[, 2], 2]
    )
  } else {
    data.frame()
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
      color = "#9aa5b1",
      linewidth = 0.3,
      alpha = alpha_edge_nonmotif
    ) +
    ggplot2::geom_segment(
      data = edge_df[highlight_edges, , drop = FALSE],
      ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
      color = "#e67e22",
      linewidth = 0.6,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, color = label,
                   size = dplyr::if_else(nodes_in_motif, motif_node_size, dim_node_nonmotif),
                   alpha = dplyr::if_else(nodes_in_motif, 1, alpha_node_nonmotif))
    ) +
    ggplot2::scale_size_identity(guide = "none") +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::ggtitle(paste0("Sample ", sample_name, " graph")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  p
}
