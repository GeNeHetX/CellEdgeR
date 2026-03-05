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
        TP = "triplet",
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
  if (!is.null(motif_key)) {
    key_prefix <- sub("_.*", "", motif_key)
    key_type <- switch(key_prefix,
      N = "node",
      E = "edge",
      T = "triangle",
      W = "wedge",
      TP = "triplet",
      TW = "triangle",
      if (!is.null(motif_layer)) motif_layer else "triangle"
    )
  } else {
    key_type <- if (!is.null(motif_layer)) motif_layer else "triangle"
  }

  if (length(edges) && !is.null(motif_labels) && key_type %in% c("triangle", "wedge", "triplet")) {
    target <- motif_labels
    edge_key <- paste(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]), sep = "_")
    edge_idx <- seq_len(nrow(edges))
    names(edge_idx) <- edge_key
    adj <- vector("list", ps$n)
    for (i in seq_len(nrow(edges))) {
      a <- edges[i, 1]
      b <- edges[i, 2]
      adj[[a]] <- c(adj[[a]], b)
      adj[[b]] <- c(adj[[b]], a)
    }
    match_unordered <- function(x, y) {
      length(x) == length(y) && all(sort(x) == sort(y))
    }
    match_wedge <- function(center_lab, leaf_labs, target_labs) {
      if (length(target_labs) != 3) return(FALSE)
      if (!identical(center_lab, target_labs[1])) return(FALSE)
      all(sort(leaf_labs) == sort(target_labs[2:3]))
    }

    want_triangle <- key_type %in% c("triangle", "triplet")
    want_wedge <- key_type %in% c("wedge", "triplet")

    if (want_triangle) {
      for (i in seq_len(ps$n)) {
        ni <- adj[[i]]
        if (length(ni) < 2) next
        ni <- ni[ni > i]
        if (length(ni) < 2) next
        for (j in ni) {
          nj <- adj[[j]]
          if (!length(nj)) next
          common <- intersect(ni[ni > j], nj)
          if (!length(common)) next
          for (k in common) {
            labs <- c(ps$labels_chr[i], ps$labels_chr[j], ps$labels_chr[k])
            if (match_unordered(labs, target)) {
              nodes_in_motif[c(i, j, k)] <- TRUE
              e1 <- edge_idx[paste(pmin(i, j), pmax(i, j), sep = "_")]
              e2 <- edge_idx[paste(pmin(i, k), pmax(i, k), sep = "_")]
              e3 <- edge_idx[paste(pmin(j, k), pmax(j, k), sep = "_")]
              idx <- c(e1, e2, e3)
              idx <- idx[!is.na(idx)]
              if (length(idx)) highlight_edges[idx] <- TRUE
            }
          }
        }
      }
    }

    if (want_wedge) {
      for (j in seq_len(ps$n)) {
        nj <- adj[[j]]
        if (length(nj) < 2) next
        comb <- utils::combn(nj, 2)
        for (cidx in seq_len(ncol(comb))) {
          i <- comb[1, cidx]
          k <- comb[2, cidx]
          # open wedge: i and k not connected
          if (!is.na(edge_idx[paste(pmin(i, k), pmax(i, k), sep = "_")])) next
          labs <- c(ps$labels_chr[j], ps$labels_chr[i], ps$labels_chr[k])
          ok <- if (key_type == "wedge") {
            match_wedge(labs[1], labs[2:3], target)
          } else {
            match_unordered(labs, target)
          }
          if (ok) {
            nodes_in_motif[c(i, j, k)] <- TRUE
            e1 <- edge_idx[paste(pmin(i, j), pmax(i, j), sep = "_")]
            e2 <- edge_idx[paste(pmin(j, k), pmax(j, k), sep = "_")]
            idx <- c(e1, e2)
            idx <- idx[!is.na(idx)]
            if (length(idx)) highlight_edges[idx] <- TRUE
          }
        }
      }
    }
  } else if (length(edges) && !is.null(motif_labels)) {
    la <- ps$labels_chr[edges[, 1]]
    lb <- ps$labels_chr[edges[, 2]]
    if (key_type == "edge" && !is.null(motif_pairs)) {
      pairs <- ifelse(la <= lb, paste(la, lb, sep = "_"), paste(lb, la, sep = "_"))
      highlight_edges <- pairs == motif_pairs
    } else {
      highlight_edges <- la %in% motif_labels & lb %in% motif_labels
    }
    if (any(highlight_edges)) {
      nodes_in_motif[unique(c(edges[highlight_edges, 1], edges[highlight_edges, 2]))] <- TRUE
    }
  } else if (!is.null(motif_labels) && key_type == "node") {
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

resolve_plot_results_table <- function(x, result_type, coef, model, fdr_method) {
  result_type <- match.arg(result_type, c("edges", "triplets_cooccurrence", "triplets_topology"))
  model <- match.arg(model, c("full", "null"))

  if (is.data.frame(x)) {
    out <- x
  } else if (inherits(x, "cellEdgeR_obj")) {
    out <- switch(result_type,
      edges = top_edges(x, coef = coef, model = model, n = Inf, fdr_method = fdr_method),
      triplets_cooccurrence = top_triplets(x, strategy = "cooccurrence", coef = coef, model = model, n = Inf, fdr_method = fdr_method),
      triplets_topology = top_triplets(x, strategy = "topology", coef = coef, model = "full", n = Inf, fdr_method = fdr_method)
    )
  } else {
    stop("x must be either a results data.frame or a cellEdgeR_obj.")
  }

  if (!is.data.frame(out) || !nrow(out)) {
    stop("No motif results available to plot.")
  }
  if (!"motif" %in% names(out)) {
    rn <- rownames(out)
    if (is.null(rn) || !length(rn)) stop("Results table must include a 'motif' column.")
    out$motif <- rn
  }
  need_cols <- c("logFC", "PValue")
  missing_cols <- setdiff(need_cols, names(out))
  if (length(missing_cols)) {
    stop("Results table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  out$motif <- as.character(out$motif)
  out$logFC <- as.numeric(out$logFC)
  out$PValue <- as.numeric(out$PValue)
  if (!"FDR" %in% names(out) || all(!is.finite(as.numeric(out$FDR)))) {
    out$FDR <- stats::p.adjust(out$PValue, method = fdr_method)
  } else {
    out$FDR <- as.numeric(out$FDR)
    miss_fdr <- !is.finite(out$FDR)
    if (any(miss_fdr)) {
      out$FDR <- stats::p.adjust(out$PValue, method = fdr_method)
    }
  }
  if (!"motif_type" %in% names(out)) {
    out$motif_type <- infer_motif_type(out$motif)
  } else {
    out$motif_type <- as.character(out$motif_type)
    miss_type <- is.na(out$motif_type) | out$motif_type == ""
    if (any(miss_type)) out$motif_type[miss_type] <- infer_motif_type(out$motif[miss_type])
  }
  out
}

prepare_volcano_dataframe <- function(tbl, fdr_cutoff, logFC_cutoff) {
  out <- tbl
  finite_p <- out$PValue[is.finite(out$PValue) & out$PValue > 0]
  min_pos <- if (length(finite_p)) min(finite_p) else 1e-300
  p_plot <- out$PValue
  p_plot[!is.finite(p_plot) | p_plot <= 0] <- min_pos
  out$neg_log10_p <- -log10(p_plot)

  sig <- is.finite(out$FDR) & out$FDR <= fdr_cutoff
  cutoff <- abs(as.numeric(logFC_cutoff)[1])
  dep <- sig & is.finite(out$logFC) & out$logFC < -cutoff
  enr <- sig & is.finite(out$logFC) & out$logFC > cutoff
  cls <- rep("Not significant", nrow(out))
  cls[dep] <- "Depleted"
  cls[enr] <- "Enriched"
  out$volcano_group <- factor(cls, levels = c("Depleted", "Not significant", "Enriched"))
  out
}

select_volcano_labels <- function(df, label_mode, label_motifs, label_n, label_top_by) {
  label_mode <- match.arg(label_mode, c("none", "manual", "top_n", "top_bottom_n"))
  label_top_by <- match.arg(label_top_by, c("PValue", "FDR", "abs_logFC"))
  if (identical(label_mode, "none")) return(character())

  if (identical(label_mode, "manual")) {
    if (is.null(label_motifs) || !length(label_motifs)) {
      stop("label_mode = \"manual\" requires label_motifs.")
    }
    label_motifs <- unique(as.character(label_motifs))
    missing <- setdiff(label_motifs, df$motif)
    if (length(missing)) {
      warning("Ignoring motifs not found in results: ", paste(missing, collapse = ", "))
    }
    return(intersect(label_motifs, df$motif))
  }

  label_n <- as.integer(label_n[1])
  if (!is.finite(label_n) || is.na(label_n) || label_n < 1) {
    stop("label_n must be >= 1 for label_mode = \"", label_mode, "\".")
  }

  keep <- is.finite(df$logFC) & is.finite(df$PValue)
  d <- df[keep, , drop = FALSE]
  if (!nrow(d)) return(character())

  if (identical(label_mode, "top_n")) {
    ord <- switch(label_top_by,
      PValue = order(d$PValue, -abs(d$logFC), na.last = NA),
      FDR = order(d$FDR, d$PValue, na.last = NA),
      abs_logFC = order(-abs(d$logFC), d$PValue, na.last = NA)
    )
    return(unique(utils::head(d$motif[ord], label_n)))
  }

  dep <- d[order(d$logFC, d$PValue, na.last = NA), , drop = FALSE]
  enr <- d[order(-d$logFC, d$PValue, na.last = NA), , drop = FALSE]
  unique(c(utils::head(dep$motif, label_n), utils::head(enr$motif, label_n)))
}

build_motif_volcano_plot <- function(df, labels, label_map = NULL, fdr_cutoff, logFC_cutoff, point_size, point_alpha, title, subtitle) {
  d <- df[is.finite(df$logFC) & is.finite(df$neg_log10_p), , drop = FALSE]
  if (!nrow(d)) stop("No finite motif values to plot.")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = neg_log10_p)) +
    ggplot2::geom_point(
      ggplot2::aes(color = volcano_group),
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = c(
        Depleted = "#2b8cbe",
        "Not significant" = "#bdbdbd",
        Enriched = "#d7301f"
      ),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 fold change",
      y = "-log10(PValue)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )

  cutoff <- abs(as.numeric(logFC_cutoff)[1])
  if (is.finite(cutoff) && cutoff > 0) {
    p <- p +
      ggplot2::geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed", linewidth = 0.4, color = "#636363")
  }

  x_rng <- range(d$logFC, na.rm = TRUE, finite = TRUE)
  x_note <- if (length(x_rng) == 2 && is.finite(x_rng[1]) && is.finite(x_rng[2]) &&
    x_rng[1] <= 0 && x_rng[2] >= 0) {
    0
  } else {
    mean(x_rng)
  }
  p <- p + ggplot2::annotate("text",
    x = x_note, y = Inf, hjust = 0.5, vjust = 1.3, size = 3.1,
    label = paste0("FDR <= ", format(fdr_cutoff, digits = 2))
  )

  if (length(labels)) {
    label_df <- d[d$motif %in% labels, , drop = FALSE]
    if (nrow(label_df)) {
      if (!is.null(label_map) && length(label_map)) {
        if (is.null(names(label_map)) || any(names(label_map) == "")) {
          stop("label_map must be a named character vector when provided.")
        }
        label_txt <- as.character(label_map[label_df$motif])
        miss <- !is.finite(match(label_df$motif, names(label_map))) | is.na(label_txt) | label_txt == ""
        label_txt[miss] <- label_df$motif[miss]
        label_df$label_text <- label_txt
      } else {
        label_df$label_text <- label_df$motif
      }
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = label_text),
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(label = label_text),
          size = 3,
          vjust = -0.4,
          check_overlap = TRUE
        )
      }
    }
  }

  p
}

pick_side_motifs <- function(df, n_side, left_motifs = NULL, right_motifs = NULL) {
  n_side <- as.integer(n_side[1])
  if (!is.finite(n_side) || is.na(n_side) || n_side < 0) stop("n_side must be >= 0.")
  base <- df[is.finite(df$logFC) & is.finite(df$PValue), , drop = FALSE]

  pick_manual <- function(keys, side) {
    if (is.null(keys)) return(NULL)
    keys <- unique(as.character(keys))
    idx <- match(keys, base$motif)
    missing <- keys[is.na(idx)]
    if (length(missing)) warning("Ignoring ", side, " motifs not found: ", paste(missing, collapse = ", "))
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(base[FALSE, , drop = FALSE])
    base[idx, , drop = FALSE]
  }

  pick_auto <- function(side, exclude = character()) {
    if (n_side == 0 || !nrow(base)) return(base[FALSE, , drop = FALSE])
    pool <- base[!base$motif %in% exclude, , drop = FALSE]
    if (!nrow(pool)) return(pool)
    if (identical(side, "left")) {
      pool <- pool[order(pool$logFC, pool$PValue, na.last = NA), , drop = FALSE]
    } else {
      pool <- pool[order(-pool$logFC, pool$PValue, na.last = NA), , drop = FALSE]
    }
    utils::head(pool, n_side)
  }

  left <- pick_manual(left_motifs, "left")
  right <- pick_manual(right_motifs, "right")

  if (is.null(left)) {
    right_keys <- if (is.null(right)) character() else right$motif
    left <- pick_auto("left", exclude = right_keys)
  }
  if (is.null(right)) {
    right <- pick_auto("right", exclude = left$motif)
  }
  list(left = left, right = right)
}

truncate_motif_labels <- function(x, n_chars = 3) {
  x <- as.character(x)
  if (length(n_chars) != 1 || is.null(n_chars)) {
    stop("side_label_truncate must be a single integer >= 1 or NA for no truncation.")
  }
  if (is.na(n_chars)) return(x)
  n_chars <- as.integer(n_chars)
  if (!is.finite(n_chars) || n_chars < 1) {
    stop("side_label_truncate must be >= 1 or NA for no truncation.")
  }
  substr(x, 1L, n_chars)
}

motif_layout_template <- function(motif, motif_type, side_label_truncate = 3) {
  labs <- strsplit(strip_prefix(motif), "_", fixed = TRUE)[[1]]
  if (!length(labs)) labs <- motif
  labs <- truncate_motif_labels(labs, n_chars = side_label_truncate)
  mt <- motif_type
  if (is.na(mt) || mt == "") mt <- infer_motif_type(motif)
  if (is.na(mt) || mt == "") mt <- "triplet"

  if (identical(mt, "node")) {
    if (length(labs) < 1) labs <- c(motif)
    return(list(
      nodes = data.frame(id = 1L, x = 0, y = 0, label = labs[1], stringsAsFactors = FALSE),
      edges = data.frame(from = integer(0), to = integer(0), linetype = character(0), stringsAsFactors = FALSE)
    ))
  }
  if (identical(mt, "edge")) {
    if (length(labs) < 2) labs <- c(labs, rep("?", 2 - length(labs)))
    return(list(
      nodes = data.frame(
        id = c(1L, 2L),
        x = c(-0.25, 0.25),
        y = c(0, 0),
        label = labs[1:2],
        stringsAsFactors = FALSE
      ),
      edges = data.frame(from = 1L, to = 2L, linetype = "solid", stringsAsFactors = FALSE)
    ))
  }
  if (length(labs) < 3) labs <- c(labs, rep("?", 3 - length(labs)))

  if (identical(mt, "wedge")) {
    return(list(
      nodes = data.frame(
        id = c(1L, 2L, 3L),
        x = c(0, -0.25, 0.25),
        y = c(0.18, -0.18, -0.18),
        label = labs[1:3],
        stringsAsFactors = FALSE
      ),
      edges = data.frame(
        from = c(1L, 1L),
        to = c(2L, 3L),
        linetype = c("solid", "solid"),
        stringsAsFactors = FALSE
      )
    ))
  }

  tri_edges <- data.frame(
    from = c(1L, 1L, 2L),
    to = c(2L, 3L, 3L),
    linetype = c("solid", "solid", if (identical(mt, "triplet")) "22" else "solid"),
    stringsAsFactors = FALSE
  )
  list(
    nodes = data.frame(
      id = c(1L, 2L, 3L),
      x = c(0, -0.24, 0.24),
      y = c(0.25, -0.18, -0.18),
      label = labs[1:3],
      stringsAsFactors = FALSE
    ),
    edges = tri_edges
  )
}

build_side_panel_data <- function(tbl, label_map = NULL, side_text_mode = c("numbered", "verbose"), side_label_truncate = 3) {
  side_text_mode <- match.arg(side_text_mode)
  if (is.null(tbl) || !nrow(tbl)) {
    return(list(
      nodes = data.frame(),
      edges = data.frame(),
      labels = data.frame(),
      n = 0L
    ))
  }
  scale_xy <- 0.95
  k <- nrow(tbl)
  nodes_all <- list()
  edges_all <- list()
  label_all <- list()

  for (i in seq_len(k)) {
    row <- tbl[i, , drop = FALSE]
    lay <- motif_layout_template(row$motif, row$motif_type, side_label_truncate = side_label_truncate)
    yc <- k - i + 1

    nodes <- lay$nodes
    nodes$motif <- row$motif
    nodes$x <- nodes$x * scale_xy
    nodes$y <- yc + nodes$y * scale_xy
    # Keep text readable by moving top-node labels above and lower-node labels below.
    nodes$label_x <- nodes$x
    nodes$label_y <- nodes$y + ifelse(nodes$y >= yc, 0.11, -0.11)
    nodes_all[[i]] <- nodes

    edges <- lay$edges
    if (nrow(edges)) {
      edges$x <- nodes$x[match(edges$from, nodes$id)]
      edges$y <- nodes$y[match(edges$from, nodes$id)]
      edges$xend <- nodes$x[match(edges$to, nodes$id)]
      edges$yend <- nodes$y[match(edges$to, nodes$id)]
      edges$motif <- row$motif
      edges_all[[i]] <- edges
    }

    if (!is.null(label_map) && length(label_map) && !is.null(names(label_map))) {
      key_label <- as.character(label_map[[row$motif]])
      if (!length(key_label) || is.na(key_label) || key_label == "") key_label <- as.character(i)
    } else {
      key_label <- as.character(i)
    }
    if (identical(side_text_mode, "numbered")) {
      label_all[[i]] <- data.frame(
        x = 0.58,
        y = yc,
        label = key_label,
        stringsAsFactors = FALSE
      )
    } else {
      fdr_txt <- if (is.finite(row$FDR)) format(signif(row$FDR, 2), scientific = TRUE) else "NA"
      label_all[[i]] <- data.frame(
        x = 0.6,
        y = yc,
        label = paste0(row$motif, "  ", sprintf("%+.2f", row$logFC), "  FDR=", fdr_txt),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    nodes = do.call(rbind, nodes_all),
    edges = if (length(edges_all)) do.call(rbind, edges_all) else data.frame(),
    labels = do.call(rbind, label_all),
    n = k
  )
}

plot_motif_side_panel <- function(tbl, panel_title, label_map = NULL, side_text_mode = c("numbered", "verbose"), side_label_truncate = 3) {
  side_text_mode <- match.arg(side_text_mode)
  dd <- build_side_panel_data(
    tbl = tbl,
    label_map = label_map,
    side_text_mode = side_text_mode,
    side_label_truncate = side_label_truncate
  )
  if (!dd$n) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No motifs selected", size = 3.5, color = "#6b6b6b") +
        ggplot2::xlim(-1, 1) +
        ggplot2::ylim(-1, 1) +
        ggplot2::ggtitle(panel_title) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5))
    )
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dd$edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = linetype),
      color = "#455a64",
      linewidth = 0.7,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = dd$nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      fill = "#ffffff",
      color = "#263238",
      size = 4.1,
      stroke = 0.8
    ) +
    ggplot2::geom_text(
      data = dd$nodes,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      size = 2.85
    ) +
    ggplot2::scale_linetype_identity()

  if (identical(side_text_mode, "numbered")) {
    p <- p +
      ggplot2::geom_text(
        data = dd$labels,
        ggplot2::aes(x = x, y = y, label = label),
        hjust = 0,
        size = 3.8,
        fontface = "bold",
        color = "#111111"
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-0.45, 0.9),
        ylim = c(0.45, dd$n + 0.55),
        clip = "off"
      )
  } else {
    p <- p +
      ggplot2::geom_text(
        data = dd$labels,
        ggplot2::aes(x = x, y = y, label = label),
        hjust = 0,
        size = 3.05,
        color = "#263238"
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-0.45, 2.2),
        ylim = c(0.45, dd$n + 0.55),
        clip = "off"
      )
  }

  p +
    ggplot2::ggtitle(panel_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )
}

#' Volcano plot for motif differential results
#'
#' Plot differential motif results (from [top_edges()] or [top_triplets()]) in a
#' standard volcano layout, with flexible motif labeling options.
#'
#' @param x A results data frame (e.g., output of [top_edges()] or [top_triplets()]),
#'   or a `cellEdgeR_obj` containing `motif_edger()` results.
#' @param result_type Which results to fetch when `x` is a `cellEdgeR_obj`:
#'   `"edges"`, `"triplets_cooccurrence"`, or `"triplets_topology"`.
#' @param coef Coefficient name or index passed to [top_edges()] / [top_triplets()]
#'   when `x` is a `cellEdgeR_obj`.
#' @param model Model (`"full"` or `"null"`) used for `"edges"` and `"triplets_cooccurrence"`
#'   when `x` is a `cellEdgeR_obj`.
#' @param fdr_method Multiple testing correction method passed to `p.adjust` when needed.
#' @param fdr_cutoff FDR threshold used for point coloring.
#' @param logFC_cutoff Absolute logFC threshold used for point coloring.
#' @param label_mode Label selection mode: `"none"`, `"manual"`, `"top_n"`, `"top_bottom_n"`.
#' @param label_motifs Motif keys to label when `label_mode = "manual"`.
#' @param label_n Number of motifs to label for `"top_n"` and `"top_bottom_n"`.
#' @param label_top_by Ranking criterion for `"top_n"` labels: `"PValue"`, `"FDR"`, `"abs_logFC"`.
#' @param point_size Point size in the volcano.
#' @param point_alpha Point alpha in the volcano.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @return A `ggplot` object.
#' @export
plot_motif_volcano <- function(
  x,
  result_type = c("edges", "triplets_cooccurrence", "triplets_topology"),
  coef = NULL,
  model = c("full", "null"),
  fdr_method = "BH",
  fdr_cutoff = 0.05,
  logFC_cutoff = 0,
  label_mode = c("none", "manual", "top_n", "top_bottom_n"),
  label_motifs = NULL,
  label_n = 10,
  label_top_by = c("PValue", "FDR", "abs_logFC"),
  point_size = 1.8,
  point_alpha = 0.8,
  title = NULL,
  subtitle = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  result_type <- match.arg(result_type)
  model <- match.arg(model)
  label_mode <- match.arg(label_mode)
  label_top_by <- match.arg(label_top_by)
  tbl <- resolve_plot_results_table(x, result_type = result_type, coef = coef, model = model, fdr_method = fdr_method)
  df <- prepare_volcano_dataframe(tbl, fdr_cutoff = fdr_cutoff, logFC_cutoff = logFC_cutoff)
  labels <- select_volcano_labels(df, label_mode = label_mode, label_motifs = label_motifs, label_n = label_n, label_top_by = label_top_by)

  if (is.null(title)) {
    title <- switch(result_type,
      edges = "Motif volcano: node/edge",
      triplets_cooccurrence = "Motif volcano: triplet co-occurrence",
      triplets_topology = "Motif volcano: triplet topology"
    )
  }
  build_motif_volcano_plot(
    df = df,
    labels = labels,
    label_map = NULL,
    fdr_cutoff = fdr_cutoff,
    logFC_cutoff = logFC_cutoff,
    point_size = point_size,
    point_alpha = point_alpha,
    title = title,
    subtitle = subtitle
  )
}

#' Triptych volcano with ordered motif side panels
#'
#' Create a 3-panel figure with the most depleted motifs on the left, the volcano in the middle,
#' and the most enriched motifs on the right. Side motifs are displayed in fixed top-to-bottom order.
#'
#' @param x A results data frame (e.g., output of [top_edges()] or [top_triplets()]),
#'   or a `cellEdgeR_obj` containing `motif_edger()` results.
#' @param result_type Which results to fetch when `x` is a `cellEdgeR_obj`:
#'   `"edges"`, `"triplets_cooccurrence"`, or `"triplets_topology"`.
#' @param coef Coefficient name or index passed to [top_edges()] / [top_triplets()]
#'   when `x` is a `cellEdgeR_obj`.
#' @param model Model (`"full"` or `"null"`) used for `"edges"` and `"triplets_cooccurrence"`
#'   when `x` is a `cellEdgeR_obj`.
#' @param fdr_method Multiple testing correction method passed to `p.adjust` when needed.
#' @param n_side Number of automatically selected motifs on each side.
#' @param left_motifs Optional motif keys to force on the depleted (left) panel.
#' @param right_motifs Optional motif keys to force on the enriched (right) panel.
#' @param label_side_motifs Logical; if `TRUE`, motifs shown on side panels are labeled on the volcano.
#' @param side_text_mode Side-panel annotation mode. `"numbered"` (default) uses compact
#'   numeric IDs shared with volcano labels; `"verbose"` uses the previous motif/logFC/FDR text.
#' @param side_label_truncate Number of characters to keep for cell labels inside side-panel motifs.
#'   Set to `NA` for no truncation.
#' @param label_mode Additional volcano label mode: `"none"`, `"manual"`, `"top_n"`, `"top_bottom_n"`.
#' @param label_motifs Motif keys to label when `label_mode = "manual"`.
#' @param label_n Number of motifs to label for `"top_n"` and `"top_bottom_n"`.
#' @param label_top_by Ranking criterion for `"top_n"` labels: `"PValue"`, `"FDR"`, `"abs_logFC"`.
#' @param fdr_cutoff FDR threshold used for volcano point coloring.
#' @param logFC_cutoff Absolute logFC threshold used for volcano point coloring.
#' @param point_size Point size in the volcano.
#' @param point_alpha Point alpha in the volcano.
#' @param panel_widths Relative widths for left/center/right panels.
#' @param title Optional volcano title.
#' @param subtitle Optional volcano subtitle.
#' @return If package `patchwork` is installed, returns a combined patchwork plot.
#'   Otherwise draws the 3 panels and returns a list with components:
#'   `depleted_panel`, `volcano`, `enriched_panel`.
#' @export
plot_motif_volcano_triptych <- function(
  x,
  result_type = c("edges", "triplets_cooccurrence", "triplets_topology"),
  coef = NULL,
  model = c("full", "null"),
  fdr_method = "BH",
  n_side = 5,
  left_motifs = NULL,
  right_motifs = NULL,
  label_side_motifs = TRUE,
  side_text_mode = c("numbered", "verbose"),
  side_label_truncate = 3,
  label_mode = c("none", "manual", "top_n", "top_bottom_n"),
  label_motifs = NULL,
  label_n = 6,
  label_top_by = c("PValue", "FDR", "abs_logFC"),
  fdr_cutoff = 0.05,
  logFC_cutoff = 0,
  point_size = 1.8,
  point_alpha = 0.8,
  panel_widths = c(1.35, 2.4, 1.35),
  title = NULL,
  subtitle = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  result_type <- match.arg(result_type)
  model <- match.arg(model)
  side_text_mode <- match.arg(side_text_mode)
  label_mode <- match.arg(label_mode)
  label_top_by <- match.arg(label_top_by)
  trunc_is_na <- length(side_label_truncate) == 1 && is.na(side_label_truncate)
  trunc_is_num <- is.numeric(side_label_truncate) && length(side_label_truncate) == 1 &&
    is.finite(side_label_truncate) && side_label_truncate >= 1
  if (!(trunc_is_na || trunc_is_num)) {
    stop("side_label_truncate must be >= 1 or NA for no truncation.")
  }
  if (!is.numeric(panel_widths) || length(panel_widths) != 3 || any(!is.finite(panel_widths)) || any(panel_widths <= 0)) {
    stop("panel_widths must be a numeric vector of length 3 with positive values.")
  }

  tbl <- resolve_plot_results_table(x, result_type = result_type, coef = coef, model = model, fdr_method = fdr_method)
  df <- prepare_volcano_dataframe(tbl, fdr_cutoff = fdr_cutoff, logFC_cutoff = logFC_cutoff)
  side <- pick_side_motifs(df, n_side = n_side, left_motifs = left_motifs, right_motifs = right_motifs)

  extra_labels <- select_volcano_labels(
    df,
    label_mode = label_mode,
    label_motifs = label_motifs,
    label_n = label_n,
    label_top_by = label_top_by
  )

  side_motif_order <- unique(c(side$left$motif, side$right$motif))
  if (length(side_motif_order)) {
    side_label_map <- if (identical(side_text_mode, "numbered")) {
      stats::setNames(as.character(seq_along(side_motif_order)), side_motif_order)
    } else {
      stats::setNames(side_motif_order, side_motif_order)
    }
  } else {
    side_label_map <- character()
  }

  label_targets <- character()
  label_map <- character()
  if (isTRUE(label_side_motifs) && length(side_label_map)) {
    label_targets <- c(label_targets, names(side_label_map))
    label_map <- side_label_map
  }
  if (length(extra_labels)) {
    label_targets <- c(label_targets, extra_labels)
    extra_map <- stats::setNames(extra_labels, extra_labels)
    if (!length(label_map)) {
      label_map <- extra_map
    } else {
      new_extra <- setdiff(names(extra_map), names(label_map))
      if (length(new_extra)) label_map <- c(label_map, extra_map[new_extra])
    }
  }
  all_labels <- unique(label_targets)
  if (!length(label_map)) label_map <- NULL

  if (is.null(title)) {
    title <- switch(result_type,
      edges = "Motif volcano: node/edge",
      triplets_cooccurrence = "Motif volcano: triplet co-occurrence",
      triplets_topology = "Motif volcano: triplet topology"
    )
  }
  volcano_plot <- build_motif_volcano_plot(
    df = df,
    labels = all_labels,
    label_map = label_map,
    fdr_cutoff = fdr_cutoff,
    logFC_cutoff = logFC_cutoff,
    point_size = point_size,
    point_alpha = point_alpha,
    title = title,
    subtitle = subtitle
  )
  depleted_plot <- plot_motif_side_panel(
    tbl = side$left,
    panel_title = "Most depleted motifs",
    label_map = side_label_map,
    side_text_mode = side_text_mode,
    side_label_truncate = side_label_truncate
  )
  enriched_plot <- plot_motif_side_panel(
    tbl = side$right,
    panel_title = "Most enriched motifs",
    label_map = side_label_map,
    side_text_mode = side_text_mode,
    side_label_truncate = side_label_truncate
  )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(
      patchwork::wrap_plots(
        depleted_plot,
        volcano_plot,
        enriched_plot,
        nrow = 1,
        widths = panel_widths
      )
    )
  }

  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 1, ncol = 3, widths = grid::unit(panel_widths, "null"))
  vp <- grid::viewport(layout = lay)
  grid::pushViewport(vp)
  print(depleted_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(volcano_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(enriched_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid::popViewport()
  warning("Package patchwork is not installed; drew the triptych and returned component plots as a list.")
  invisible(list(
    depleted_panel = depleted_plot,
    volcano = volcano_plot,
    enriched_panel = enriched_plot
  ))
}
