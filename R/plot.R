#' Boxplot of normalized motif counts
#'
#' Visualize normalized motif counts (e.g., from [normalize_motif_counts()]) across samples for a given motif.
#'
#' @param norm_counts List returned by [normalize_motif_counts()].
#' @param motif_key Character motif identifier (e.g., `"A_B_C"`).
#' @param layer Which layer to plot; one of `"size1"`, `"size2"`, `"size3"`, or `"wedges"` if present.
#' @param sample_df Optional data frame of sample metadata; rownames must match the motif columns.
#' @param group_var Optional column name in `sample_df` to use for grouping/coloring the boxplot.
#' @return A `ggplot` object.
#' @export
plot_motif_box <- function(norm_counts, motif_key, layer = c("size3", "size2", "size1", "wedges"),
                           sample_df = NULL, group_var = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  if (!is.list(norm_counts)) {
    stop("norm_counts must be the list returned by normalize_motif_counts().")
  }
  layer <- match.arg(layer, several.ok = FALSE)
  if (is.null(norm_counts[[layer]])) stop("Layer ", layer, " not found in norm_counts.")
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

#' Plot a sample graph with highlighted motifs
#'
#' Draw the Delaunay edges and cell coordinates for a single sample, optionally highlighting
#' a motif by labels (and edges connecting those labels).
#'
#' @param graph_obj Output of [build_cell_graphs()] with stored coordinates.
#' @param sample_id Sample name to plot.
#' @param max_edge_len Optional numeric threshold to prune long edges for display; set `Inf` to keep all.
#' @param highlight_labels Optional character vector of labels to emphasize.
#' @param motif_key Optional motif identifier (e.g., `"A_B_C"`); when provided, the labels in the key
#'   are highlighted and edges connecting those labels are accentuated.
#' @param motif_layer Layer for the motif key; one of `"size3"`, `"size2"`, `"size1"`, or `"wedges"`.
#' @param cells_by_sample Optional named list of raw sample data frames; used only when \code{graph_obj}
#'   lacks stored coordinates (backward compatibility with older objects).
#' @return A `ggplot` object.
#' @export
plot_sample_graph <- function(graph_obj, sample_id, max_edge_len = Inf, highlight_labels = NULL,
                              motif_key = NULL, motif_layer = c("size3", "size2", "size1", "wedges"),
                              cells_by_sample = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting. Install it with install.packages('ggplot2').")
  }
  if (!inherits(graph_obj, "cellEdgeR_graphs")) stop("graph_obj must come from build_cell_graphs().")
  if (!sample_id %in% graph_obj$samples) stop("sample_id not found in graph_obj.")
  ps <- graph_obj$per_sample[[sample_id]]
  if (is.null(ps$xy)) {
    if (!is.null(cells_by_sample) && sample_id %in% names(cells_by_sample)) {
      df <- standardize_sample_df(cells_by_sample[[sample_id]], sample_id)
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
  motif_layer <- match.arg(motif_layer)
  motif_labels <- NULL
  motif_pairs <- NULL
  if (!is.null(motif_key)) {
    motif_labels <- strsplit(motif_key, "_")[[1]]
    if (motif_layer == "size2") {
      motif_pairs <- paste(sort(motif_labels[seq_len(min(2, length(motif_labels)))]), collapse = "_")
    }
    highlight_labels <- unique(c(highlight_labels, motif_labels))
  }

  nodes <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    label = ps$labels_chr,
    highlight = if (is.null(highlight_labels)) FALSE else ps$labels_chr %in% highlight_labels,
    stringsAsFactors = FALSE
  )
  edges <- ps$edges
  if (length(edges) && is.finite(max_edge_len)) {
    keep <- ps$edge_len <= max_edge_len
    edges <- edges[keep, , drop = FALSE]
  }

  highlight_edges <- logical(nrow(edges))
  if (length(edges) && !is.null(motif_labels)) {
    la <- ps$labels_chr[edges[, 1]]
    lb <- ps$labels_chr[edges[, 2]]
    if (motif_layer == "size2" && !is.null(motif_pairs)) {
      pairs <- ifelse(la <= lb, paste(la, lb, sep = "_"), paste(lb, la, sep = "_"))
      highlight_edges <- pairs == motif_pairs
    } else {
      highlight_edges <- la %in% motif_labels & lb %in% motif_labels
    }
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
      alpha = 0.5
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
      ggplot2::aes(x = x, y = y, color = label, size = highlight, alpha = highlight)
    ) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.7, `TRUE` = 1), guide = "none") +
    ggplot2::ggtitle(paste0("Sample ", sample_id, " graph")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  p
}
