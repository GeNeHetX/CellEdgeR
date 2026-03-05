# Shared helpers for HochSchulz_2022_Melanoma validation scripts.
# Keep in sync with analysis/hochschulz_2022_melanoma_validation.Rmd.

load_hochschulz <- function(data_type, full_dataset) {
  options(timeout = max(600, getOption("timeout")))
  tryCatch(
    imcdatasets::HochSchulz_2022_Melanoma(data_type = data_type, full_dataset = full_dataset),
    error = function(e) {
      message("Initial download failed: ", conditionMessage(e))
      message("Clearing ExperimentHub cache for EH7824 and retrying...")
      try({
        eh <- ExperimentHub::ExperimentHub()
        try(ExperimentHub::removeCache(eh, "EH7824"), silent = TRUE)
        if (requireNamespace("BiocFileCache", quietly = TRUE)) {
          rid <- BiocFileCache::bfcquery(eh@bfc, "EH7824")
          if (nrow(rid)) BiocFileCache::bfcremove(eh@bfc, rid$rid)
        }
      }, silent = TRUE)
      imcdatasets::HochSchulz_2022_Melanoma(data_type = data_type, full_dataset = full_dataset)
    }
  )
}

first_matching <- function(cols, candidates) {
  idx <- match(tolower(candidates), tolower(cols))
  if (all(is.na(idx))) return(NA_character_)
  cols[idx[which(!is.na(idx))[1]]]
}

pick_group_col <- function(df, candidates) {
  for (cand in candidates) {
    if (!cand %in% names(df)) next
    vals <- df[[cand]]
    n_lev <- length(unique(vals))
    if (n_lev >= 2 && n_lev <= 6) return(cand)
  }
  NA_character_
}

prepare_cells_by_sample <- function(spe, label_col, image_col, x_col, y_col,
                                    min_cells, max_cells, max_images) {
  cd <- as.data.frame(colData(spe))
  if (!all(c(label_col, image_col) %in% names(cd))) {
    stop("Missing label/image columns. Set label_col and image_col explicitly.")
  }

  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords)) {
    if (!all(c(x_col, y_col) %in% names(cd))) {
      stop("Missing spatial coordinates. Set x_col/y_col explicitly.")
    }
    coords <- as.matrix(cd[, c(x_col, y_col)])
  }

  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    label = cd[[label_col]],
    image = cd[[image_col]],
    stringsAsFactors = FALSE
  )

  df <- df %>%
    group_by(image) %>%
    filter(n() >= min_cells) %>%
    ungroup()

  if (!is.null(max_cells)) {
    df <- df %>%
      group_by(image) %>%
      slice_sample(n = min(n(), max_cells)) %>%
      ungroup()
  }

  if (!is.null(max_images)) {
    keep <- df %>% distinct(image) %>% slice_sample(n = min(n(), max_images))
    df <- df %>% filter(image %in% keep$image)
  }

  cells_by_sample <- split(df[, c("x", "y", "label")], df$image)
  list(cells_by_sample = cells_by_sample, cell_df = df)
}

plot_volcano <- function(tbl, title) {
  ggplot(tbl, aes(x = logFC, y = -log10(PValue))) +
    geom_point(alpha = 0.5, size = 1) +
    labs(title = title, x = "logFC", y = "-log10(PValue)") +
    theme_minimal()
}

permute_graph_labels <- function(graph_obj, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- graph_obj
  out$per_sample_graph <- lapply(out$per_sample_graph, function(ps) {
    idx <- sample.int(ps$n)
    ps$labels_id <- ps$labels_id[idx]
    ps$labels_chr <- ps$labels_chr[idx]
    ps
  })
  out
}

build_triplet_cache <- function(counts) {
  tri <- counts$raw_count$triangle
  wedge <- counts$raw_count$wedge
  if (is.null(tri) || is.null(wedge)) return(NULL)
  tri_keys <- rownames(tri)
  wedge_keys <- rownames(wedge)
  build_triplet_map <- getFromNamespace("build_triplet_map", "CellEdgeR")
  collapse_triplet_counts <- getFromNamespace("collapse_triplet_counts", "CellEdgeR")
  collapse_triplet_offsets <- getFromNamespace("collapse_triplet_offsets", "CellEdgeR")
  triplet_map <- build_triplet_map(tri_keys, wedge_keys, prefix = "TP")
  trip_counts <- collapse_triplet_counts(tri, wedge, triplet_map, counts$sample_name)
  tri_off <- counts$offsets$volume$triangle
  wedge_off <- counts$offsets$volume$wedge
  trip_off <- collapse_triplet_offsets(tri_off, wedge_off, triplet_map, counts$sample_name)
  trip_norm <- trip_counts / exp(trip_off)
  list(counts = trip_counts, offsets = trip_off, norm = trip_norm)
}

run_parallel <- function(n_cores) {
  isTRUE(n_cores > 1) && .Platform$OS.type != "windows"
}

par_apply <- function(X, FUN, n_cores, ...) {
  if (run_parallel(n_cores)) {
    parallel::mclapply(X, FUN, mc.cores = n_cores, mc.set.seed = FALSE, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

safe_pull <- function(x, name) {
  if (is.list(x) && !is.null(x[[name]])) return(x[[name]])
  data.frame()
}

collect_layer_tables <- function(res, coef_name = NULL) {
  edges_all <- top_edges(res, coef = coef_name)
  nodes_tbl <- edges_all[edges_all$motif_type == "node", , drop = FALSE]
  edges_tbl <- edges_all[edges_all$motif_type == "edge", , drop = FALSE]
  co_tbl <- top_triplets(res, strategy = "cooccurrence", coef = coef_name)
  topo_tbl <- top_triplets(res, strategy = "topology", coef = coef_name)
  wedge_tbl <- topo_tbl[topo_tbl$motif_type == "wedge", , drop = FALSE]
  tri_tbl <- topo_tbl[topo_tbl$motif_type == "triangle", , drop = FALSE]
  list(
    node = nodes_tbl,
    edge = edges_tbl,
    triplet_cooccurrence = co_tbl,
    wedge = wedge_tbl,
    triangle = tri_tbl
  )
}

stack_layer_pvals <- function(layer_tbls, perm_id = NA_integer_) {
  out <- lapply(names(layer_tbls), function(layer) {
    tbl <- layer_tbls[[layer]]
    if (!nrow(tbl)) return(NULL)
    data.frame(
      layer = layer,
      motif = tbl$motif,
      logFC = tbl$logFC,
      PValue = tbl$PValue,
      perm = perm_id,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

make_qq_df <- function(pvals, layer) {
  p <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (!length(p)) return(data.frame())
  obs <- -log10(sort(p))
  exp <- -log10(ppoints(length(p)))
  data.frame(exp = exp, obs = obs, layer = layer)
}

subset_graphs_by_sample <- function(graphs, samples_keep) {
  samples_keep <- intersect(samples_keep, graphs$sample_name)
  graphs$sample_name <- samples_keep
  graphs$per_sample_graph <- graphs$per_sample_graph[samples_keep]
  graphs
}

sample_by_group <- function(sample_df, group_col, frac, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df <- sample_df
  df$sample <- rownames(df)
  df <- df %>%
    group_by(.data[[group_col]]) %>%
    slice_sample(n = max(1, floor(n() * frac))) %>%
    ungroup()
  df$sample
}

make_pair_df <- function(base_tbl, alt_tbl, layer) {
  if (!nrow(base_tbl) || !nrow(alt_tbl)) return(NULL)
  merged <- inner_join(
    base_tbl %>% select(motif, logFC, PValue),
    alt_tbl %>% select(motif, logFC, PValue),
    by = "motif",
    suffix = c("_base", "_alt")
  )
  merged$layer <- layer
  merged$logP_base <- -log10(merged$PValue_base)
  merged$logP_alt <- -log10(merged$PValue_alt)
  merged
}

summarize_pairs <- function(pair_df, top_n = 20) {
  if (is.null(pair_df) || !nrow(pair_df)) return(NULL)
  top_base <- head(pair_df$motif[order(pair_df$PValue_base, na.last = TRUE)], top_n)
  top_alt <- head(pair_df$motif[order(pair_df$PValue_alt, na.last = TRUE)], top_n)
  denom <- length(union(top_base, top_alt))
  top_jacc <- if (denom > 0) length(intersect(top_base, top_alt)) / denom else NA_real_
  data.frame(
    layer = unique(pair_df$layer),
    n = nrow(pair_df),
    cor_logFC = suppressWarnings(cor(pair_df$logFC_base, pair_df$logFC_alt, use = "complete.obs")),
    cor_logP = suppressWarnings(cor(pair_df$logP_base, pair_df$logP_alt, use = "complete.obs")),
    top_jaccard = top_jacc,
    stringsAsFactors = FALSE
  )
}

compare_layer_sets <- function(base_layers, alt_layers, top_n = 20) {
  pair_list <- list()
  summary_list <- list()
  for (layer in names(base_layers)) {
    pair_df <- make_pair_df(base_layers[[layer]], alt_layers[[layer]], layer)
    if (is.null(pair_df) || !nrow(pair_df)) next
    pair_list[[layer]] <- pair_df
    summary_list[[layer]] <- summarize_pairs(pair_df, top_n = top_n)
  }
  list(
    pairs = bind_rows(pair_list),
    summary = bind_rows(summary_list)
  )
}

crop_cells_by_space <- function(cells_by_sample, keep_frac) {
  keep_frac <- max(min(keep_frac, 1), 0)
  if (keep_frac == 1) return(cells_by_sample)
  q <- (1 - keep_frac) / 2
  lapply(cells_by_sample, function(df) {
    if (!nrow(df)) return(df)
    x_lo <- as.numeric(stats::quantile(df$x, probs = q, na.rm = TRUE))
    x_hi <- as.numeric(stats::quantile(df$x, probs = 1 - q, na.rm = TRUE))
    y_lo <- as.numeric(stats::quantile(df$y, probs = q, na.rm = TRUE))
    y_hi <- as.numeric(stats::quantile(df$y, probs = 1 - q, na.rm = TRUE))
    df[df$x >= x_lo & df$x <= x_hi & df$y >= y_lo & df$y <= y_hi, , drop = FALSE]
  })
}

random_subsample_cells <- function(cells_by_sample, keep_frac, seed = NULL) {
  keep_frac <- max(min(keep_frac, 1), 0)
  if (keep_frac == 1) return(cells_by_sample)
  if (!is.null(seed)) set.seed(seed)
  lapply(cells_by_sample, function(df) {
    if (!nrow(df)) return(df)
    n_keep <- max(1, floor(nrow(df) * keep_frac))
    df %>% slice_sample(n = min(n_keep, nrow(df)))
  })
}

filter_cells_by_min <- function(cells_by_sample, min_cells) {
  keep <- vapply(cells_by_sample, nrow, integer(1)) >= min_cells
  cells_by_sample[keep]
}
