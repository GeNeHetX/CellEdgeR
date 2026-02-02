test_that("motif counting returns expected shapes and counts", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )

  graphs <- build_cell_graphs(cells, verbose = FALSE)
  res <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)

  expect_equal(res$label_levels, c("A", "B"))
  expect_equal(res$sample_name, c("s1", "s2"))

  y1 <- as.matrix(res$raw_count$node)
  expect_equal(drop(y1["N_A", "s1"]), 2)
  expect_equal(drop(y1["N_B", "s2"]), 2)

  y2 <- as.matrix(res$raw_count$edge)
  expect_equal(drop(y2["E_A_A", "s1"]), 1)
  expect_equal(drop(y2["E_A_B", "s1"]), 2)
  expect_equal(drop(y2["E_A_B", "s2"]), 2)
  expect_equal(drop(y2["E_B_B", "s2"]), 1)

  y3 <- as.matrix(res$raw_count$triangle)
  expect_equal(drop(y3["T_A_A_B", "s1"]), 1)
  expect_equal(drop(y3["T_A_B_B", "s2"]), 1)

  expect_equal(res$exposure$edges, c(s1 = 3, s2 = 3))
  expect_equal(res$exposure$triangles, c(s1 = 1, s2 = 1))
})

test_that("wedge can be returned alongside triangles", {
  cells <- list(
    s1 = data.frame(
      x = c(0, 2, 4, 1),
      y = c(0, 0, 0, 1),
      label = c("A", "B", "C", "A")
    )
  )

  graphs <- build_cell_graphs(cells, verbose = FALSE)
  res <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)

  expect_named(res$raw_count, c("node", "edge", "triangle", "wedge"))
  expect_setequal(rownames(res$raw_count$wedge), c("W_A_A_C", "W_B_A_C"))
  expect_equal(drop(as.matrix(res$raw_count$wedge)["W_B_A_C", "s1"]), 1)
  expect_equal(res$exposure$wedge, c(s1 = 2))
})

test_that("prebuilt graphs can be reused across motif runs", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)

  expect_true(inherits(graphs, "cellEdgeR_obj"))
  res_a <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
  res_b <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  res_positional <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)

  expect_equal(res_a$sample_name, res_b$sample_name)
  expect_equal(res_a$sample_name, res_positional$sample_name)
  expect_true(res_a$exposure$edges["s1"] <= res_b$exposure$edges["s1"])
})

test_that("edgeR pipeline runs with BH correction", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motif_obj$sample_name)

  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ condition",
    verbose = FALSE
  )

  expect_true(inherits(res, "cellEdgeR_obj"))
  expect_true(is.list(res$edger$strategies))
  expect_true(all(c("volume", "ancova") %in% names(res$edger$strategies)))
  tbl <- top_motifs_simple(res, coef = "conditiontreated")
  expect_true(is.data.frame(tbl))
  expect_true(all(c("motif", "motif_type", "logFC", "PValue", "FDR") %in% names(tbl)))
  expect_true(all(tbl$motif_type %in% c("node", "edge")))
})

test_that("motif_edger supports intercept-only design", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  sample_df <- data.frame(intercept = rep(1, length(motif_obj$sample_name)), row.names = motif_obj$sample_name)
  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ 1",
    verbose = FALSE
  )
  tbl <- top_motifs_simple(res, model = "null")
  expect_true(all(is.na(tbl$PValue) | is.numeric(tbl$PValue)))
})

test_that("get_motif_values returns motif and submotif values", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  tri_key <- rownames(motif_obj$raw_count$triangle)[1]
  vals <- get_motif_values(
    motif_obj,
    motif_key = tri_key,
    value = "raw"
  )
  expect_true(is.data.frame(vals))
  expect_equal(rownames(vals), motif_obj$sample_name)
  expect_true(tri_key %in% colnames(vals))

  norm_vals <- get_motif_values(
    motif_obj,
    motif_key = tri_key,
    value = "norm"
  )
  expect_true(tri_key %in% colnames(norm_vals))

  sub_vals <- get_motif_values(
    motif_obj,
    motif_key = tri_key,
    include_submotifs = TRUE,
    value = "raw"
  )
  expect_equal(colnames(sub_vals)[1], tri_key)
  expect_true(any(grepl("^E_", colnames(sub_vals))))
})

test_that("wedge is modeled when not merged", {
  cells <- list(
    s1 = data.frame(
      x = c(0, 2, 4, 1),
      y = c(0, 0, 0, 1),
      label = c("A", "B", "C", "A")
    ),
    s2 = data.frame(
      x = c(0, 2, 4, 1),
      y = c(0, 0, 0, 1),
      label = c("A", "B", "C", "A")
    )
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  sample_df <- data.frame(group = c("g1", "g2"), row.names = motif_obj$sample_name)
  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ group",
    verbose = FALSE
  )
  tbl <- top_motifs_triplet(res, strategy = "volume", coef = "groupg2")
  expect_true(any(tbl$motif_type == "wedge"))
  if (any(tbl$motif_type == "wedge")) {
    expect_true(all(grepl("^W_", tbl$motif[tbl$motif_type == "wedge"])))
  }
})

test_that("triplet_mode merge collapses 3-node motifs", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
    s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  sample_df <- data.frame(group = c("g1", "g2"), row.names = motif_obj$sample_name)
  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ group",
    triplet_mode = "merge",
    strategies = "volume",
    verbose = FALSE
  )
  expect_true(any(res$edger$motif_info$motif_type == "triplet"))
  expect_false(any(res$edger$motif_info$motif_type %in% c("triangle", "wedge")))
  tbl <- top_motifs_triplet(res, strategy = "volume", coef = "groupg2", triplet_mode = "merge")
  expect_true(all(tbl$motif_type == "triplet"))
})

test_that("triplet_mode closure uses triplet_force covariate", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
    s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  sample_df <- data.frame(group = c("g1", "g2"), row.names = motif_obj$sample_name)
  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ group",
    triplet_mode = "closure",
    strategies = "ancova",
    verbose = FALSE
  )
  expect_true("triplet_force" %in% res$edger$strategies$ancova$coef_names)
})

test_that("erosion drops boundary cells without retriangulating", {
  cells <- list(
    s1 = data.frame(
      x = c(0, 1, 0),
      y = c(0, 0, 1),
      label = c("A", "B", "C"),
      boundary = c(TRUE, FALSE, FALSE)
    )
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = NA_real_, verbose = FALSE)
  expect_equal(as.numeric(motif_obj$raw_count$node["N_A", "s1"]), 0)
  expect_equal(as.numeric(motif_obj$raw_count$node["N_B", "s1"]), 1)
  expect_equal(as.numeric(motif_obj$raw_count$node["N_C", "s1"]), 1)
  edge_key <- "E_B_C"
  if (edge_key %in% rownames(motif_obj$raw_count$edge)) {
    expect_equal(as.numeric(motif_obj$raw_count$edge[edge_key, "s1"]), 1)
  }
})

test_that("normalized counts follow motifs offsets", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  norm <- motif_obj$norm_counts$volume

  y_nodes <- as.matrix(motif_obj$raw_count$node)
  dge_nodes <- edgeR::calcNormFactors(edgeR::DGEList(counts = y_nodes), method = "TMM")
  norm_factors <- dge_nodes$samples$norm.factors
  if (!is.null(names(norm_factors))) norm_factors <- norm_factors[motif_obj$sample_name]
  expected_cells <- matrix(
    log(pmax(as.numeric(motif_obj$exposure$cells), 1)) + log(norm_factors),
    nrow = max(1, nrow(motif_obj$raw_count$node)),
    ncol = length(motif_obj$sample_name),
    byrow = TRUE,
    dimnames = list(rownames(motif_obj$raw_count$node), motif_obj$sample_name)
  )
  expect_equal(as.matrix(norm$node), as.matrix(motif_obj$raw_count$node) / exp(expected_cells))
})

test_that("motif_space_size reports combinatorial counts", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
    s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motif_obj <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  sample_df <- data.frame(group = c("g1", "g2"), row.names = motif_obj$sample_name)
  res <- motif_edger(
    cellgraph = motif_obj,
    sample_df = sample_df,
    design_formula = "~ group",
    triplet_mode = "merge",
    strategies = "volume",
    verbose = FALSE
  )
  space <- motif_space_size(res)
  expect_equal(space$labels, 3)
  expect_equal(space$triplet_mode, "merge")
  expect_true(isTRUE(space$include_wedge))
  expect_equal(space$counts$n_possible[space$counts$layer == "triplet"], choose(5, 3))
})

test_that("geometry-based triangulation tolerates duplicated and collinear points", {
  cells <- list(
    s1 = data.frame(
      x = c(0, 0, 1, 2, 2),
      y = c(0, 0, 0, 0, 1),
      label = c("A", "B", "A", "B", "C")
    ),
    s2 = data.frame(
      x = c(0, 1, 2, 2, 3),
      y = c(1, 1, 1, 2, 2),
      label = c("C", "B", "B", "A", "A")
    )
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  res <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
  expect_true(is.list(res))
  expect_true(all(res$exposure$edges >= 0))
  expect_true(all(res$exposure$triangles >= 0))
  expect_true(all(res$exposure$wedge >= 0))
})

test_that("motif objects can be merged to reproduce full counts", {
  demo <- make_demo_samples(n_cells = 30, n_labels = 3, n_samples = 6, seed = 123)
  group_a <- names(demo)[1:3]
  group_b <- names(demo)[4:6]
  demo[group_b] <- lapply(demo[group_b], function(df) {
    df$label <- ifelse(df$label == "c", "a", df$label)
    df
  })

  full <- count_motifs_graphs(
    build_cell_graphs(demo, verbose = FALSE),
    max_edge_len = NA_real_,
    include_wedge = TRUE,
    verbose = FALSE
  )
  obj_a <- count_motifs_graphs(
    build_cell_graphs(demo[group_a], verbose = FALSE),
    max_edge_len = NA_real_,
    include_wedge = TRUE,
    verbose = FALSE
  )
  obj_b <- count_motifs_graphs(
    build_cell_graphs(demo[group_b], verbose = FALSE),
    max_edge_len = NA_real_,
    include_wedge = TRUE,
    verbose = FALSE
  )

  merged <- merge_motif_objs(obj_a, obj_b, verbose = FALSE)

  expect_equal(merged$sample_name, full$sample_name)
  expect_equal(merged$label_levels, full$label_levels)
  expect_equal(merged$lab_to_id, full$lab_to_id)
  expect_equal(merged$raw_count, full$raw_count)
  expect_equal(merged$exposure, full$exposure)
  expect_equal(merged$offsets, full$offsets)
  expect_equal(merged$norm_counts, full$norm_counts)
  expect_equal(merged$relative_counts, full$relative_counts)
  expect_equal(merged$offset_part_id, full$offset_part_id)
  expect_equal(merged$offset_part_values, full$offset_part_values)
  expect_equal(merged$parameters, full$parameters)
  expect_equal(merged$per_sample_graph, full$per_sample_graph)
})
