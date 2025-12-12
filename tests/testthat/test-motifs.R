test_that("motif counting returns expected shapes and counts", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )

  res <- count_motifs_graphs(cells_by_sample = cells, max_edge_len = 3, verbose = FALSE)

  expect_equal(res$label_levels, c("A", "B"))
  expect_equal(res$samples, c("s1", "s2"))

  y1 <- as.matrix(res$counts$size1)
  expect_equal(drop(y1["A", "s1"]), 2)
  expect_equal(drop(y1["B", "s2"]), 2)

  y2 <- as.matrix(res$counts$size2)
  expect_equal(drop(y2["A_A", "s1"]), 1)
  expect_equal(drop(y2["A_B", "s1"]), 2)
  expect_equal(drop(y2["A_B", "s2"]), 2)
  expect_equal(drop(y2["B_B", "s2"]), 1)

  y3 <- as.matrix(res$counts$size3)
  expect_equal(drop(y3["A_A_B", "s1"]), 1)
  expect_equal(drop(y3["A_B_B", "s2"]), 1)

  expect_equal(res$exposure$edges, c(s1 = 3, s2 = 3))
  expect_equal(res$exposure$triangles, c(s1 = 1, s2 = 1))
})

test_that("wedges can be returned alongside triangles", {
  cells <- list(
    s1 = data.frame(
      x = c(0, 2, 4, 1),
      y = c(0, 0, 0, 1),
      label = c("A", "B", "C", "A")
    )
  )

  res <- count_motifs_graphs(cells, max_edge_len = NA_real_, include_wedges = TRUE, verbose = FALSE)

  expect_named(res$counts, c("size1", "size2", "size3", "wedges"))
  expect_setequal(rownames(res$counts$wedges), c("A_A_C", "B_A_C"))
  expect_equal(drop(as.matrix(res$counts$wedges)["B_A_C", "s1"]), 1)
  expect_equal(res$exposure$wedges, c(s1 = 2))
})

test_that("prebuilt graphs can be reused across motif runs", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)

  expect_s3_class(graphs, "cellEdgeR_graphs")
  res_a <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
  res_b <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedges = TRUE, verbose = FALSE)
  res_positional <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)

  expect_equal(res_a$samples, res_b$samples)
  expect_equal(res_a$samples, res_positional$samples)
  expect_true(res_a$exposure$edges["s1"] <= res_b$exposure$edges["s1"])
})

test_that("edgeR pipeline runs with DAG correction", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  motif_obj <- count_motifs_graphs(cells, max_edge_len = 3, verbose = FALSE)
  sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motif_obj$samples)

  res <- motif_edger(
    motif_obj = motif_obj,
    sample_df = sample_df,
    design_formula = "~ condition",
    coef = "conditiontreated",
    verbose = FALSE
  )

  expect_named(res$results, c("size1", "size2", "size3"))
  check_cols <- function(df) is.null(df) || all(c("motif", "logFC", "PValue", "FDR_BH") %in% names(df))
  expect_true(all(vapply(res$results, check_cols, logical(1))))
  expect_equal(res$fdr$method, "bh")
  expect_null(res$dag)

  res_dag <- motif_edger(
    motif_obj = motif_obj,
    sample_df = sample_df,
    design_formula = "~ condition",
    coef = "conditiontreated",
    verbose = FALSE,
    fdr_method = "dagger"
  )

  expect_equal(res_dag$fdr$method, "dagger")
  expect_true(is.list(res_dag$dag))
  expect_true(is.matrix(res_dag$dag$edges) || length(res_dag$dag$edges) == 0)
  expect_true(all(vapply(res_dag$results, function(df) is.null(df) || "FDR_DAG" %in% names(df), logical(1))))
})

test_that("normalized counts follow motifs offsets", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  motif_obj <- count_motifs_graphs(cells, max_edge_len = 3, verbose = FALSE)
  norm <- normalize_motif_counts(motif_obj, pseudo = 0.5)

  expected_cells <- matrix(
    log(pmax(as.numeric(motif_obj$exposure$cells), 1)),
    nrow = max(1, nrow(motif_obj$counts$size1)),
    ncol = length(motif_obj$samples),
    byrow = TRUE,
    dimnames = list(rownames(motif_obj$counts$size1), motif_obj$samples)
  )
  expect_equal(as.matrix(norm$size1), as.matrix(motif_obj$counts$size1) / exp(expected_cells))
})

test_that("normalize_motif_counts can return likelihood ratios", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  motif_obj <- count_motifs_graphs(cells, max_edge_len = 3, verbose = FALSE)
  norm <- normalize_motif_counts(motif_obj, pseudo = 0.5, return_lr = TRUE)
  expect_true(is.list(norm))
  expect_equal(dim(norm$size3), dim(motif_obj$counts$size3))
  expect_true(all(is.finite(as.matrix(norm$size3))))
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
  res <- count_motifs_graphs(cells, max_edge_len = NA_real_, include_wedges = TRUE, verbose = FALSE)
  expect_true(is.list(res))
  expect_true(all(res$exposure$edges >= 0))
  expect_true(all(res$exposure$triangles >= 0))
  expect_true(all(res$exposure$wedges >= 0))
})
