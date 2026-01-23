skip_if_not_installed("ggplot2")

test_that("plot_motif_box returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motifs <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  norm <- motifs$norm_counts$volume
  p <- plot_motif_box(norm, motif_key = "T_A_A_B")
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_graph returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  p <- plot_sample_graph(graphs, sample_name = "s1", motif_key = "T_A_A_B")
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_graph can recover coordinates when missing", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  graphs$per_sample_graph[[1]]$xy <- NULL
  p <- plot_sample_graph(graphs, sample_name = "s1", motif_key = "T_A_A_B",
    cells_by_sample = cells)
  expect_s3_class(p, "ggplot")
})

test_that("plot_motif_slope_test returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motifs <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  p <- plot_motif_slope_test(motifs)
  expect_s3_class(p, "ggplot")
})

test_that("plot_motif_artifact_check returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  motifs <- count_motifs_graphs(graphs, max_edge_len = 3, verbose = FALSE)
  p <- plot_motif_artifact_check(motifs)
  expect_s3_class(p, "ggplot")
})
