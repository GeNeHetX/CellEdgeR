skip_if_not_installed("ggplot2")

test_that("plot_motif_box returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  motifs <- count_motifs_graphs(cells, max_edge_len = 3, verbose = FALSE)
  norm <- normalize_motif_counts(motifs, pseudo = 0.5)
  p <- plot_motif_box(norm, motif_key = "A|A|B", layer = "size3")
  expect_s3_class(p, "ggplot")
})

test_that("plot_sample_graph returns a ggplot", {
  cells <- list(
    s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
    s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
  )
  graphs <- build_cell_graphs(cells, verbose = FALSE)
  p <- plot_sample_graph(graphs, sample_id = "s1", motif_key = "A|A|B", motif_layer = "size3")
  expect_s3_class(p, "ggplot")
})
