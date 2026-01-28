# Build cached cell graphs for motif counting

Builds a constrained Delaunay triangulation for each sample and caches
the raw edges so multiple motif-counting passes can reuse the same
triangulation.

## Usage

``` r
build_cell_graphs(cells_by_sample, n_cores = 1, verbose = TRUE, max_labels = 100)
```

## Arguments

- cells_by_sample:

  Named list of data frames where the first two columns are numeric
  `x, y` coordinates and the third column holds the cell-type label.
  Names of the list are treated as sample names.

- n_cores:

  Integer number of cores to request; values greater than 1 trigger
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix-alikes, while 1 runs sequentially.

- verbose:

  Logical; print progress messages for standardization and graph
  construction.

- max_labels:

  Maximum number of unique labels allowed across all samples (fails
  early if exceeded).

## Details

This helper performs the expensive Delaunay triangulation step and
records every edge's length. Pass the returned `cellEdgeR_graphs` object
to `count_motifs_graphs(graph_obj = ...)` to compute motif counts
without re-running the triangulation, even when the pruning or wedge
options change.

## Value

A list of class `cellEdgeR_graphs` containing:

- `sample_name`:

  Character vector of sample names.

- `label_levels/lab_to_id`:

  Sorted unique cell-type labels and a global label-to-integer lookup.

- `per_sample_graph`:

  Per-sample entries storing `n` (node count), `labels_chr`,
  `labels_id`, `edges`, edge lengths `edge_len`, and coordinates `xy`.

- `raw_count/exposure/offsets/norm_counts/relative_counts`:

  Empty placeholders until
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md)
  is called.

- `edger`:

  Empty placeholder until
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
  is called.

- `parameters`:

  Build parameters (placeholders for `max_edge_len`, `offset_pseudo`,
  `built_from`).

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
```
