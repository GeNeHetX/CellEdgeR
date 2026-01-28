# Plot a sample graph with highlighted motifs

Draw the Delaunay edges and cell coordinates for a single sample,
optionally highlighting a motif by labels (and edges connecting those
labels).

## Usage

``` r
plot_sample_graph(
  graph_obj,
  sample_name,
  max_edge_len = Inf,
  highlight_labels = NULL,
  motif_key = NULL,
  motif_layer = NULL,
  cells_by_sample = NULL,
  motif_node_size = 3,
  dim_node_nonmotif = 0.6,
  alpha_node_nonmotif = 0.4,
  alpha_edge_nonmotif = 0.3
)
```

## Arguments

- graph_obj:

  Output of
  [`build_cell_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/build_cell_graphs.md)
  with stored coordinates.

- sample_name:

  Sample name to plot.

- max_edge_len:

  Optional numeric threshold to prune long edges for display; set `Inf`
  to keep all.

- highlight_labels:

  Optional character vector of labels to emphasize.

- motif_key:

  Optional motif identifier (e.g., `"E_A_B"`); when provided, the labels
  in the key are highlighted and edges connecting those labels are
  accentuated.

- motif_layer:

  Layer for the motif key; when `NULL`, inferred from the motif prefix.

- cells_by_sample:

  Optional named list of raw sample data frames; used only when
  `graph_obj` lacks stored coordinates (backward compatibility with
  older objects).

- motif_node_size:

  Size for nodes participating in the highlighted motif.

- dim_node_nonmotif:

  Factor to shrink nodes that are not in the highlighted motif.

- alpha_node_nonmotif:

  Alpha for nodes that are not in the highlighted motif.

- alpha_edge_nonmotif:

  Alpha for edges that are not in the highlighted motif.

## Value

A `ggplot` object.
