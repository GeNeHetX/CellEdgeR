# Boxplot of normalized motif counts (internal helper)

Visualize normalized motif counts (e.g., from
[`get_norm_counts()`](https://GeNeHetX.github.io/CellEdgeR/reference/get_norm_counts.md))
across samples for a given motif.

## Usage

``` r
plot_motif_box(
  norm_counts,
  motif_key,
  layer = NULL,
  sample_df = NULL,
  group_var = NULL,
  offset_mode = NULL
)
```

## Arguments

- norm_counts:

  List of normalized counts (e.g., `cellgraph$norm_counts` from
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md));
  when multiple offset modes are present, select one via `offset_mode`.

- motif_key:

  Character motif identifier with prefix (e.g., `"E_A_B"`).

- layer:

  Which layer to plot; when `NULL`, inferred from the motif prefix
  (`N_`, `E_`, `T_`, or `W_`).

- sample_df:

  Optional data frame of sample metadata; rownames must match the motif
  columns.

- group_var:

  Optional column name in `sample_df` to use for grouping/coloring the
  boxplot.

## Value

A `ggplot` object.
