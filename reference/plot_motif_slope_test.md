# Slope test for geometric scaling

Compare observed triangle counts to the structural expectation from
sub-motif offsets using a log-log regression.

## Usage

``` r
plot_motif_slope_test(
  cellgraph,
  offset_mode = "hier_null",
  layer = "triangle",
  motif_key = NULL,
  log_base = exp(1),
  pseudocount = NULL
)
```

## Arguments

- cellgraph:

  Output of
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md)
  (or
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md))
  with stored offsets.

- offset_mode:

  Offset set to use for the expectation; defaults to `"hier_null"`.

- layer:

  Motif layer to plot; defaults to `"triangle"`.

- motif_key:

  Optional character vector of motif keys to include; defaults to all
  motifs.

- log_base:

  Logarithm base for the regression; defaults to natural log.

- pseudocount:

  Positive value added before log transform; defaults to
  `cellgraph$parameters$offset_pseudo` or 1.

## Value

A `ggplot` object with the fitted slope annotation.
