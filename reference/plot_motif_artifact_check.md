# Artifact check: residuals versus edge density

Plot per-sample log residuals from hierarchical offsets against total
edge density.

## Usage

``` r
plot_motif_artifact_check(
  cellgraph,
  offset_mode = "hier_null",
  layer = "triangle",
  motif_key = NULL,
  log_base = 2,
  pseudocount = NULL,
  smooth = TRUE
)
```

## Arguments

- cellgraph:

  Output of
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md)
  (or
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md))
  with stored relative counts.

- offset_mode:

  Offset set to use for residuals; defaults to `"hier_null"`.

- layer:

  Motif layer to plot; defaults to `"triangle"`.

- motif_key:

  Optional character vector of motif keys to include; defaults to all
  motifs.

- log_base:

  Logarithm base for the x-axis; defaults to 2 to match residuals.

- pseudocount:

  Positive value added before log transform; defaults to
  `cellgraph$parameters$offset_pseudo` or 1.

- smooth:

  Logical; add a linear trend line.

## Value

A `ggplot` object.
