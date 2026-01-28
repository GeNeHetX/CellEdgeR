# Retrieve normalized counts as a single data frame

Pulls all normalized motifs from a cellgraph (`cellEdgeR_obj`) into one
data frame by row-binding layers.

## Usage

``` r
get_norm_counts(cellgraph, offset_mode = NULL, log2 = TRUE, pseudocount = 1)
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` returned by
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md).

- offset_mode:

  Which offset mode to pull normalized counts from; defaults to the
  first available when multiple are stored.

- log2:

  Logical; if `TRUE` (default), apply `log2(count + pseudocount)` to the
  normalized counts.

- pseudocount:

  Small positive value added before log transformation (default 1 to
  match offset construction).

## Value

A matrix/data frame with samples as rownames and motifs as column names.
