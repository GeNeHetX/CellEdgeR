# Count the number of possible tested motif labels

Compute the combinatorial number of possible edge, wedge, and triangle
label motifs implied by a `cellEdgeR_obj`. Counts reflect label
combinations, not graph isomorphism classes. Node motifs and merged
triplets are not part of the tested motif space.

## Usage

``` r
motif_space_size(cellgraph, include_wedge = NULL)
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` (typically after
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)).

- include_wedge:

  Logical; whether to include wedge motifs in the count. Defaults to
  `cellgraph$parameters$include_wedge` when available.

## Value

A list with `labels`, `include_wedge`, `counts` (data frame), and
`total` (sum of possible tested motifs).

## See also

[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md),
[`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md),
[`top_motifs2()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
  s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
)
obj <- count_motifs_graphs(build_cell_graphs(cells), include_wedge = TRUE)
#> Samples: s1, s2
#> Building Delaunay graphs (sequential)…
#> Erosion enabled but no boundary masks provided; counting all cells.
#> Edge pruning disabled (max_edge_len is NA/NULL/<= 0).
#> Counting node motifs (cells by label)…
#> Counting edge motifs (unordered label pairs)…
#> Counting triangle motifs (unordered label triplets) in C++…
#> Also collecting wedge (open triplets)…
#> Counts ready: |labels|=3, singles=3, pairs=3, triangles=1, wedge=0
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: edgeR dispersion estimation failed for model residuals; returning zeros.
meta <- data.frame(group = c("g1", "g2"), row.names = obj$sample_name)
obj <- motif_edger(obj, meta, "~ group")
#> Filtering low-count motifs with edgeR::filterByExpr: kept 0 / 4.

motif_space_size(obj)
#> $labels
#> [1] 3
#> 
#> $include_wedge
#> [1] TRUE
#> 
#> $counts
#>      layer n_possible
#> 1     edge          6
#> 2 triangle         10
#> 3    wedge         18
#> 
#> $total
#> [1] 34
#> 
```
