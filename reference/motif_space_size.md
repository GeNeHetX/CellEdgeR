# Count possible motif labels

Compute the combinatorial number of possible motif labelings implied by
a `cellEdgeR_obj`. Counts reflect label combinations, not graph
isomorphism classes.

## Usage

``` r
motif_space_size(cellgraph, triplet_mode = NULL, include_wedge = NULL)
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` (typically after
  [`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)).

- triplet_mode:

  Triplet handling mode (`separate`, `merge`, or `closure`). Defaults to
  the mode stored in `cellgraph$edger$triplet_mode`.

- include_wedge:

  Logical; whether to include wedge motifs in the count. Defaults to
  `cellgraph$parameters$include_wedge` when available.

## Value

A list with `labels`, `triplet_mode`, `include_wedge`, `counts` (a data
frame of layers and counts), and `total` (sum of possible motifs).

## See also

[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md),
[`top_motifs_simple`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs_simple.md),
[`top_motifs_triplet`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs_triplet.md)

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
#> Warning: edgeR fit failed for edge-derived offsets (null); falling back to counts.
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: No residual df: setting dispersion to NA
#> Warning: edgeR dispersion estimation failed for model residuals; returning zeros.
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: edgeR fit failed for model residuals; returning zeros.
#> Warning: No residual df: setting dispersion to NA
#> Warning: edgeR dispersion estimation failed for model residuals; returning zeros.
meta <- data.frame(group = c("g1", "g2"), row.names = obj$sample_name)
obj <- motif_edger(obj, meta, "~ group", triplet_mode = "merge")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> edgeR GLM fit failed for volume/null model; tests will be empty.
#> Fitting ancova models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Ancova dispersion could not be estimated; returning NA results.

motif_space_size(obj)
#> $labels
#> [1] 3
#> 
#> $triplet_mode
#> [1] "merge"
#> 
#> $include_wedge
#> [1] TRUE
#> 
#> $counts
#>     layer n_possible
#> 1    node          3
#> 2    edge          6
#> 3 triplet         10
#> 
#> $total
#> [1] 19
#> 
```
