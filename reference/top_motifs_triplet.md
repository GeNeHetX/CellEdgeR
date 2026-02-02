# Extract top 3-node motifs

Return differential results for 3-node motifs only. When
`triplet_mode = "merge"`, results are returned for unordered triplets.
Otherwise, triangles and wedges are returned separately (when
available).

## Usage

``` r
top_motifs_triplet(cellgraph, strategy = c("ancova", "volume"),
  coef = NULL, model = c("full", "null"), n = Inf, fdr_method = "BH",
  triplet_mode = NULL)
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` with `edger` results.

- strategy:

  Which strategy to use: `volume` or `ancova`.

- coef:

  Coefficient name or index; defaults to the first non-intercept
  coefficient.

- model:

  Which stored model to use for edgeR strategies: `full` or `null`.

- n:

  Number of motifs to return; defaults to all.

- fdr_method:

  Multiple testing correction method for `p.adjust` (default `BH`).

- triplet_mode:

  Triplet handling mode (`separate`, `merge`, or `closure`). Defaults to
  the mode stored in `cellgraph$edger$triplet_mode`.

## Value

A data frame with columns: motif, motif_type, logFC, PValue, FDR,
model_used.

## See also

[`top_motifs_simple`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs_simple.md),
[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
  s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE)
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
sample_df <- data.frame(group = c("g1", "g2"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ group", triplet_mode = "merge")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> edgeR GLM fit failed for volume/null model; tests will be empty.
#> Fitting ancova models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Ancova dispersion could not be estimated; returning NA results.

top_motifs_triplet(res, strategy = "volume", coef = "groupg2", triplet_mode = "merge")
#> Warning: No edgeR tests stored for volume coef: groupg2. Returning NA results.
#>      motif motif_type logFC PValue FDR model_used
#> 7 TP_A_B_C    triplet    NA     NA  NA     volume
```
