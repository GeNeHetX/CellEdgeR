# Extract top motifs from stored strategies

Returns a ranked data frame of motifs using strategy results stored by
[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md).

## Usage

``` r
top_motifs(cellgraph, strategy = c("hybrid", "ancova", "volume",
  "hierarchical"), coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH", append_strategies = NULL)
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` with `edger` results.

- strategy:

  Which strategy to use: `hybrid`, `ancova`, `volume`, or `hierarchical`
  (defaults to `hybrid`; hybrid applies volume to node/edge motifs and
  ancova to triangle/wedge motifs).

- coef:

  Coefficient name or index; defaults to the first non-intercept
  coefficient, or the intercept when only an intercept is present.

- model:

  Which stored model to use for edgeR strategies: `full` or `null`.

- n:

  Number of motifs to return; defaults to all.

- fdr_method:

  Multiple testing correction method for `p.adjust` (default `BH`).

- append_strategies:

  Optional vector of strategies to append raw PValues/logFC for. Use
  `TRUE` to append all other strategies.

## Value

A data frame with columns: motif, motif\\type, logFC, PValue, FDR, and
model\\used. Additional `logFC_<strategy>` and `PValue_<strategy>`
columns are added when requested. For `hybrid`, FDR is computed
separately for node/edge motifs and triangle/wedge motifs.

## See also

[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3)
#> Counting node motifs (cells by label)…
#> Counting edge motifs (unordered label pairs)…
#> Counting triangle motifs (unordered label triplets) in C++…
#> Counts ready: |labels|=2, singles=2, pairs=3, triangles=2
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ condition")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> Fitting edgeR (QL) for hierarchical offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for hierarchical/full model; tests will be empty.
#> Fitting ancova models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Ancova dispersion could not be estimated; returning NA results.
top_motifs(res, strategy = "volume", coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>     motif motif_type logFC PValue FDR model_used
#> 1     N_A       node    NA     NA  NA     volume
#> 2     N_B       node    NA     NA  NA     volume
#> 3   E_A_A       edge    NA     NA  NA     volume
#> 4   E_A_B       edge    NA     NA  NA     volume
#> 5   E_B_B       edge    NA     NA  NA     volume
#> 6 T_A_A_B   triangle    NA     NA  NA     volume
#> 7 T_A_B_B   triangle    NA     NA  NA     volume
```
