# Top motif results

Convenience helpers for retrieving ranked motif results from
[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md).
`top_edges()` returns edge motifs only. `top_motifs2()` returns wedge
and triangle motifs with the corresponding edge-submotif statistics
added as columns. `motif_results()` returns both tables as a named list.
`top_triplets()` is a deprecated compatibility alias for
`top_motifs2()`.

## Usage

``` r
top_edges(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_motifs2(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_triplets(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH", strategy = NULL)

motif_results(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` with `edger` results.

- coef:

  Coefficient name or index; defaults to the first non-intercept
  coefficient.

- model:

  Which stored model to use: `full` or `null`.

- n:

  Number of motifs to return; defaults to all.

- fdr_method:

  Multiple testing correction method for `p.adjust` (default `BH`).

- strategy:

  Deprecated compatibility argument for `top_triplets()`; ignored.

## Value

A data frame with columns `motif`, `motif_type`, `logFC`, `PValue`,
`FDR`, and `model_used`. `top_motifs2()` also includes `edge12`,
`edge13`, `edge23` and their edge-level `logFC`, `PValue`, and `FDR`
values.

## See also

[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md),
[`motif_space_size()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_space_size.md)

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3, include_wedge = TRUE)
#> Erosion enabled but no boundary masks provided; counting all cells.
#> Counting node motifs (cells by label)…
#> Counting edge motifs (unordered label pairs)…
#> Counting triangle motifs (unordered label triplets) in C++…
#> Also collecting wedge (open triplets)…
#> Counts ready: |labels|=2, singles=2, pairs=3, triangles=2, wedge=0
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)

res <- motif_edger(motifs, sample_df, "~ condition")
#> Filtering low-count motifs with edgeR::filterByExpr: kept 0 / 5.
top_edges(res, coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#> [1] motif      motif_type logFC      PValue     FDR        model_used
#> <0 rows> (or 0-length row.names)
top_motifs2(res, coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>  [1] motif                  motif_type             logFC                 
#>  [4] PValue                 FDR                    model_used            
#>  [7] edge12                 edge13                 edge23                
#> [10] edge12_logFC           edge13_logFC           edge23_logFC          
#> [13] edge12_PValue          edge13_PValue          edge23_PValue         
#> [16] edge12_FDR             edge13_FDR             edge23_FDR            
#> [19] submotif_min_FDR       submotif_max_abs_logFC
#> <0 rows> (or 0-length row.names)
motif_results(res, coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#> $edges
#> [1] motif      motif_type logFC      PValue     FDR        model_used
#> <0 rows> (or 0-length row.names)
#> 
#> $motifs2
#>  [1] motif                  motif_type             logFC                 
#>  [4] PValue                 FDR                    model_used            
#>  [7] edge12                 edge13                 edge23                
#> [10] edge12_logFC           edge13_logFC           edge23_logFC          
#> [13] edge12_PValue          edge13_PValue          edge23_PValue         
#> [16] edge12_FDR             edge13_FDR             edge23_FDR            
#> [19] submotif_min_FDR       submotif_max_abs_logFC
#> <0 rows> (or 0-length row.names)
#> 
```
