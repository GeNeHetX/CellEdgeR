# Extract top simple motifs (nodes and edges)

Return volume-based differential results for node and edge motifs only.
Use
[`top_motifs_triplet`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs_triplet.md)
for 3-node motifs.

## Usage

``` r
top_motifs_simple(cellgraph, coef = NULL, model = c("full", "null"),
  n = Inf, fdr_method = "BH")
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` with `edger` results.

- coef:

  Coefficient name or index; defaults to the first non-intercept
  coefficient.

- model:

  Which stored model to use for edgeR strategies: `full` or `null`.

- n:

  Number of motifs to return; defaults to all.

- fdr_method:

  Multiple testing correction method for `p.adjust` (default `BH`).

## Value

A data frame with columns: motif, motif_type, logFC, PValue, FDR,
model_used.

## See also

[`top_motifs_triplet`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs_triplet.md),
[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3)
#> Erosion enabled but no boundary masks provided; counting all cells.
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
#> Fitting ancova models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Ancova dispersion could not be estimated; returning NA results.

top_motifs_simple(res, coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>   motif motif_type logFC PValue FDR model_used
#> 1   N_A       node    NA     NA  NA     volume
#> 2   N_B       node    NA     NA  NA     volume
#> 3 E_A_A       edge    NA     NA  NA     volume
#> 4 E_A_B       edge    NA     NA  NA     volume
#> 5 E_B_B       edge    NA     NA  NA     volume
```
