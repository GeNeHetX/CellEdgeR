# Differential motif testing with edgeR

Fit a single volume-offset edgeR QL model across edge, wedge, and
triangle motifs. Node counts are used internally to build graph volumes
but are not tested by `motif_edger()`. Wedges and triangles are kept as
separate motifs; there is no merged triplet test. Low-count motifs are
filtered with
[`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
before model fitting, so motifs carried by too few samples are not
tested.

## Usage

``` r
motif_edger(cellgraph, sample_df, design_formula, verbose = TRUE)
```

## Arguments

- cellgraph:

  Output list from
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md).

- sample_df:

  Data frame with sample metadata; rownames must match
  `cellgraph$sample_name`.

- design_formula:

  Formula or formula string passed to `model.matrix`, e.g.
  `~ condition + batch`.

- verbose:

  Logical; print progress.

## Value

The input `cellgraph` augmented with `edger`, containing one stored
strategy named `volume`. Use
[`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
for edge results and
[`top_motifs2()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
for wedge/triangle results with edge-submotif statistics. Filter details
are stored in `cellgraph$edger$filter`.

## Details

The model is volume based by default. It stacks edge, triangle, and
wedge motifs and fits edgeR quasi-likelihood tests using the volume
offsets stored by
[`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md).
The output result sets are `edges` and `motifs2`. Motifs filtered by
[`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
do not receive p-values or FDR values.

## See also

[`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md),
[`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md),
[`top_motifs2()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)

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
head(top_edges(res, coef = "conditiontreated"))
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#> [1] motif      motif_type logFC      PValue     FDR        model_used
#> <0 rows> (or 0-length row.names)
head(top_motifs2(res, coef = "conditiontreated"))
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>  [1] motif                  motif_type             logFC                 
#>  [4] PValue                 FDR                    model_used            
#>  [7] edge12                 edge13                 edge23                
#> [10] edge12_logFC           edge13_logFC           edge23_logFC          
#> [13] edge12_PValue          edge13_PValue          edge23_PValue         
#> [16] edge12_FDR             edge13_FDR             edge23_FDR            
#> [19] submotif_min_FDR       submotif_max_abs_logFC
#> <0 rows> (or 0-length row.names)
```
