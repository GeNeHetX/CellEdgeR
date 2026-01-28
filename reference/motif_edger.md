# Differential motif testing with multiple strategies

Fits multiple edgeR-based strategies across stacked motifs and stores
results for later extraction.

## Usage

``` r
motif_edger(cellgraph, sample_df, design_formula, verbose = TRUE,
  merge_triplets = FALSE, strategies = c("volume", "hierarchical", "ancova"))
```

## Arguments

- cellgraph:

  Output of
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md),
  containing motif counts and offsets.

- sample_df:

  Data frame with sample metadata; row names must match
  `cellgraph\$sample_name`.

- design_formula:

  Formula string compatible with
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  (e.g. `~ condition + batch`).

- verbose:

  Logical; print progress for edgeR.

- merge_triplets:

  If `TRUE`, triangles and wedges are merged under a shared TW prefix
  before fitting.

- strategies:

  Character vector of strategies to run; defaults to `volume`,
  `hierarchical`, and `ancova`.

## Details

Three strategies are supported by default:

- `volume`: volume offsets (Chung-Lu baseline)

- `hierarchical`: sub-motif offsets using fitted edges (hier\\null)

- `ancova`: volume offsets plus edge-force covariates per triangle/wedge

Use
[`top_motifs`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs.md)
to extract ranked tables with logFC and FDR.

## Value

The input cellgraph augmented with `edger`, containing:

- `strategies`:

  Named list of strategy results (`volume`, `hierarchical`, `ancova`).

- `motif_info`:

  Data frame of motif keys and motif types for joins.

- `sample_df`:

  Sample metadata used to build the model matrix.

- `merge_triplets`:

  Whether triangles+wedge were merged.

## See also

[`count_motifs_graphs`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md),
[`top_motifs`](https://GeNeHetX.github.io/CellEdgeR/reference/top_motifs.md)

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
head(top_motifs(res, strategy = "volume", coef = "conditiontreated"))
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>     motif motif_type logFC PValue FDR model_used
#> 1     N_A       node    NA     NA  NA     volume
#> 2     N_B       node    NA     NA  NA     volume
#> 3   E_A_A       edge    NA     NA  NA     volume
#> 4   E_A_B       edge    NA     NA  NA     volume
#> 5   E_B_B       edge    NA     NA  NA     volume
#> 6 T_A_A_B   triangle    NA     NA  NA     volume
```
