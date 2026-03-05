# Differential motif testing with multiple strategies

Fits multiple edgeR-based strategies across stacked motifs and stores
results for later extraction.

## Usage

``` r
motif_edger(cellgraph, sample_df, design_formula, verbose = TRUE,
  triplet_mode = c("both", "separate", "merge", "closure"),
  strategies = c("volume", "submotif_adj"))
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

- triplet_mode:

  How to handle 3-node motifs: `both` computes co-occurrence (merge) and
  topology (closure) in a single call (default). `separate` keeps
  triangle and wedge motifs, `merge` combines wedges+triangles into
  unordered triplet motifs, and `closure` models triangles with
  total-triplet adjustment. `merge` and `closure` require
  `count_motifs_graphs(..., include_wedge = TRUE)`.

- strategies:

  Character vector of strategies to run; defaults to `volume` and
  submotif-adjusted (`submotif_adj`).

## Details

Two strategies are supported by default:

- `volume`: volume offsets (Chung-Lu baseline)

- `submotif_adj`: submotif-adjusted model using volume offsets plus
  lower-order covariates for 3-node motifs

Use
[`top_edges`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
for node/edge motifs and
[`top_triplets`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
for 3-node motifs.

## Value

The input cellgraph augmented with `edger`, containing:

- `strategies`:

  Named list of strategy results (`volume`, `submotif_adj`).

- `motif_info`:

  Data frame of motif keys and motif types for joins.

- `sample_df`:

  Sample metadata used to build the model matrix.

- `triplet_mode`:

  Triplet handling mode used for 3-node motifs (may be a named list when
  `both`).

## See also

[`count_motifs_graphs`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md),
[`top_edges`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md),
[`top_triplets`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)

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
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ condition")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> Fitting submotif-adjusted models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Submotif-adjusted dispersion could not be estimated; returning NA results.
head(top_edges(res, coef = "conditiontreated"))
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>   motif motif_type logFC PValue FDR model_used
#> 1   N_A       node    NA     NA  NA     volume
#> 2   N_B       node    NA     NA  NA     volume
#> 3 E_A_A       edge    NA     NA  NA     volume
#> 4 E_A_B       edge    NA     NA  NA     volume
#> 5 E_B_B       edge    NA     NA  NA     volume
head(top_triplets(res, strategy = "topology", coef = "conditiontreated"))
#>           motif motif_type logFC PValue FDR model_used
#> T_A_A_B T_A_A_B   triangle    NA     NA  NA   topology
#> T_A_B_B T_A_B_B   triangle    NA     NA  NA   topology
```
