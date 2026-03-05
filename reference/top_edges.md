# Top motif results

Convenience helpers for retrieving ranked motif results from
[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md).
`top_edges()` returns node/edge motifs. `top_triplets()` returns 3-node
motifs using either co-occurrence (merged triplets, volume offsets) or
topology (submotif-adjusted, closure mode).

## Usage

``` r
top_edges(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_triplets(cellgraph, strategy = c("cooccurrence", "topology"),
  coef = NULL, model = c("full", "null"), n = Inf, fdr_method = "BH")
```

## Details

**Co-occurrence (strategy = `cooccurrence`):** wedges and triangles are
merged into unordered triplets (`TP_*`) and modeled with volume offsets
only. This tests whether a label triplet is over/under represented as a
whole, regardless of closure. Computed by default when
[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
is run with `include_wedge = TRUE`; explicitly use
`triplet_mode = "merge"` to recompute only this mode.

**Topology (strategy = `topology`):** wedges and triangles are modeled
separately in closure mode. Wedges are *submotif-adjusted* using
edge-derived covariates (edge-force), and triangles are tested for
closure using total triples as a covariate (submotif-adjusted). Computed
by default when
[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
is run with `include_wedge = TRUE`; explicitly use
`triplet_mode = "closure", strategies = "submotif_adj"` to recompute
only this mode.

**Model selection (full vs null):** when `strategy = "cooccurrence"`,
`model = "full"` uses the full design formula, and `model = "null"` uses
an intercept-only model (useful for diagnostics or when you want the
baseline rate). Topology uses the submotif-adjusted model and ignores
`model`.

## Arguments

- cellgraph:

  A `cellEdgeR_obj` with `edger` results.

- coef:

  Coefficient name or index; defaults to the first non-intercept
  coefficient.

- model:

  Which stored model to use for edgeR strategies: `full` or `null`.
  `full` uses the provided design formula; `null` uses an intercept-only
  model. Ignored when `strategy = "topology"`.

- n:

  Number of motifs to return; defaults to all.

- fdr_method:

  Multiple testing correction method for `p.adjust` (default `BH`).

- strategy:

  Triplet strategy: `cooccurrence` (volume offsets, merged triplets) or
  `topology` (submotif-adjusted, closure mode).

## Value

A data frame with columns: motif, motif\\type, logFC, PValue, FDR, and
model\\used.

## See also

[`motif_edger`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md),
[`motif_space_size`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_space_size.md)

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

res_merge <- motif_edger(motifs, sample_df, "~ condition", triplet_mode = "merge")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> Fitting submotif-adjusted models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Submotif-adjusted dispersion could not be estimated; returning NA results.
res_close <- motif_edger(motifs, sample_df, "~ condition", triplet_mode = "closure")
#> Fitting edgeR (QL) for volume offsets...
#> Warning: No residual df: cannot estimate dispersion
#> edgeR dispersion estimation failed for volume/full model; tests will be empty.
#> Fitting submotif-adjusted models (per motif)...
#> Warning: No residual df: cannot estimate dispersion
#> Warning: Submotif-adjusted dispersion could not be estimated; returning NA results.

top_edges(res_merge, coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>   motif motif_type logFC PValue FDR model_used
#> 1   N_A       node    NA     NA  NA     volume
#> 2   N_B       node    NA     NA  NA     volume
#> 3 E_A_A       edge    NA     NA  NA     volume
#> 4 E_A_B       edge    NA     NA  NA     volume
#> 5 E_B_B       edge    NA     NA  NA     volume
top_triplets(res_merge, strategy = "cooccurrence", coef = "conditiontreated")
#> Warning: No edgeR tests stored for volume coef: conditiontreated. Returning NA results.
#>      motif motif_type logFC PValue FDR   model_used
#> 6 TP_A_A_B    triplet    NA     NA  NA cooccurrence
#> 7 TP_A_B_B    triplet    NA     NA  NA cooccurrence
top_triplets(res_close, strategy = "topology", coef = "conditiontreated")
#>           motif motif_type logFC PValue FDR model_used
#> T_A_A_B T_A_A_B   triangle    NA     NA  NA   topology
#> T_A_B_B T_A_B_B   triangle    NA     NA  NA   topology
```
