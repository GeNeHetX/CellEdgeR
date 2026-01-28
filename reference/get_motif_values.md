# Retrieve motif values (raw, normalized)

Pulls raw counts or normalized counts for a motif (or its lower-order
submotifs) into a samples-by-motifs data frame.

## Usage

``` r
get_motif_values(cellgraph, motif_key = NULL, include_submotifs = FALSE,
  value = c("raw", "norm"))
```

## Arguments

- cellgraph:

  A `cellEdgeR_obj` returned by
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md).

- motif_key:

  Character vector of motif keys (e.g., `E_A_B`, `T_A_B_C`). If `NULL`,
  return all motifs.

- include_submotifs:

  Logical; if `TRUE`, include lower-order submotifs implied by the motif
  key (edges for triangles/wedges, nodes for edges).

- value:

  Which value to return: `raw` or `norm`.

## Value

A data frame with sample names as row names and one column per motif.
When `include_submotifs = TRUE`, requested motifs appear first, followed
by unique submotifs.

## Details

Normalized values are computed directly from the volume offsets so they
match the offsets used by
[`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
(including node TMM adjustments when present). Call `get_motif_values()`
multiple times to retrieve different value types. Submotifs are derived
structurally (edges for triangles/wedges, nodes for edges) and then
filtered to motifs present in the object.

## See also

[`count_motifs_graphs`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md),
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
tri_key <- rownames(motifs$raw_count$triangle)[1]
get_motif_values(motifs, motif_key = tri_key, value = "raw")
#>    T_A_A_B
#> s1       1
#> s2       0
get_motif_values(motifs, motif_key = tri_key, value = "norm")
#>    T_A_A_B
#> s1   1.125
#> s2   0.000
```
