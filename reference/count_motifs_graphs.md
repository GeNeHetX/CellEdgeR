# Count cell-type motifs on spatial graphs

Reuses a constrained Delaunay triangulation, prunes edges that exceed
`max\_edge\_len`, and counts cell-type singletons, unordered pairs
(edges), and unordered triplets (triangles). Triangle counts are
accelerated via the C++ helper exposed in CellEdgeR.

## Usage

``` r
count_motifs_graphs(graph_obj, max_edge_len = NA_real_, include_wedge = FALSE, verbose = TRUE,
  offset_pseudo = 1, n_cores = 1)
```

## Arguments

- graph_obj:

  Output of
  [`build_cell_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/build_cell_graphs.md);
  the triangulation is reused and only pruning/wedge parameters are
  applied here.

- max_edge_len:

  Numeric edge-length threshold; Delaunay edges longer than this are
  dropped before counting. Set to `NA`, `NULL`, or `<= 0` to keep the
  full Delaunay graph.

- include_wedge:

  Logical flag; when `TRUE` also counts open three-node paths (wedge)
  and records them under `counts\$wedge` plus `exposure\$wedge`.

- verbose:

  Logical; print progress messages for graph construction, counting
  stages, and (when enabled) wedge collection.

- offset_pseudo:

  Small positive constant used inside offsets.

- n_cores:

  Parallelism hint for motif counting; values greater than 1 trigger
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix-alikes, while 1 runs sequentially.

## Details

All samples share a global label vocabulary so that counts align across
matrices. Call
[`build_cell_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/build_cell_graphs.md)
to create the reusable triangulation, then reuse it across different
pruning/wedge configurations. Counts are sparse Matrix objects (node for
cells, edge for pairs, triangle for triplets, plus wedge when
requested). Enabling wedge adds `counts\$wedge` and `exposure\$wedge`;
the latter also reports wedge totals alongside edge/triangle exposures.

## Value

A named list with these entries:

- `sample_name`:

  Character vector of sample names (names of `cells_by_sample`).

- `label_levels/lab_to_id`:

  Sorted unique cell-type labels and integer lookup.

- `per_sample_graph`:

  Cached per-sample graph structures (edges, lengths, labels,
  coordinates).

- `raw_count`:

  Sparse matrices for `node`, `edge`, `triangle`, and optional `wedge`
  motifs (rows = motifs, cols = samples).

- `exposure`:

  Totals used for offsets: `cells`, `edges`, `triangles`, `volumes` per
  label×sample, and `wedge` when requested.

- `offsets`:

  Offset sets (e.g., `volume`, `hier_null`), each containing
  log-expected matrices per layer. Node offsets are adjusted with TMM
  factors; other layers use the structural offsets only.

- `norm_counts`:

  Normalized counts per offset set (counts / `exp(offset)`).

- `relative_counts`:

  edgeR intercept-only log2 residuals per offset set.

- `offset_part_id`:

  For each offset set, the components contributing to each motif's
  offset.

- `offset_part_values`:

  Numeric values referenced by `offset_part_id` (e.g., volumes, 2m, edge
  posteriors).

- `edger`:

  edgeR fits/tests once
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
  is run.

- `parameters`:

  Run parameters (e.g., `max_edge_len`, `include_wedge`,
  `offset_pseudo`, available `offset_modes`, layer names,
  node_tmm_offsets flag).

## Examples

``` r
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
str(counts_full$raw_count)
#> List of 4
#>  $ node    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:4] 0 1 0 1
#>   .. ..@ p       : int [1:3] 0 2 4
#>   .. ..@ Dim     : int [1:2] 2 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:2] "N_A" "N_B"
#>   .. .. ..$ : chr [1:2] "s1" "s2"
#>   .. ..@ x       : num [1:4] 2 1 1 2
#>   .. ..@ factors : list()
#>  $ edge    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:4] 0 1 1 2
#>   .. ..@ p       : int [1:3] 0 2 4
#>   .. ..@ Dim     : int [1:2] 3 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:3] "E_A_A" "E_A_B" "E_B_B"
#>   .. .. ..$ : chr [1:2] "s1" "s2"
#>   .. ..@ x       : num [1:4] 1 2 2 1
#>   .. ..@ factors : list()
#>  $ triangle:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:2] 0 1
#>   .. ..@ p       : int [1:3] 0 1 2
#>   .. ..@ Dim     : int [1:2] 2 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:2] "T_A_A_B" "T_A_B_B"
#>   .. .. ..$ : chr [1:2] "s1" "s2"
#>   .. ..@ x       : num [1:2] 1 1
#>   .. ..@ factors : list()
#>  $ wedge   :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int(0) 
#>   .. ..@ p       : int [1:3] 0 0 0
#>   .. ..@ Dim     : int [1:2] 0 2
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:2] "s1" "s2"
#>   .. ..@ x       : num(0) 
#>   .. ..@ factors : list()
```
