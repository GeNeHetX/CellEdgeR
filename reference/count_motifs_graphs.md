# Count cell-type motifs on Delaunay graphs

Reuse a prebuilt Delaunay triangulation, prune edges that exceed
`max_edge_len`, and count cell-type singletons, unordered pairs (edges),
unordered triplets (triangles), and optionally wedge. Triangle counts
are accelerated via the C++ helper exposed in CellEdgeR.

## Usage

``` r
count_motifs_graphs(graph_obj, max_edge_len = NA, include_wedge = FALSE,
  verbose = TRUE, offset_pseudo = 1, n_cores = 1, erosion = TRUE,
  erosion_cells = NULL)
```

## Arguments

- graph_obj:

  Output of
  [`build_cell_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/build_cell_graphs.md);
  the triangulation is reused and only pruning/wedge parameters are
  applied here.

- max_edge_len:

  Numeric threshold; Delaunay edges longer than this are dropped. Set to
  `NA`, `NULL`, or `<= 0` to skip pruning.

- include_wedge:

  Logical; if `TRUE`, returns open-triplet (wedge) counts alongside
  triangles.

- verbose:

  Logical; print progress.

- offset_pseudo:

  Small positive constant used inside offsets.

- n_cores:

  Parallelism hint for motif counting; values greater than 1 trigger
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix-alikes, while 1 runs sequentially.

- erosion:

  Logical; if `TRUE` (default), exclude motifs that touch boundary
  cells. Boundary cells can be provided as a logical column (e.g.
  `boundary`/`is_boundary`) in the input data or via `erosion_cells`.

- erosion_cells:

  Optional list of boundary masks/indices by sample name. Each entry can
  be a logical vector (length = number of cells) or integer indices to
  exclude.

## Value

A cellgraph object (class `cellEdgeR_obj`) with the original graph info
plus:

- `raw_count`:

  Sparse matrices for `node`, `edge`, `triangle`, and optional `wedge`
  motifs (rows = motifs, cols = samples).

- `exposure`:

  Totals used in offsets: `cells`, `edges`, `triangles`, `volumes` per
  label by sample, `center_pairs`, plus `wedge` and `triples` when
  requested.

- `offsets`:

  Volume offsets containing log-expected matrices per layer. Node
  offsets are retained for normalization and submotif inspection; tested
  layers use structural volume offsets.

- `norm_counts`:

  Normalized counts for the volume offset (counts / `exp(offset)`).

- `relative_counts`:

  edgeR intercept-only log2 residuals for the volume offset.

- `offset_part_id`:

  Components contributing to each motif's volume offset.

- `offset_part_values`:

  Numeric values referenced by `offset_part_id` (e.g. volumes and `2m`).

- `edger`:

  edgeR fits/tests once
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
  is run.

- `parameters`:

  Run parameters including `max_edge_len`, `include_wedge`,
  `offset_pseudo`, `offset_modes`, `offset_version`, layer names, and
  `node_tmm_offsets`.

## Examples

``` r
demo <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(demo, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
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
