# Merge motif objects

Combine two `cellEdgeR_obj` outputs from
[`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md)
into a single object with recomputed offsets and normalized counts.

## Usage

``` r
merge_motif_objs(motif_obj_a, motif_obj_b, verbose = TRUE)
```

## Arguments

- motif_obj_a,motif_obj_b:

  `cellEdgeR_obj` objects returned by
  [`count_motifs_graphs()`](https://GeNeHetX.github.io/CellEdgeR/reference/count_motifs_graphs.md).

- verbose:

  Logical; print progress while recomputing offsets.

## Details

The two objects must have disjoint sample names and overlapping label
sets; if the label sets are disjoint the merge stops with an error.
Motif key names are validated against each object's label levels before
merging. Per-sample graphs are retained only when present in both
inputs; label IDs are remapped to the merged label vocabulary.

## Value

A merged `cellEdgeR_obj` containing the combined samples, with offsets,
normalized counts, relative counts, and offset metadata recomputed on
the merged data.

## Examples

``` r
demo <- make_demo_samples(seed = 1, n_samples = 4)
graphs_full <- build_cell_graphs(demo, verbose = FALSE)
full <- count_motifs_graphs(graphs_full, max_edge_len = NA_real_, verbose = FALSE)

grp <- names(demo)
obj_a <- count_motifs_graphs(build_cell_graphs(demo[grp[1:2]], verbose = FALSE), verbose = FALSE)
obj_b <- count_motifs_graphs(build_cell_graphs(demo[grp[3:4]], verbose = FALSE), verbose = FALSE)
merged <- merge_motif_objs(obj_a, obj_b, verbose = FALSE)
```
