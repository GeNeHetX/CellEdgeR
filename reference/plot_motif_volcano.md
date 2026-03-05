# Volcano plot for motif differential results

Plot differential motif results (from
[`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
or
[`top_triplets()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md))
in a standard volcano layout, with flexible motif labeling options.

## Usage

``` r
plot_motif_volcano(
  x,
  result_type = c("edges", "triplets_cooccurrence", "triplets_topology"),
  coef = NULL,
  model = c("full", "null"),
  fdr_method = "BH",
  fdr_cutoff = 0.05,
  logFC_cutoff = 0,
  label_mode = c("none", "manual", "top_n", "top_bottom_n"),
  label_motifs = NULL,
  label_n = 10,
  label_top_by = c("PValue", "FDR", "abs_logFC"),
  point_size = 1.8,
  point_alpha = 0.8,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- x:

  A results data frame (e.g., output of
  [`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
  or
  [`top_triplets()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)),
  or a `cellEdgeR_obj` containing
  [`motif_edger()`](https://GeNeHetX.github.io/CellEdgeR/reference/motif_edger.md)
  results.

- result_type:

  Which results to fetch when `x` is a `cellEdgeR_obj`: `"edges"`,
  `"triplets_cooccurrence"`, or `"triplets_topology"`.

- coef:

  Coefficient name or index passed to
  [`top_edges()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
  /
  [`top_triplets()`](https://GeNeHetX.github.io/CellEdgeR/reference/top_edges.md)
  when `x` is a `cellEdgeR_obj`.

- model:

  Model (`"full"` or `"null"`) used for `"edges"` and
  `"triplets_cooccurrence"` when `x` is a `cellEdgeR_obj`.

- fdr_method:

  Multiple testing correction method passed to `p.adjust` when needed.

- fdr_cutoff:

  FDR threshold used for point coloring.

- logFC_cutoff:

  Absolute logFC threshold used for point coloring.

- label_mode:

  Label selection mode: `"none"`, `"manual"`, `"top_n"`,
  `"top_bottom_n"`.

- label_motifs:

  Motif keys to label when `label_mode = "manual"`.

- label_n:

  Number of motifs to label for `"top_n"` and `"top_bottom_n"`.

- label_top_by:

  Ranking criterion for `"top_n"` labels: `"PValue"`, `"FDR"`,
  `"abs_logFC"`.

- point_size:

  Point size in the volcano.

- point_alpha:

  Point alpha in the volcano.

- title:

  Optional plot title.

- subtitle:

  Optional plot subtitle.

## Value

A `ggplot` object.
