# Triptych volcano with ordered motif side panels

Create a 3-panel figure with the most depleted motifs on the left, the
volcano in the middle, and the most enriched motifs on the right. Side
motifs are displayed in fixed top-to-bottom order.

## Usage

``` r
plot_motif_volcano_triptych(
  x,
  result_type = c("edges", "triplets_cooccurrence", "triplets_topology"),
  coef = NULL,
  model = c("full", "null"),
  fdr_method = "BH",
  n_side = 5,
  left_motifs = NULL,
  right_motifs = NULL,
  label_side_motifs = TRUE,
  side_text_mode = c("numbered", "verbose"),
  side_label_truncate = 3,
  label_mode = c("none", "manual", "top_n", "top_bottom_n"),
  label_motifs = NULL,
  label_n = 6,
  label_top_by = c("PValue", "FDR", "abs_logFC"),
  fdr_cutoff = 0.05,
  logFC_cutoff = 0,
  point_size = 1.8,
  point_alpha = 0.8,
  panel_widths = c(1.35, 2.4, 1.35),
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

- n_side:

  Number of automatically selected motifs on each side.

- left_motifs:

  Optional motif keys to force on the depleted (left) panel.

- right_motifs:

  Optional motif keys to force on the enriched (right) panel.

- label_side_motifs:

  Logical; if `TRUE`, motifs shown on side panels are labeled on the
  volcano.

- side_text_mode:

  Side-panel annotation mode. `"numbered"` (default) uses compact
  numeric IDs shared with volcano labels; `"verbose"` uses the previous
  motif/logFC/FDR text.

- side_label_truncate:

  Number of characters to keep for cell labels inside side-panel motifs.
  Set to `NA` for no truncation.

- label_mode:

  Additional volcano label mode: `"none"`, `"manual"`, `"top_n"`,
  `"top_bottom_n"`.

- label_motifs:

  Motif keys to label when `label_mode = "manual"`.

- label_n:

  Number of motifs to label for `"top_n"` and `"top_bottom_n"`.

- label_top_by:

  Ranking criterion for `"top_n"` labels: `"PValue"`, `"FDR"`,
  `"abs_logFC"`.

- fdr_cutoff:

  FDR threshold used for volcano point coloring.

- logFC_cutoff:

  Absolute logFC threshold used for volcano point coloring.

- point_size:

  Point size in the volcano.

- point_alpha:

  Point alpha in the volcano.

- panel_widths:

  Relative widths for left/center/right panels.

- title:

  Optional volcano title.

- subtitle:

  Optional volcano subtitle.

## Value

If package `patchwork` is installed, returns a combined patchwork plot.
Otherwise draws the 3 panels and returns a list with components:
`depleted_panel`, `volcano`, `enriched_panel`.
