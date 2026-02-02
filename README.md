# CellEdgeR
Differential analysis of cell-type interaction motifs from spatial coordinates.

## Overview
CellEdgeR builds per-sample Delaunay graphs, counts motifs (nodes, edges, triangles, optional wedges), and fits edgeR models.
Offsets are computed once and stored in the cellgraph (offset modes: `volume`, `hier_null`):
- Node offsets: `log(total cells) + log(TMM)` when available.
- Edge/triangle/wedge offsets: structural Chung-Lu (volume-based).
Motif keys are prefixed by layer: `N_` (nodes), `E_` (edges), `T_` (triangles), `W_` (wedge), and `TP_` when triangles/wedges are merged into unordered triplets.

## Installation
Install dependencies and then install the package:

```r
install.packages(c("geometry", "Matrix"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")

devtools::install_github("GeNeHetX/CellEdgeR")
library(CellEdgeR)
```

## Typical workflow

1. **Prepare your data**
   Provide a *named* list of samples. Each element is a data frame with numeric `x`, `y` coordinates in the first two
   columns and a label vector (cell type, cluster, etc.) in the third column. Optionally add a logical
   `boundary`/`is_boundary` column to flag boundary cells for erosion (excluded by default).

   ```r
   samples_list <- make_demo_samples(n_cells=10000,seed = 42) # enriches T_a_b_c in the second half by default

   metadata_df <- data.frame(
     condition = rep(c("ctrl", "treated"), each = length(samples_list) / 2),
     batch = sample(c("batch1", "batch2"), length(samples_list), replace = TRUE),
     row.names = names(samples_list)
   )

   ```

2. **Build graphs, then count motifs (same cellgraph object)**
   ```r
   cellgraph <- build_cell_graphs(samples_list, n_cores = 4, verbose = TRUE)
   cellgraph <- count_motifs_graphs(
     cellgraph,
     max_edge_len = 50,     # set NA/NULL/<=0 to disable pruning
     include_wedge = TRUE,
     n_cores = 4
   )
   ```

   `build_cell_graphs()` returns a `cellEdgeR_obj` with per-sample Delaunay graphs. `count_motifs_graphs()` reuses that
   object and adds motif counts, exposures, offsets, normalized counts, relative counts, and run parameters. Offsets,
   normalized counts, and relative counts are stored per offset mode (e.g., `cellgraph$norm_counts$volume`,
   `cellgraph$norm_counts$hier_null`).

3. **Fit differential models (offsets already computed)**
   ```r
   cellgraph <- motif_edger(
     cellgraph = cellgraph,
     sample_df = metadata_df,  # rownames must match sample names
     design_formula = "~ condition + batch"
   )

   simple_tbl <- top_motifs_simple(cellgraph, coef = "conditiontreated")
   triplet_tbl <- top_motifs_triplet(cellgraph, strategy = "ancova", coef = "conditiontreated")
   head(simple_tbl)
   head(triplet_tbl)
   ```

   `motif_edger()` runs volume and ancova strategies by default and stores results in
   `cellgraph$edger$strategies`. Use `top_motifs_simple()` for node/edge motifs (volume offsets)
   and `top_motifs_triplet()` for 3-node motifs (volume or ancova). For triplets, use
   `triplet_mode = "merge"` to collapse wedges+triangles into unordered triplet motifs, or
   `triplet_mode = "closure"` to model wedges separately and test triangle closure using total
   triples (open+closed) as a covariate. Both require `count_motifs_graphs(..., include_wedge = TRUE)`.

4. **Inspect/visualize motifs (optional)**

   ```r
   norm_df <- get_motif_values(cellgraph, value = "norm")
   head(norm_df)
   ```

   ```r
   library(ggplot2)

   motif_key <- "T_a_b_c"
   df_plot <- data.frame(sample = rownames(norm_df), value = norm_df[[motif_key]])
   df_plot <- merge(df_plot, metadata_df, by.x = "sample", by.y = "row.names", all.x = TRUE)
   ggplot(df_plot, aes(x = condition, y = value, fill = condition)) +
     geom_boxplot() +
     geom_jitter(width = 0.1, alpha = 0.6) +
     ylab(paste0("Normalized (log2) count: ", motif_key)) +
     xlab("condition")
   ```

   ```r
   plot_sample_graph(
     cellgraph,
     sample_name = "s1",
     max_edge_len = 50,
     motif_key = "E_a_b",
     motif_node_size=0.5,
     dim_node_nonmotif = 0.1,
     alpha_node_nonmotif = 0.2,
     alpha_edge_nonmotif = 0.2
   )
   ```

## Choosing a pruning threshold
The `max_edge_len` threshold controls which Delaunay edges are retained before motif counting. To pick a cutoff,
inspect the edge-length distribution and drop the sparse tail while preserving biologically relevant neighbors.

```r
edge_lengths <- unlist(lapply(cellgraph$per_sample_graph, function(ps) ps$edge_len))
summary(edge_lengths)

ggplot2::ggplot(data.frame(edge_len = edge_lengths), ggplot2::aes(edge_len)) +
  ggplot2::geom_density(fill = "#4477aa", alpha = 0.3, color = "#223355") +
  ggplot2::scale_x_log10("Edge length (log scale)") +
  ggplot2::ylab("Density") +
  ggplot2::ggtitle("Edge length distribution (density, log-x)") +
  ggplot2::theme_minimal()
```

## Tips
- Reuse the `graph_obj` when experimenting with pruning or the wedge option (`include_wedge`) to avoid rerunning the
  triangulation.
- If you want erosion, provide a boundary mask per sample (column `boundary`/`is_boundary` or `erosion_cells`); motifs
  touching those cells are excluded without re-triangulating.
- Record the `max_edge_len` used for each run so downstream results stay traceable.
- For quick inspection, use `get_motif_values(cellgraph, value = "norm")` and apply your own plotting/statistics.
- `get_motif_values(..., value = "norm")` uses the `volume` offsets only, matching the hybrid modeling offsets.
- `motif_space_size()` reports the combinatorial number of possible motif labels given the stored label set and triplet mode.

## Glossary
- `make_demo_samples()`: Create a synthetic named list of samples for quick testing (defaults to a boosted `T_a_b_c` motif in the second half; set `boost_motif = FALSE` to disable).
- `build_cell_graphs()`: Build per-sample Delaunay graphs from coordinates and labels.
- `count_motifs_graphs()`: Count node/edge/triangle (and optional wedge) motifs and compute offsets/normalization.
- `merge_motif_objs()`: Merge two motif objects and recompute offsets/normalized counts.
- `motif_edger()`: Fit differential motif models (volume, ancova) and store results.
- `top_motifs_simple()`: Return ranked node/edge motifs (volume offsets).
- `top_motifs_triplet()`: Return ranked 3-node motifs (volume or ancova; separate or merged).
- `motif_space_size()`: Count possible motif labels based on label set and triplet mode.
- `get_motif_values()`: Return raw/normalized motif values (optionally including lower-order submotifs; normalized uses `volume` offsets).
- `plot_sample_graph()`: Plot one sample's graph and highlight motifs.
- `cellEdgeR_obj`: Main object storing `sample_name`, counts, exposures, offsets, and edgeR results.
- `cellEdgeR_graphs`: Graph-only object with `per_sample_graph` ready for motif counting.
