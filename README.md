# CellEdgeR
Differential analysis of cell-type interaction motifs built from spatial coordinates.

## Overview
`CellEdgeR` builds spatial graphs via Delaunay triangulation, counts cell-type motifs (nodes, edges, triangles, optional wedge), and fits edgeR models. Offsets are precomputed once per run and stored in the same cellgraph. Motif names are prefixed by layer: `N_` (nodes), `E_` (edges), `T_` (triangles), `W_` (wedge), and `TW_` when triangles/wedge are merged for modeling.

## Installation
Add the package dependencies before installing `CellEdgeR`:

```r
install.packages(c("geometry", "Matrix"))
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
```

Then install and load `CellEdgeR` as you would any other R package:

```r
devtools::install_github("GeNeHetX/CellEdgeR")
library(CellEdgeR)
```

## Typical workflow (two simple steps)
1. **Prepare your data**  
   Provide a *named* list of samples where each element is a data frame with numeric `x`, `y` coordinates in the first two columns and a label vector (cell type, cluster, etc.) in the third column. Inputs are validated (numeric coords, non-NA labels, reasonable number of labels) with clear errors.

   To generate a quick synthetic dataset for testing, call the helper bundled with the package:

   ```r
   library(CellEdgeR)
   samples_list <- make_demo_samples(seed = 42)
   head(samples_list[[1]])

   metadata_df <- data.frame(
     condition = rep(c("ctrl", "treated"), each = length(samples_list) / 2),
     batch = sample(c("batch1", "batch2"),length(samples_list),replace=T),
     row.names = names(samples_list)
   )
   ```

2. **Build graphs, then count motifs (same cellgraph object)**  
   ```r
   cellgraph <- build_cell_graphs(samples_list, n_cores = 4, verbose = TRUE)
   cellgraph <- count_motifs_graphs(cellgraph,      # graph must be built first
     max_edge_len = 50,          # set NA/NULL/<=0 to disable pruning
     include_wedge = TRUE,
     n_cores = 4                 # optional parallel counting
   )
   ```
   `build_cell_graphs()` returns a cellgraph (`cellEdgeR_obj`) with per-sample Delaunay graphs. `count_motifs_graphs()` reuses that object and adds motif `counts`, `exposure` (edges, triangles, volumes), `offsets`, `norm_counts` (offset ratios), `relative_counts` (edgeR intercept-only log2 residuals), and records the parameters used. Offsets, normalized counts, and residuals are stored per offset mode (e.g., `cellgraph$norm_counts$volume`, `cellgraph$norm_counts$hier_null`). Node offsets include TMM factors computed from node counts (log(total cells) + log(TMM)), while edge/triangle/wedge offsets remain structural. Motif names carry prefixes (`N_`, `E_`, `T_`, `W_`, optionally `TW_`), so numeric labels won’t be mangled.

3. **Fit differential models (offsets already computed)**
   ```r
   cellgraph <- motif_edger(
     cellgraph = cellgraph,
     sample_df = metadata_df,  # rownames must match sample names
     design_formula = "~ condition + batch",
     merge_triplets = FALSE         # set TRUE to merge triangles+wedge as TW_*
  )
  diff_tbl <- top_motifs(cellgraph, strategy = "volume", coef = "conditiontreated")
  head(diff_tbl)
  ```
  `motif_edger()` runs multiple strategies by default (`volume`, `hierarchical`, `ancova`) and stores results in `cellgraph$edger$strategies`. `top_motifs()` returns a data frame ordered by p-value with logFC, FDR (BH by default), and `model_used`; it defaults to the `hybrid` strategy (volume for node/edge motifs, ancova for triangle/wedge motifs with separate FDR per group), but you can use `strategy` to select the model and `model = "null"` for the intercept-only edgeR fits.

4. **Inspect/visualize motifs (optional)**  
   - Pull all normalized motifs into one data frame (samples as rows, motifs as columns):  
     ```r
     norm_df <- get_norm_counts(cellgraph, offset_mode = "volume")
     head(norm_df)
     ```
   - Boxplot a normalized motif across samples using ggplot2:  
     ```r
     library(ggplot2)
     motif_key <- "T_A_A_B"
     df_plot <- data.frame(sample = rownames(norm_df), value = norm_df[[motif_key]])
     df_plot <- merge(df_plot, metadata_df, by.x = "sample", by.y = "row.names", all.x = TRUE)
     ggplot(df_plot, aes(x = condition, y = value, fill = condition)) +
       geom_boxplot() +
       geom_jitter(width = 0.1, alpha = 0.6) +
       ylab(paste0("Normalized (log2) count: ", motif_key)) +
       xlab("condition")
     ```
   - Plot a sample graph and highlight a motif (lighter/smaller non-motif nodes/edges are adjustable):  
     ```r
     plot_sample_graph(
       cellgraph, sample_name = "s1", max_edge_len = 50,
       motif_key = "T_A_A_B",
       dim_node_nonmotif = 0.5, alpha_node_nonmotif = 0.3, alpha_edge_nonmotif = 0.25
     )
     ```

## Choosing a pruning threshold
The `max_edge_len` threshold in the workflow above controls which Delaunay edges are retained before motif counting. To pick a sensible cutoff, inspect the edge-length distribution from the built graphs and pick a value that removes the sparse tails while leaving the bulk of biologically relevant neighbors intact.

Example analysis:

```r
edge_lengths <- unlist(lapply(graphs$per_sample_graph, function(ps) ps$edge_len))
summary(edge_lengths)

ggplot2::ggplot(data.frame(edge_len = edge_lengths), ggplot2::aes(edge_len)) +
  ggplot2::geom_density(fill = "#4477aa", alpha = 0.3, color = "#223355") +
  ggplot2::scale_x_log10("Edge length (log scale)") +
  ggplot2::ylab("Density") +
  ggplot2::ggtitle("Edge length distribution (density, log-x)") +
  ggplot2::theme_minimal()
  
```

Plotting the histogram (or a density/ECDF) lets you visualize where most edges fall and suggests a `max_edge_len` value that trims only the longest, likely spurious connections. Re-run `count_motifs_graphs()` with that threshold and assess how sensitive the motif counts (and downstream tests) are to the change.

## Tips
- Reuse the `graph_obj` whenever experimenting with pruning or the wedge option (`include_wedge`) to avoid rerunning the costly triangulation.
- For reproducibility, record the `max_edge_len` used when generating `motif_obj` so downstream models can be traced back to the same graph.
- If you need normalized values for quick inspection, use `get_norm_counts(cellgraph, offset_mode = "volume")` (or another stored offset mode) to pull a samples-by-motifs matrix; apply your own plotting/statistics on that output.

## Glossary
- `make_demo_samples()`: Create a synthetic named list of samples for quick testing.
- `build_cell_graphs()`: Build per-sample Delaunay graphs from coordinates and labels.
- `count_motifs_graphs()`: Count node/edge/triangle (and optional wedge) motifs and compute offsets/normalization.
- `merge_motif_objs()`: Merge two motif objects and recompute offsets/normalized counts.
- `motif_edger()`: Fit differential motif models (volume, hierarchical, ancova) and store results.
- `top_motifs()`: Return a ranked data frame of motifs for a chosen strategy/coef.
- `get_motif_values()`: Return raw/normalized motif values (optionally including hierarchical submotifs).
- `get_norm_counts()`: Return a samples-by-motifs matrix of normalized counts.
- `plot_sample_graph()`: Plot one sample's graph and highlight motifs.
- `cellEdgeR_obj`: Main object storing `sample_name`, counts, exposures, offsets, and edgeR results.
- `cellEdgeR_graphs`: Graph-only object with `per_sample_graph` ready for motif counting.
