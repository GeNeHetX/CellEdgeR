# CellEdgeR
Differential analysis of cell-type interaction motifs built from spatial coordinates.

## Overview
`CellEdgeR` builds spatial graphs via Delaunay triangulation, counts cell-type motifs (single cells, edges, triangles, and optional wedges), and finally fits edgeR models to test for differential motif enrichment. The workflow separates graph construction from the motif counting/fitting steps so you can experiment with pruning thresholds without repeatedly triangulating the raw samples. Motif names are prefixed by layer for clarity: `N_` (nodes), `E_` (edges), `T_` (triangles), `W_` (wedges), and optionally `TW_` when triangles/wedges are merged for modeling.

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
   Provide a *named* list of samples where each element is a data frame with numeric `x`, `y` coordinates in the first two columns and a label vector (cell type, cluster, etc.) in the third column. The helper validates structure, numeric coords, non-NA labels, and caps unique labels (default 100) with clear errors if not met.

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

2. **Build graphs (validated) and count motifs in one clean flow**  
   ```r
   graphs <- build_cell_graphs(samples_list, n_cores = 4, verbose = TRUE)
   motifs <- count_motifs_graphs(
     graph_obj = graphs,
     max_edge_len = 50,       # set NA/NULL/<=0 to disable pruning
     include_wedges = TRUE,
     offset_pseudo = 0.5
   )
   ```
   The returned `motifs` object contains raw counts (`counts`), exposures, offsets (`offsets`), and normalized counts (`norm_counts`, counts divided by the offsets). Motif names carry prefixes (`N_`, `E_`, `T_`, `W_`), so you can target motifs unambiguously, even when labels are numeric.

3. **Fit differential models**
   ```r
   results <- motif_edger(
     motif_obj = motifs,
     sample_df = metadata_df,  # rownames must match sample IDs
     design_formula = "~ condition + batch",
     coef_variable = "condition",   # which variable
     coef_level    = "treated",     # which level (optional)
     merge_triplets = FALSE         # set TRUE to merge triangles+wedges as TW_*
   )
   ```
   BH correction is applied per layer by default; `fdr_method = "dagger"` adds the DAG-based filter. Offsets are reused from the `motif_obj` so testing is consistent with normalization.

4. **Visualize motifs (optional)**  
   - Boxplot a normalized motif across samples (prefix inferred):  
     ```r
     plot_motif_box(motifs$norm_counts, motif_key = "T_A_A_B", sample_df = metadata_df, group_var = "condition")
     ```
   - Plot a sample graph and highlight a motif (lighter/smaller non-motif nodes/edges are adjustable):  
     ```r
     plot_sample_graph(
       graphs, sample_id = "s1", max_edge_len = 50,
       motif_key = "T_A_A_B",
       dim_node_nonmotif = 0.5, alpha_node_nonmotif = 0.3, alpha_edge_nonmotif = 0.25
     )
     ```

## Choosing a pruning threshold
The `max_edge_len` threshold in the workflow above controls which Delaunay edges are retained before motif counting. To pick a sensible cutoff, inspect the edge-length distribution from the built graphs and pick a value that removes the sparse tails while leaving the bulk of biologically relevant neighbors intact.

Example analysis:

```r
edge_lengths <- unlist(lapply(graphs$per_sample, function(ps) ps$edge_len))
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
- Reuse the `graph_obj` whenever experimenting with pruning or the wedge option (`include_wedges`) to avoid rerunning the costly triangulation.
- For reproducibility, record the `max_edge_len` used when generating `motif_obj` so downstream models can be traced back to the same graph.
- If you need counts normalized before modeling, call `normalize_motif_counts()` on the `motif_obj` prior to inspection.
