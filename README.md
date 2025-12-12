# CellEdgeR
Differential analysis of cell-type interaction motifs built from spatial coordinates.

## Overview
`CellEdgeR` builds spatial graphs via Delaunay triangulation, counts cell-type motifs (single cells, edges, triangles, and optional wedges), and finally fits edgeR models to test for differential motif enrichment. The workflow separates graph construction from the motif counting/fitting steps so you can experiment with pruning thresholds without repeatedly triangulating the raw samples.

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

## Typical workflow
1. **Prepare your data**  
   Provide a *named* list of samples where each element is a data frame with numeric `x`, `y` coordinates in the first two columns and a label vector (cell type, cluster, etc.) in the third column.

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

2. **Build the intact graphs**  
   ```r
   graphs <- build_cell_graphs(samples_list, n_cores = 4, verbose = TRUE)
   ```  
   This step can take some time for large datasets; the triangulations are stored in the returned object so you can repeatedly prune edges or recompute motifs without rebuilding. Keep `verbose = TRUE` to monitor progress.

3. **Count motifs with optional pruning (choose `max_edge_len`)**  
   Before pruning, inspect edge lengths from `build_cell_graphs()` to identify a cutoff that removes the long, likely spurious edges while keeping the bulk of neighbors. Example analysis is shown in the next section.
   ```r
   motifs <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 50, include_wedges = TRUE)
   ```  
   Set `max_edge_len` to `NA`/`NULL`/`<= 0` to skip pruning. Use the stored `graph_obj` to avoid rerunning the expensive triangulation.

4. **Fit differential models**
   ```r
   results <- motif_edger(
     motif_obj = motifs,
     sample_df = metadata_df,  # data frame with sample metadata (rownames = sample IDs)
     design_formula = "~ condition + batch",
     coef = "conditiontreated"
   )
   ```  
   Provide sample metadata (rownames matching the sample IDs) and specify the model terms you want to test. Choose `fdr_method = "dagger"` only if you require the hierarchical DAG correction; otherwise the default BH correction is applied per layer.

5. **Visualize motifs (optional)**  
   <!-- - Boxplot of normalized counts for a specific motif across samples (grouped by condition):  
     ```r
     norm <- normalize_motif_counts(motifs, pseudo = 0.5)
     plot_motif_box(norm, motif_key = "a_a_b", layer = "size3", sample_df = metadata_df, group_var = "condition")
     ```
     Or switch to likelihood ratios instead of normalized counts:
     ```r
     lr <- normalize_motif_counts(motifs, pseudo = 0.5, return_lr = TRUE)
     plot_motif_box(lr, motif_key = "a_a_b", layer = "size3", sample_df = metadata_df, group_var = "condition")
     ```
   Or use a custom plot -->
   Use a custom plot
   ```r
   library(ggpubr)
      
   normcnt <- edgeRnorm_motif_counts(motifs)
   mergedf=data.frame((t(as.matrix(normcnt$size3))[rownames(metadata_df),]),metadata_df)
   ggboxplot(mergedf,x="batch",fill="condition",y="a_a_b")

   ```

   - Plot a sample graph and highlight a motif (nodes and edges). If your `graph_obj` was built with an older CellEdgeR that did not store coordinates, also pass the original `cells_by_sample` list:  
     ```r
     plot_sample_graph(graphs, sample_id = "s1", max_edge_len = 50,
                       motif_key = "a_a_b", motif_layer = "size3",
                       cells_by_sample = samples_list)


     ```
     - Or plot multiple ones:
     ```r
     library(ggpubr)
      ggarrange(plotlist=lapply(paste0("s",c(1:3,22:24)),plot_sample_graph,graph_obj=graphs, max_edge_len = 50,motif_key = "a_a_b", motif_layer = "size3"))
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
