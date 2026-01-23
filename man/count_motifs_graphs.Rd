\name{count_motifs_graphs}
\alias{count_motifs_graphs}
\title{Count cell-type motifs on spatial graphs}
\description{
Reuses a constrained Delaunay triangulation, prunes edges that exceed \code{max\_edge\_len},
and counts cell-type singletons, unordered pairs (edges), and unordered triplets (triangles).
Triangle counts are accelerated via the C++ helper exposed in \pkg{CellEdgeR}.
}
\usage{
count_motifs_graphs(graph_obj, max_edge_len = NA_real_, include_wedge = FALSE, verbose = TRUE,
  offset_pseudo = 1, n_cores = 1)
}
\arguments{
  \item{graph_obj}{Output of \code{build_cell_graphs()}; the triangulation is reused and only pruning/wedge parameters are applied here.}
  \item{max_edge_len}{Numeric edge-length threshold; Delaunay edges longer than this are dropped before counting.
    Set to \code{NA}, \code{NULL}, or \code{<= 0} to keep the full Delaunay graph.}
  \item{include_wedge}{Logical flag; when \code{TRUE} also counts open three-node paths (wedge) and records them under \code{counts\$wedge} plus \code{exposure\$wedge}.}
  \item{verbose}{Logical; print progress messages for graph construction, counting stages, and (when enabled) wedge collection.}
  \item{offset_pseudo}{Small positive constant used inside offsets.}
  \item{n_cores}{Parallelism hint for motif counting; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.}
}
\details{
All samples share a global label vocabulary so that counts align across matrices.
Call \code{build_cell_graphs()} to create the reusable triangulation, then reuse it across different pruning/wedge configurations.
Counts are sparse \pkg{Matrix} objects (node for cells, edge for pairs, triangle for triplets, plus wedge when requested).
Enabling wedge adds \code{counts\$wedge} and \code{exposure\$wedge}; the latter also reports wedge totals
alongside edge/triangle exposures.
}
\value{
A named list with these entries:
\describe{
  \item{\code{sample_name}}{Character vector of sample names (names of \code{cells_by_sample}).}
  \item{\code{label_levels/lab_to_id}}{Sorted unique cell-type labels and integer lookup.}
  \item{\code{per_sample_graph}}{Cached per-sample graph structures (edges, lengths, labels, coordinates).}
  \item{\code{raw_count}}{Sparse matrices for \code{node}, \code{edge}, \code{triangle}, and optional \code{wedge} motifs (rows = motifs, cols = samples).}
  \item{\code{exposure}}{Totals used for offsets: \code{cells}, \code{edges}, \code{triangles}, \code{volumes} per label×sample, and \code{wedge} when requested.}
  \item{\code{offsets}}{Offset sets (e.g., \code{volume}, \code{hier_null}), each containing log-expected matrices per layer.
  Node offsets are adjusted with TMM factors; other layers use the structural offsets only.}
  \item{\code{norm_counts}}{Normalized counts per offset set (counts / \code{exp(offset)}).}
  \item{\code{relative_counts}}{edgeR intercept-only log2 residuals per offset set.}
  \item{\code{offset_part_id}}{For each offset set, the components contributing to each motif's offset.}
  \item{\code{offset_part_values}}{Numeric values referenced by \code{offset_part_id} (e.g., volumes, 2m, edge posteriors).}
  \item{\code{edger}}{edgeR fits/tests once \code{motif_edger()} is run.}
  \item{\code{parameters}}{Run parameters (e.g., \code{max_edge_len}, \code{include_wedge}, \code{offset_pseudo}, available \code{offset_modes}, layer names, node_tmm_offsets flag).}
}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
str(counts_full$raw_count)
}
