\name{build_cell_graphs}
\alias{build_cell_graphs}
\title{Build cached cell graphs for motif counting}
\description{
Builds a constrained Delaunay triangulation for each sample and caches the raw edges so multiple
motif-counting passes can reuse the same triangulation.
}
\usage{
build_cell_graphs(cells_by_sample, n_cores = 1, verbose = TRUE, max_labels = 100)
}
\arguments{
  \item{cells_by_sample}{Named list of data frames where the first two columns are numeric \code{x, y} coordinates
    and the third column holds the cell-type label. Names of the list are treated as sample names.}
  \item{n_cores}{Integer number of cores to request; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.}
  \item{verbose}{Logical; print progress messages for standardization and graph construction.}
  \item{max_labels}{Maximum number of unique labels allowed across all samples (fails early if exceeded).}
}
\details{
This helper performs the expensive Delaunay triangulation step and records every edge's length.
Pass the returned \code{cellEdgeR_graphs} object to \code{count_motifs_graphs(graph_obj = ...)} to compute
motif counts without re-running the triangulation, even when the pruning or wedge options change.
}
\value{
A list of class \code{cellEdgeR_graphs} containing:
\describe{
  \item{\code{sample_name}}{Character vector of sample names.}
  \item{\code{label_levels/lab_to_id}}{Sorted unique cell-type labels and a global label-to-integer lookup.}
  \item{\code{per_sample_graph}}{Per-sample entries storing \code{n} (node count), \code{labels_chr}, \code{labels_id}, \code{edges}, edge lengths \code{edge_len}, and coordinates \code{xy}.}
  \item{\code{raw_count/exposure/offsets/norm_counts/relative_counts}}{Empty placeholders until \code{count_motifs_graphs()} is called.}
  \item{\code{edger}}{Empty placeholder until \code{motif_edger()} is called.}
  \item{\code{parameters}}{Build parameters (placeholders for \code{max_edge_len}, \code{offset_pseudo}, \code{built_from}).}
}}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
}
