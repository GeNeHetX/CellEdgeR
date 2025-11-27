\name{build_cell_graphs}
\alias{build_cell_graphs}
\title{Build cached cell graphs for motif counting}
\description{
Builds a constrained Delaunay triangulation for each sample and caches the raw edges so multiple
motif-counting passes can reuse the same triangulation.
}
\usage{
build_cell_graphs(cells_by_sample, n_cores = 1, verbose = TRUE)
}
\arguments{
  \item{cells_by_sample}{Named list of data frames where the first two columns are numeric \code{x, y} coordinates
    and the third column holds the cell-type label. Names of the list are treated as sample identifiers.}
  \item{n_cores}{Integer number of cores to request; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.}
  \item{verbose}{Logical; print progress messages for standardization and graph construction.}
}
\details{
This helper performs the expensive Delaunay triangulation step and records every edge's length.
Pass the returned \code{cellEdgeR_graphs} object to \code{count_motifs_graphs(graph_obj = ...)} to compute
motif counts without re-running the triangulation, even when the pruning or wedge options change.
}
\value{
A list of class \code{cellEdgeR_graphs} containing:
\describe{
  \item{\code{samples}}{Character vector of sample IDs.}
  \item{\code{label_levels}}{Sorted unique cell-type labels spanning all samples.}
  \item{\code{lab_to_id}}{Global label-to-integer lookup.}
  \item{\code{per_sample}}{List with entries for each sample that store \code{edges}, their lengths \code{edge_len}, and the label IDs/characters.}
}}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
}
