\name{count_motifs_graphs}
\alias{count_motifs_graphs}
\title{Count cell-type motifs on spatial graphs}
\description{
Builds a constrained Delaunay triangulation for each sample, prunes edges that exceed \code{max\_edge\_len},
and counts cell-type singletons, unordered pairs (edges), and unordered triplets (triangles).
Triangle counts are accelerated via the C++ helper exposed in \pkg{CellEdgeR}.
}
\usage{
count_motifs_graphs(cells_by_sample = NULL, graph_obj = NULL, max_edge_len = NA_real_, n_cores = 1, include_wedges = FALSE, verbose = TRUE)
}
\arguments{
  \item{cells_by_sample}{Named list of data frames. Each data frame must have numeric \code{x} and \code{y} coordinates
    plus a \code{label} column indicating the cell type. Names of the list are treated as sample identifiers.
    Provide this argument when you do not already have a prebuilt graph object.}
  \item{graph_obj}{Optional output of \code{build_cell_graphs()}; the triangulation is reused and \code{cells_by_sample} is ignored when provided.}
  \item{max_edge_len}{Numeric edge-length threshold; Delaunay edges longer than this are dropped before counting.
    Set to \code{NA}, \code{NULL}, or \code{<= 0} to keep the full Delaunay graph.}
  \item{n_cores}{Integer number of cores to request when \code{cells_by_sample} is supplied; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.}
  \item{include_wedges}{Logical flag; when \code{TRUE} also counts open three-node paths (wedges) and records them under \code{counts\$wedges} plus \code{exposure\$wedges}.}
  \item{verbose}{Logical; print progress messages for graph construction, counting stages, and (when enabled) wedge collection.}
}
\details{
All samples share a global label vocabulary so that counts align across matrices.
Call \code{build_cell_graphs()} if you want to reuse the same triangulation across multiple
pruning/wedge configurations; otherwise \code{cells_by_sample} is used directly and the graphs are re-built.
Counts are sparse \pkg{Matrix} objects (size 1 for cells, size 2 for edges, size 3 for triangles, plus wedges when requested).
Enabling wedges adds \code{counts\$wedges} and \code{exposure\$wedges}; the latter also reports wedge totals
alongside edge/triangle exposures.
}
\value{
A named list with these entries:
\describe{
  \item{\code{samples}}{Character vector of sample IDs (names of \code{cells_by_sample}).}
  \item{\code{label_levels}}{Sorted unique cell-type labels spanning all samples.}
  \item{\code{counts}}{List with \code{size1}, \code{size2}, and \code{size3} sparse matrices plus \code{counts\$wedges} when requested.}
  \item{\code{exposure}}{List with \code{cells}, \code{edges}, \code{triangles}, and (with \code{include_wedges = TRUE}) \code{wedges}.}
  \item{\code{meta}}{Named list with parameters such as \code{max_edge_len} and \code{include_wedges}.}
}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedges = TRUE, verbose = FALSE)
str(counts_full$counts)
}
