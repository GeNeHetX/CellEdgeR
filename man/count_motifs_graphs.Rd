\name{count_motifs_graphs}
\alias{count_motifs_graphs}
\title{Count cell-type motifs on Delaunay graphs}
\description{
Reuse a prebuilt Delaunay triangulation, prune edges that exceed
\code{max_edge_len}, and count cell-type singletons, unordered pairs (edges),
unordered triplets (triangles), and optionally wedge. Triangle counts are
accelerated via the C++ helper exposed in \pkg{CellEdgeR}.
}
\usage{
count_motifs_graphs(graph_obj, max_edge_len = NA, include_wedge = FALSE,
  verbose = TRUE, offset_pseudo = 1, n_cores = 1, erosion = TRUE,
  erosion_cells = NULL)
}
\arguments{
  \item{graph_obj}{Output of \code{build_cell_graphs()}; the triangulation is reused and only pruning/wedge parameters are applied here.}
  \item{max_edge_len}{Numeric threshold; Delaunay edges longer than this are dropped. Set to \code{NA}, \code{NULL}, or \code{<= 0} to skip pruning.}
  \item{include_wedge}{Logical; if \code{TRUE}, returns open-triplet (wedge) counts alongside triangles.}
  \item{verbose}{Logical; print progress.}
  \item{offset_pseudo}{Small positive constant used inside offsets.}
  \item{n_cores}{Parallelism hint for motif counting; values greater than 1 trigger \code{parallel::mclapply} on Unix-alikes, while 1 runs sequentially.}
  \item{erosion}{Logical; if \code{TRUE} (default), exclude motifs that touch boundary cells. Boundary cells can be provided as a logical column (e.g. \code{boundary}/\code{is_boundary}) in the input data or via \code{erosion_cells}.}
  \item{erosion_cells}{Optional list of boundary masks/indices by sample name. Each entry can be a logical vector (length = number of cells) or integer indices to exclude.}
}
\value{
A cellgraph object (class \code{cellEdgeR_obj}) with the original graph info plus:
\describe{
  \item{\code{raw_count}}{Sparse matrices for \code{node}, \code{edge}, \code{triangle}, and optional \code{wedge} motifs (rows = motifs, cols = samples).}
  \item{\code{exposure}}{Totals used in offsets: \code{cells}, \code{edges}, \code{triangles}, \code{volumes} per label by sample, \code{center_pairs}, plus \code{wedge} and \code{triples} when requested.}
  \item{\code{offsets}}{Volume offsets containing log-expected matrices per layer. Node offsets are retained for normalization and submotif inspection; tested layers use structural volume offsets.}
  \item{\code{norm_counts}}{Normalized counts for the volume offset (counts / \code{exp(offset)}).}
  \item{\code{relative_counts}}{edgeR intercept-only log2 residuals for the volume offset.}
  \item{\code{offset_part_id}}{Components contributing to each motif's volume offset.}
  \item{\code{offset_part_values}}{Numeric values referenced by \code{offset_part_id} (e.g. volumes and \code{2m}).}
  \item{\code{edger}}{edgeR fits/tests once \code{\link[=motif_edger]{motif_edger()}} is run.}
  \item{\code{parameters}}{Run parameters including \code{max_edge_len}, \code{include_wedge}, \code{offset_pseudo}, \code{offset_modes}, \code{offset_version}, layer names, and \code{node_tmm_offsets}.}
}
}
\examples{
demo <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(demo, verbose = FALSE)
counts <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = FALSE)
counts_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedge = TRUE, verbose = FALSE)
str(counts_full$raw_count)
}
