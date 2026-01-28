\name{get_motif_values}
\alias{get_motif_values}
\title{Retrieve motif values (raw, normalized)}
\description{
Pulls raw counts or normalized counts for a motif (or its lower-order submotifs) into a samples-by-motifs data frame.
}
\usage{
get_motif_values(cellgraph, motif_key = NULL, include_submotifs = FALSE,
  value = c("raw", "norm"))
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} returned by \code{count_motifs_graphs()}.}
  \item{motif_key}{Character vector of motif keys (e.g., \code{E_A_B}, \code{T_A_B_C}). If \code{NULL}, return all motifs.}
  \item{include_submotifs}{Logical; if \code{TRUE}, include lower-order submotifs implied by the motif key (edges for triangles/wedges, nodes for edges).}
  \item{value}{Which value to return: \code{raw} or \code{norm}.}
}
\value{
A data frame with sample names as row names and one column per motif. When
\code{include_submotifs = TRUE}, requested motifs appear first, followed by unique submotifs.
}
\details{
Normalized values are computed directly from the volume offsets so they match the offsets used by
\code{motif_edger()} (including node TMM adjustments when present). Call \code{get_motif_values()}
multiple times to retrieve different value types. Submotifs are derived structurally (edges for
triangles/wedges, nodes for edges) and then filtered to motifs present in the object.
}
\seealso{
\code{\link{count_motifs_graphs}}, \code{\link{motif_edger}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3)
tri_key <- rownames(motifs$raw_count$triangle)[1]
get_motif_values(motifs, motif_key = tri_key, value = "raw")
get_motif_values(motifs, motif_key = tri_key, value = "norm")
}
