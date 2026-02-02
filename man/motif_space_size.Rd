\name{motif_space_size}
\alias{motif_space_size}
\title{Count possible motif labels}
\description{
Compute the combinatorial number of possible motif labelings implied by a
\code{cellEdgeR_obj}. Counts reflect label combinations, not graph isomorphism classes.
}
\usage{
motif_space_size(cellgraph, triplet_mode = NULL, include_wedge = NULL)
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} (typically after \code{\link{motif_edger}}).}
  \item{triplet_mode}{Triplet handling mode (\code{separate}, \code{merge}, or \code{closure}).
  Defaults to the mode stored in \code{cellgraph$edger$triplet_mode}.}
  \item{include_wedge}{Logical; whether to include wedge motifs in the count. Defaults to
  \code{cellgraph$parameters$include_wedge} when available.}
}
\value{
A list with \code{labels}, \code{triplet_mode}, \code{include_wedge}, \code{counts}
(a data frame of layers and counts), and \code{total} (sum of possible motifs).
}
\seealso{
\code{\link{motif_edger}}, \code{\link{top_motifs_simple}}, \code{\link{top_motifs_triplet}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
  s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
)
obj <- count_motifs_graphs(build_cell_graphs(cells), include_wedge = TRUE)
meta <- data.frame(group = c("g1", "g2"), row.names = obj$sample_name)
obj <- motif_edger(obj, meta, "~ group", triplet_mode = "merge")

motif_space_size(obj)
}
