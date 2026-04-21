\name{motif_space_size}
\alias{motif_space_size}
\title{Count the number of possible tested motif labels}
\description{
Compute the combinatorial number of possible edge, wedge, and triangle label
motifs implied by a \code{cellEdgeR_obj}. Counts reflect label combinations, not
graph isomorphism classes. Node motifs and merged triplets are not part of the
tested motif space.
}
\usage{
motif_space_size(cellgraph, include_wedge = NULL)
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} (typically after \code{\link[=motif_edger]{motif_edger()}}).}
  \item{include_wedge}{Logical; whether to include wedge motifs in the count.
  Defaults to \code{cellgraph$parameters$include_wedge} when available.}
}
\value{
A list with \code{labels}, \code{include_wedge}, \code{counts} (data frame), and
\code{total} (sum of possible tested motifs).
}
\seealso{
\code{\link[=motif_edger]{motif_edger()}}, \code{\link[=top_edges]{top_edges()}}, \code{\link[=top_edges]{top_motifs2()}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
  s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
)
obj <- count_motifs_graphs(build_cell_graphs(cells), include_wedge = TRUE)
meta <- data.frame(group = c("g1", "g2"), row.names = obj$sample_name)
obj <- motif_edger(obj, meta, "~ group")

motif_space_size(obj)
}
