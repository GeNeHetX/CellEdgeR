\name{motif_edger}
\alias{motif_edger}
\title{Differential motif testing with edgeR}
\description{
Fit a single volume-offset edgeR QL model across edge, wedge, and triangle
motifs. Node counts are used internally to build graph volumes but are not
tested by \code{motif_edger()}. Wedges and triangles are kept as separate motifs;
there is no merged triplet test. Low-count motifs are filtered with
\code{edgeR::filterByExpr()} before model fitting, so motifs carried by too few
samples are not tested.
}
\usage{
motif_edger(cellgraph, sample_df, design_formula, verbose = TRUE)
}
\arguments{
  \item{cellgraph}{Output list from \code{\link[=count_motifs_graphs]{count_motifs_graphs()}}.}
  \item{sample_df}{Data frame with sample metadata; rownames must match \code{cellgraph$sample_name}.}
  \item{design_formula}{Formula or formula string passed to \code{model.matrix}, e.g. \code{~ condition + batch}.}
  \item{verbose}{Logical; print progress.}
}
\value{
The input \code{cellgraph} augmented with \code{edger}, containing one stored
strategy named \code{volume}. Use \code{\link[=top_edges]{top_edges()}} for edge results and
\code{\link[=top_edges]{top_motifs2()}} for wedge/triangle results with edge-submotif statistics.
Filter details are stored in \code{cellgraph$edger$filter}.
}
\details{
The model is volume based by default. It stacks edge, triangle, and wedge
motifs and fits edgeR quasi-likelihood tests using the volume offsets stored by
\code{count_motifs_graphs()}. The output result sets are \code{edges} and
\code{motifs2}. Motifs filtered by \code{edgeR::filterByExpr()} do not receive
p-values or FDR values.
}
\seealso{
\code{\link[=count_motifs_graphs]{count_motifs_graphs()}}, \code{\link[=top_edges]{top_edges()}}, \code{\link[=top_edges]{top_motifs2()}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3, include_wedge = TRUE)
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ condition")
head(top_edges(res, coef = "conditiontreated"))
head(top_motifs2(res, coef = "conditiontreated"))
}
