\name{top_motifs}
\alias{top_motifs}
\title{Extract top motifs from stored strategies}
\description{
Returns a ranked data frame of motifs using strategy results stored by \code{\link{motif_edger}}.
}
\usage{
top_motifs(cellgraph, strategy = c("hybrid", "ancova", "volume",
  "hierarchical"), coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH", append_strategies = NULL)
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} with \code{edger} results.}
  \item{strategy}{Which strategy to use: \code{hybrid}, \code{ancova}, \code{volume}, or \code{hierarchical} (defaults to \code{hybrid}; hybrid applies volume to node/edge motifs and ancova to triangle/wedge motifs).}
  \item{coef}{Coefficient name or index; defaults to the first non-intercept coefficient, or the intercept when only an intercept is present.}
  \item{model}{Which stored model to use for edgeR strategies: \code{full} or \code{null}.}
  \item{n}{Number of motifs to return; defaults to all.}
  \item{fdr_method}{Multiple testing correction method for \code{p.adjust} (default \code{BH}).}
  \item{append_strategies}{Optional vector of strategies to append raw PValues/logFC for. Use \code{TRUE} to append all other strategies.}
}
\value{
A data frame with columns: motif, motif\_type, logFC, PValue, FDR, and model\_used.
Additional \code{logFC_<strategy>} and \code{PValue_<strategy>} columns are added when requested.
For \code{hybrid}, FDR is computed separately for node/edge motifs and triangle/wedge motifs.
}
\seealso{
\code{\link{motif_edger}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3)
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ condition")
top_motifs(res, strategy = "volume", coef = "conditiontreated")
}
