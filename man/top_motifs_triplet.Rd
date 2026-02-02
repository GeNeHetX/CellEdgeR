\name{top_motifs_triplet}
\alias{top_motifs_triplet}
\title{Extract top 3-node motifs}
\description{
Return differential results for 3-node motifs only. When \code{triplet_mode = "merge"},
results are returned for unordered triplets. Otherwise, triangles and wedges are
returned separately (when available).
}
\usage{
top_motifs_triplet(cellgraph, strategy = c("ancova", "volume"),
  coef = NULL, model = c("full", "null"), n = Inf, fdr_method = "BH",
  triplet_mode = NULL)
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} with \code{edger} results.}
  \item{strategy}{Which strategy to use: \code{volume} or \code{ancova}.}
  \item{coef}{Coefficient name or index; defaults to the first non-intercept coefficient.}
  \item{model}{Which stored model to use for edgeR strategies: \code{full} or \code{null}.}
  \item{n}{Number of motifs to return; defaults to all.}
  \item{fdr_method}{Multiple testing correction method for \code{p.adjust} (default \code{BH}).}
  \item{triplet_mode}{Triplet handling mode (\code{separate}, \code{merge}, or \code{closure}).
  Defaults to the mode stored in \code{cellgraph$edger$triplet_mode}.}
}
\value{
A data frame with columns: motif, motif_type, logFC, PValue, FDR, model_used.
}
\seealso{
\code{\link{top_motifs_simple}}, \code{\link{motif_edger}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C")),
  s2 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "B", "C"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = NA_real_, include_wedge = TRUE)
sample_df <- data.frame(group = c("g1", "g2"), row.names = motifs$sample_name)
res <- motif_edger(motifs, sample_df, "~ group", triplet_mode = "merge")

top_motifs_triplet(res, strategy = "volume", coef = "groupg2", triplet_mode = "merge")
}
