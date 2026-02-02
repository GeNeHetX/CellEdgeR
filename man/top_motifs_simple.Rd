\name{top_motifs_simple}
\alias{top_motifs_simple}
\title{Extract top simple motifs (nodes and edges)}
\description{
Return volume-based differential results for node and edge motifs only.
Use \code{\link{top_motifs_triplet}} for 3-node motifs.
}
\usage{
top_motifs_simple(cellgraph, coef = NULL, model = c("full", "null"),
  n = Inf, fdr_method = "BH")
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} with \code{edger} results.}
  \item{coef}{Coefficient name or index; defaults to the first non-intercept coefficient.}
  \item{model}{Which stored model to use for edgeR strategies: \code{full} or \code{null}.}
  \item{n}{Number of motifs to return; defaults to all.}
  \item{fdr_method}{Multiple testing correction method for \code{p.adjust} (default \code{BH}).}
}
\value{
A data frame with columns: motif, motif_type, logFC, PValue, FDR, model_used.
}
\seealso{
\code{\link{top_motifs_triplet}}, \code{\link{motif_edger}}
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

top_motifs_simple(res, coef = "conditiontreated")
}
