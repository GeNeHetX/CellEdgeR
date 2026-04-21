\name{top_edges}
\alias{top_edges}
\alias{top_motifs2}
\alias{top_triplets}
\alias{motif_results}
\title{Top motif results}
\description{
Convenience helpers for retrieving ranked motif results from
\code{\link[=motif_edger]{motif_edger()}}. \code{top_edges()} returns edge motifs only.
\code{top_motifs2()} returns wedge and triangle motifs with the corresponding
edge-submotif statistics added as columns. \code{motif_results()} returns both
tables as a named list. \code{top_triplets()} is a deprecated compatibility alias
for \code{top_motifs2()}.
}
\usage{
top_edges(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_motifs2(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_triplets(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH", strategy = NULL)

motif_results(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} with \code{edger} results.}
  \item{coef}{Coefficient name or index; defaults to the first non-intercept coefficient.}
  \item{model}{Which stored model to use: \code{full} or \code{null}.}
  \item{n}{Number of motifs to return; defaults to all.}
  \item{fdr_method}{Multiple testing correction method for \code{p.adjust} (default \code{BH}).}
  \item{strategy}{Deprecated compatibility argument for \code{top_triplets()}; ignored.}
}
\value{
A data frame with columns \code{motif}, \code{motif_type}, \code{logFC},
\code{PValue}, \code{FDR}, and \code{model_used}. \code{top_motifs2()} also includes
\code{edge12}, \code{edge13}, \code{edge23} and their edge-level \code{logFC},
\code{PValue}, and \code{FDR} values.
}
\seealso{
\code{\link[=motif_edger]{motif_edger()}}, \code{\link[=motif_space_size]{motif_space_size()}}
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
top_edges(res, coef = "conditiontreated")
top_motifs2(res, coef = "conditiontreated")
motif_results(res, coef = "conditiontreated")
}
