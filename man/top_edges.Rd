\name{top_edges}
\alias{top_edges}
\alias{top_triplets}
\title{Top motif results}
\description{
Convenience helpers for retrieving ranked motif results from \code{motif_edger()}.
\code{top_edges()} returns node/edge motifs. \code{top_triplets()} returns 3-node
motifs using either co-occurrence (merged triplets, volume offsets) or topology
(submotif-adjusted, closure mode).
}
\details{
\strong{Co-occurrence (strategy = \code{cooccurrence}):} wedges and triangles are
merged into unordered triplets (\code{TP_*}) and modeled with volume offsets only.
This tests whether a label triplet is over/under represented as a whole, regardless
of closure. Computed by default when \code{motif_edger()} is run with
\code{include_wedge = TRUE}; explicitly use \code{triplet_mode = "merge"} to recompute
only this mode.

\strong{Topology (strategy = \code{topology}):} wedges and triangles are modeled
separately in closure mode. Wedges are \emph{submotif-adjusted} using edge-derived
covariates (edge-force), and triangles are tested for closure using total triples
as a covariate (submotif-adjusted). Computed by default when \code{motif_edger()}
is run with \code{include_wedge = TRUE}; explicitly use
\code{triplet_mode = "closure", strategies = "submotif_adj"} to recompute only this mode.

\strong{Model selection (full vs null):} when \code{strategy = "cooccurrence"},
\code{model = "full"} uses the full design formula, and \code{model = "null"} uses
an intercept-only model (useful for diagnostics or when you want the baseline rate).
Topology uses the submotif-adjusted model and ignores \code{model}.
}
\usage{
top_edges(cellgraph, coef = NULL, model = c("full", "null"), n = Inf,
  fdr_method = "BH")

top_triplets(cellgraph, strategy = c("cooccurrence", "topology"),
  coef = NULL, model = c("full", "null"), n = Inf, fdr_method = "BH")
}
\arguments{
  \item{cellgraph}{A \code{cellEdgeR_obj} with \code{edger} results.}
  \item{coef}{Coefficient name or index; defaults to the first non-intercept coefficient.}
  \item{model}{Which stored model to use for edgeR strategies: \code{full} or \code{null}.
  \code{full} uses the provided design formula; \code{null} uses an intercept-only
  model. Ignored when \code{strategy = "topology"}.}
  \item{n}{Number of motifs to return; defaults to all.}
  \item{fdr_method}{Multiple testing correction method for \code{p.adjust} (default \code{BH}).}
  \item{strategy}{Triplet strategy: \code{cooccurrence} (volume offsets, merged triplets)
  or \code{topology} (submotif-adjusted, closure mode).}
}
\value{
A data frame with columns: motif, motif\_type, logFC, PValue, FDR, and model\_used.
}
\seealso{
\code{\link{motif_edger}}, \code{\link{motif_space_size}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
graphs <- build_cell_graphs(cells, verbose = FALSE)
motifs <- count_motifs_graphs(graphs, max_edge_len = 3, include_wedge = TRUE)
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$sample_name)

res_merge <- motif_edger(motifs, sample_df, "~ condition", triplet_mode = "merge")
res_close <- motif_edger(motifs, sample_df, "~ condition", triplet_mode = "closure")

top_edges(res_merge, coef = "conditiontreated")
top_triplets(res_merge, strategy = "cooccurrence", coef = "conditiontreated")
top_triplets(res_close, strategy = "topology", coef = "conditiontreated")
}
