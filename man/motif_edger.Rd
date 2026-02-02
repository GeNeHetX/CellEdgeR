\name{motif_edger}
\alias{motif_edger}
\title{Differential motif testing with multiple strategies}
\description{
Fits multiple edgeR-based strategies across stacked motifs and stores results for later extraction.
}
\usage{
motif_edger(cellgraph, sample_df, design_formula, verbose = TRUE,
  triplet_mode = c("separate", "merge", "closure"),
  strategies = c("volume", "ancova"))
}
\arguments{
  \item{cellgraph}{Output of \code{count_motifs_graphs()}, containing motif counts and offsets.}
  \item{sample_df}{Data frame with sample metadata; row names must match \code{cellgraph\$sample_name}.}
  \item{design_formula}{Formula string compatible with \code{stats::model.matrix()} (e.g. \code{~ condition + batch}).}
  \item{verbose}{Logical; print progress for edgeR.}
  \item{triplet_mode}{How to handle 3-node motifs: \code{separate} keeps triangle and wedge motifs (default),
  \code{merge} combines wedges+triangles into unordered triplet motifs, and \code{closure} models triangles with
  total-triplet adjustment. \code{merge} and \code{closure} require \code{count_motifs_graphs(..., include_wedge = TRUE)}.}
  \item{strategies}{Character vector of strategies to run; defaults to \code{volume} and \code{ancova}.}
}
\details{
Two strategies are supported by default:
\itemize{
  \item \code{volume}: volume offsets (Chung-Lu baseline)
  \item \code{ancova}: volume offsets plus edge-force covariates for 3-node motifs
}
Use \code{\link{top_motifs_simple}} for node/edge motifs and
\code{\link{top_motifs_triplet}} for 3-node motifs. The hybrid summary in
\code{\link{top_motifs}} combines volume (nodes/edges) and ancova (3-node motifs).
}
\value{
The input cellgraph augmented with \code{edger}, containing:
\describe{
  \item{\code{strategies}}{Named list of strategy results (\code{volume}, \code{ancova}).}
  \item{\code{motif_info}}{Data frame of motif keys and motif types for joins.}
  \item{\code{sample_df}}{Sample metadata used to build the model matrix.}
  \item{\code{triplet_mode}}{Triplet handling mode used for 3-node motifs.}
}
}
\seealso{
\code{\link{count_motifs_graphs}}, \code{\link{top_motifs}}, \code{\link{top_motifs_simple}},
\code{\link{top_motifs_triplet}}
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
head(top_motifs_simple(res, coef = "conditiontreated"))
head(top_motifs_triplet(res, strategy = "ancova", coef = "conditiontreated"))
}
