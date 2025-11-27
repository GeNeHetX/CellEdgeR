\name{motif_edger}
\alias{motif_edger}
\title{Differential motif testing with edgeR}
\description{
Fits edgeR quasi-likelihood (QL) models for counts of size 1/2/3 motifs, applies per-layer BH correction by default,
and optionally runs the historical DAGGER-style hierarchical FDR (DAG-based) correction.
}
\usage{
motif_edger(motif_obj, sample_df, design_formula, coef = NULL, pseudo = 0.5,
  alpha = 0.05, verbose = TRUE, fdr_method = c("bh", "dagger"))
}
\arguments{
  \item{motif_obj}{Output of \code{count_motifs_graphs()}, containing \code{counts} matrices and exposure totals.}
  \item{sample_df}{Data frame with sample metadata; row names must match \code{motif_obj\$samples}.}
  \item{design_formula}{Formula string compatible with \code{stats::model.matrix()} (e.g. \code{~ condition + batch}).}
  \item{coef}{Index or name of the coefficient to test (non-intercept).}
  \item{pseudo}{Small positive number added inside logarithms when building offsets to avoid \code{log(0)}.}
  \item{alpha}{FDR threshold (used only when \code{fdr_method = "dagger"}).}
  \item{verbose}{Logical; print progress for offsets, edgeR, and (when selected) DAG correction.}
  \item{fdr_method}{Multiple-testing strategy. \code{"bh"} (default) applies BH per layer; \code{"dagger"} reuses the DAG-based refinement once enough samples allow estimation.}
}
\details{
The routine builds offsets derived from expected counts and exposures, fits \pkg{edgeR} GLMs, and attaches BH-adjusted p-values for each motif.
If \code{fdr_method = "dagger"}, a DAG connecting size-1 -> size-2 -> size-3 structures is assembled, and DAGGER-style filtering augments the results with \code{FDR_DAG} columns.
}
\value{
A list with:
\describe{
  \item{\code{design}}{Model matrix used for the fit.}
  \item{\code{results}}{Named list of data frames (\code{size1}, \code{size2}, \code{size3}). Each contains \code{motif}, \code{logFC}, \code{PValue}, and \code{FDR_BH}; the DAGGER branch adds \code{FDR_DAG}.}
  \item{\code{dag}}{NULL when only BH was used; otherwise a list with \code{nodes} and \code{edges} describing the hierarchy.}
  \item{\code{fdr}}{Metadata list with \code{method} and the \code{alpha} threshold applied (only when \code{"dagger"} is selected).}
}
}
\seealso{
\code{\link{count_motifs_graphs}}
}
\examples{
cells <- list(
  s1 = data.frame(x = c(0, 1, 0), y = c(0, 0, 1), label = c("A", "A", "B")),
  s2 = data.frame(x = c(0, 2, 0), y = c(0, 0, 2), label = c("A", "B", "B"))
)
motifs <- count_motifs_graphs(cells, max_edge_len = 3)
sample_df <- data.frame(condition = c("ctrl", "treated"), row.names = motifs$samples)
res <- motif_edger(motifs, sample_df, "~ condition", coef = "conditiontreated")
str(res$results$size2)
}
