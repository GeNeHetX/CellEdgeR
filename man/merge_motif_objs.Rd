\name{merge_motif_objs}
\alias{merge_motif_objs}
\title{Merge motif objects}
\description{
Combine two \code{cellEdgeR_obj} outputs from \code{count_motifs_graphs()} into a single object
with recomputed offsets and normalized counts.
}
\usage{
merge_motif_objs(motif_obj_a, motif_obj_b, verbose = TRUE)
}
\arguments{
  \item{motif_obj_a,motif_obj_b}{\code{cellEdgeR_obj} objects returned by \code{count_motifs_graphs()}.}
  \item{verbose}{Logical; print progress while recomputing offsets.}
}
\details{
The two objects must have disjoint sample names and overlapping label sets; if the label
sets are disjoint the merge stops with an error. Motif key names are validated against
each object's label levels before merging. Per-sample graphs are retained only when
present in both inputs; label IDs are remapped to the merged label vocabulary.
}
\value{
A merged \code{cellEdgeR_obj} containing the combined samples, with offsets,
normalized counts, relative counts, and offset metadata recomputed on the merged data.
}
\examples{
demo <- make_demo_samples(seed = 1, n_samples = 4)
graphs_full <- build_cell_graphs(demo, verbose = FALSE)
full <- count_motifs_graphs(graphs_full, max_edge_len = NA_real_, verbose = FALSE)

grp <- names(demo)
obj_a <- count_motifs_graphs(build_cell_graphs(demo[grp[1:2]], verbose = FALSE), verbose = FALSE)
obj_b <- count_motifs_graphs(build_cell_graphs(demo[grp[3:4]], verbose = FALSE), verbose = FALSE)
merged <- merge_motif_objs(obj_a, obj_b, verbose = FALSE)
}
