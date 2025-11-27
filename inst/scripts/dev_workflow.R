# Development snippets for CellEdgeR
# Each block can be sourced or run interactively. Nothing destructive is executed
# without an explicit command.

# ---- Environment ----------------------------------------------------------
# setwd("/path/to/CellEdgeR")
# install.packages(c("devtools", "testthat", "roxygen2", "Rcpp"))

# install this package
# Run from the package root directory (the folder containing DESCRIPTION)
if (interactive()) {
  message("Installing CellEdgeR (run from package root)...")
  if (requireNamespace("remotes", quietly = TRUE)) {
    remotes::install_local(".", dependencies = TRUE,force=T)
  } else if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install(upgrade = "never")
  } else {
    stop("Please install 'remotes' or 'devtools' to install the package")
  }
}

# ---- Rcpp glue
 ------------------------------------------------------------
# Run after editing files in src/
if (interactive()) {
  message("Rebuilding Rcpp exports…")
  Rcpp::compileAttributes()
}

# ---- Quick sanity test ----------------------------------------------------
if (interactive()) {
  library(CellEdgeR)
  ncell=100
  mlab=3
  nsam=24

  alab=letters[1:mlab]
  blab=letters[c(1,1,1:mlab)]

  demo <- setNames(lapply(1:nsam,\(i)
    data.frame(x = runif(ncell), y = runif(ncell), label = sample(list(alab,blab)[[1+(i>(nsam/2))]],ncell,replace=T))
  ),paste0("s",1:nsam))

  graphs <- build_cell_graphs(demo, n_cores = 2, verbose = T)
  motifs <- count_motifs_graphs(graph_obj = graphs, max_edge_len = 3, verbose = T)
  motifs_full <- count_motifs_graphs(graph_obj = graphs, max_edge_len = NA_real_, include_wedges = TRUE, verbose = T)

  sample_df <- data.frame(condition = sample(c("ctrl", "treated"),nsam,replace=T),
   row.names = motifs$samples)

  res <- motif_edger(
    motif_obj = motifs,
    sample_df = sample_df,
    design_formula = "~ condition",
    coef = "conditiontreated",
    verbose = T
  )
  print(lapply(res$results, head))
}

if(interactive()){
  imodatL=readRDS(file="~/workspace/DEV/hashGraphMotifs/savingStuff/testdata/eachImotepcellLabel.rds")


    # imodatL2=lapply(imodatL,\(x){colnames(x)=c( "x", "y", "label");x})
    graphs_l <- build_cell_graphs(imodatL, n_cores = 12, verbose = T)

  saveRDS(graphs_l,file="~/workspace/DEV/hashGraphMotifs/savingStuff/testdata/celledger_ImotepGraph.rds")


    library(CellEdgeR)
    imographs_l=readRDS(file="~/workspace/DEV/hashGraphMotifs/savingStuff/testdata/celledger_ImotepGraph.rds")


    motifs <- count_motifs_graphs(graph_obj = imographs_l, max_edge_len = 50, verbose = T,include_wedges = T)

    prog=paste0("s", c(1:4,7,8,9))
    resp=paste0("s",c( 5,6,9,10,11,12))



    normcounts=  normalize_motif_counts(motifs)      

    sample_df <- data.frame(condition = c("prog","resp")[(motifs$samples %in% resp )+1],
    row.names = motifs$samples)

    res <- motif_edger(
      motif_obj = motifs,
      sample_df = sample_df,
      design_formula = "~ condition",
      verbose = T
    )



}

# ---- Testthat + check -----------------------------------------------------
if (interactive()) {
  devtools::test()          # fast unit tests
  # devtools::check()       # full check before release
}

# ---- README / docs --------------------------------------------------------
# README is plain Markdown; knit here if converted to Rmd in the future.
# rmarkdown::render("README.Rmd")

# ---- Git hygiene ----------------------------------------------------------
if (interactive()) {
  system("git status")
  # system('git add .')
  # system('git commit -m \"Add feature or fix\"')
  # system('git push origin main')
}
