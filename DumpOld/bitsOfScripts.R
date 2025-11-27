# install.packages("Rcpp")
library(Rcpp)

# Count triangles by UNORDERED label triplets on an undirected graph given as 1-based edge list.
# Also returns total triangle count (exposure).
cppFunction(code = '
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// make adjacency (0-based), symmetric
static void build_adj(int n, const IntegerVector& ei, const IntegerVector& ej,
                      std::vector< std::vector<int> >& adj) {
  adj.assign(n, std::vector<int>());
  int m = ei.size();
  for (int e=0;e<m;++e){
    int u = ei[e]-1, v = ej[e]-1;
    if (u==v) continue;
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  for (int i=0;i<n;++i){
    auto &nb = adj[i];
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
  }
}

// intersect sorted neighbor lists (only neighbors > j to avoid dups)
static int intersect_counting(const std::vector<int>& A, const std::vector<int>& B,
                              int j, std::vector<int>& out){
  out.clear();
  size_t i=0,k=0;
  while (i<A.size() && k<B.size()){
    if (A[i]==B[k]){
      if (A[i] > j) out.push_back(A[i]); // only keep > j to enforce i<j<k
      ++i; ++k;
    } else if (A[i] < B[k]) {
      ++i;
    } else {
      ++k;
    }
  }
  return (int)out.size();
}

// [[Rcpp::export]]
List count_triangle_labels(
    const IntegerVector& n_nodes,
    const IntegerVector& ei,
    const IntegerVector& ej,
    const IntegerVector& labels,    // 1..K integers, length = n_nodes
    const IntegerVector& label_ids  // 1..K (ensures compact range)
){
  int n = n_nodes[0];
  if (n <= 2) {
    return List::create(
      _["keys"] = CharacterVector(0),
      _["counts"] = IntegerVector(0),
      _["tri_total"] = 0
    );
  }
  // adjacency
  std::vector< std::vector<int> > adj;
  build_adj(n, ei, ej, adj);

  // map triangle label triplets -> count (unordered a<=b<=c)
  std::unordered_map< std::string, int > counter;
  counter.reserve(1024);

  std::vector<int> tmp;
  long long tri_total = 0;

  for (int i=0;i<n;++i){
    const auto &Ni = adj[i];
    for (size_t jj=0; jj<Ni.size(); ++jj){
      int j = Ni[jj];
      if (j <= i) continue; // enforce i<j
      const auto &Nj = adj[j];
      std::vector<int> cand;
      intersect_counting(Ni, Nj, j, cand); // neighbors k of both i and j, with k>j
      tri_total += cand.size();
      if (cand.empty()) continue;
      // label at i and j
      int li = labels[i]; // 1..K
      int lj = labels[j];
      for (int kk : cand){
        int lk = labels[kk];
        // unordered key a|b|c with a<=b<=c in LABEL ID order
        int a = li, b = lj, c = lk;
        if (a > b) std::swap(a,b);
        if (b > c) std::swap(b,c);
        if (a > b) std::swap(a,b);
        // build key string "a|b|c"
        std::string key = std::to_string(a) + "|" + std::to_string(b) + "|" + std::to_string(c);
        auto it = counter.find(key);
        if (it == counter.end()) counter.emplace(key, 1);
        else (++(it->second));
      }
    }
  }

  // dump to vectors
  CharacterVector keys; keys.reserve(counter.size());
  IntegerVector vals(counter.size());
  int idx=0;
  for (auto &kv : counter){
    keys.push_back(kv.first);
    vals[idx++] = kv.second;
  }

  return List::create(
    _["keys"] = keys,
    _["counts"] = vals,
    _["tri_total"] = (double)tri_total
  );
}
')



# Depends: deldir, Matrix, data.table, igraph (for edges only)
count_motifs_graphs <- function(
  cells_by_sample,                 # named list of data.frames: x, y, label (character or factor)
  max_edge_len,                    # prune Delaunay edges longer than this (in same units as x,y)
  use_parallel = FALSE,            # optional: parallel graph builds across samples
  n_cores = max(1, parallel::detectCores() - 1),
  verbose = TRUE
){
  requireNamespace("deldir"); requireNamespace("Matrix"); requireNamespace("data.table"); requireNamespace("igraph")
  library(data.table)

  samples <- names(cells_by_sample)
  if (is.null(samples) || length(samples)==0) stop("cells_by_sample must be a *named* list.")
  if (verbose) message("Samples: ", paste(samples, collapse = ", "))

  # 1) Standardize labels to integers 1..K globally (vectorized)
  all_labels_chr <- unlist(lapply(cells_by_sample, function(d) as.character(d$label)), use.names = FALSE)
  lab_levels <- sort(unique(all_labels_chr))
  K <- length(lab_levels)
  lab_to_id <- setNames(seq_len(K), lab_levels)

  # helper: build Delaunay, prune edges, return edge list (2-col matrix 1-based)
  build_edges <- function(df){
    tri <- deldir::deldir(df$x, df$y)
    tt  <- deldir::triangles(tri)
    if (NROW(tt)==0) return(matrix(numeric(0), ncol=2))
    idx <- as.data.frame(tt[, c("i1","i2","i3")])
    e   <- unique(rbind(idx[,1:2], idx[,c(1,3)], idx[,2:3]))
    colnames(e) <- c("from","to")
    xy <- as.matrix(df[, c("x","y")])
    len <- sqrt(rowSums( (xy[e[,1],] - xy[e[,2],])^2 ))
    e[len <= max_edge_len, , drop = FALSE]
  }

  # optional parallelism on graph building
  map_fun <- if (use_parallel && .Platform$OS.type != "windows") parallel::mclapply else lapply

  if (verbose) message("Building Delaunay + pruning (", if (use_parallel && .Platform$OS.type!="windows") paste0(n_cores, " cores") else "sequential", ")…")
  per_sample <- map_fun(samples, function(s){
    df <- as.data.frame(cells_by_sample[[s]])
    if (!all(c("x","y","label") %in% names(df))) stop("Each data.frame must contain x,y,label")
    labs_chr <- as.character(df$label)
    labs_id  <- unname(lab_to_id[labs_chr])
    edges    <- build_edges(df)
    list(sample=s, n=nrow(df), labels_id=labs_id, labels_chr=labs_chr, edges=edges)
  }, mc.cores = n_cores)

  # 2) Vectorized singles (size-1) per sample
  if (verbose) message("Counting size-1 (cells by label)…")
  sing_list <- lapply(per_sample, function(ps){
    if (ps$n==0) return(integer(0))
    tab <- tabulate(ps$labels_id, nbins = K)
    names(tab) <- lab_levels
    tab
  })
  # Stack into matrix motifs x samples
  Y1 <- do.call(cbind, sing_list)
  dimnames(Y1) <- list(lab_levels, samples)
  Y1 <- Matrix::Matrix(Y1, sparse = TRUE)  # dgCMatrix

  # 3) Vectorized pairs (size-2) per sample (table over unordered label pairs)
  if (verbose) message("Counting size-2 (edges by unordered label pair)…")
  pairs_all_keys <- character(0)
  pairs_by_sample <- vector("list", length(samples))
  names(pairs_by_sample) <- samples

  for (idx in seq_along(per_sample)) {
    ps <- per_sample[[idx]]
    s  <- ps$sample
    e  <- ps$edges
    if (length(e)==0) {
      pairs_by_sample[[s]] <- integer(0)
      next
    }
    la <- ps$labels_chr[e[,1]]
    lb <- ps$labels_chr[e[,2]]
    # unordered keys with vectorized pmin/pmax on factors not safe -> use character then sort
    key <- ifelse(la <= lb, paste(la, lb, sep="|"), paste(lb, la, sep="|"))
    tab <- sort(tapply(rep(1L, length(key)), key, sum))  # vectorized sum by key
    pairs_all_keys <- union(pairs_all_keys, names(tab))
    pairs_by_sample[[s]] <- tab
  }
  # build Y2 with all keys
  if (length(pairs_all_keys)){
    Y2 <- Matrix::Matrix(0, nrow = length(pairs_all_keys), ncol = length(samples),
                         sparse = TRUE, dimnames = list(pairs_all_keys, samples))
    for (s in samples){
      tab <- pairs_by_sample[[s]]
      if (length(tab)) {
        Y2[names(tab), s] <- as.integer(tab)
      }
    }
  } else {
    Y2 <- Matrix::Matrix(0, nrow = 0, ncol = length(samples), sparse = TRUE, dimnames = list(character(), samples))
  }

  # 4) Triangles (size-3) via Rcpp (fast, no R loops over triangles)
  if (verbose) message("Counting size-3 (triangles by unordered label triplet) in C++…")
  tris_all_keys <- character(0)
  tris_counts_list <- vector("list", length(samples))
  names(tris_counts_list) <- samples
  tri_totals <- numeric(length(samples)); names(tri_totals) <- samples

  for (idx in seq_along(per_sample)) {
    ps <- per_sample[[idx]]
    s  <- ps$sample
    if (nrow(ps$edges)==0 || ps$n < 3) {
      tris_counts_list[[s]] <- integer(0); tri_totals[s] <- 0
      next
    }
    out <- count_triangle_labels(n_nodes = ps$n,
                                 ei = ps$edges[,1],
                                 ej = ps$edges[,2],
                                 labels = ps$labels_id,
                                 label_ids = seq_len(K))
    keys <- as.character(out$keys)
    if (length(keys)){
      # convert "id|id|id" back to label strings (vectorized)
      # split once, map ids to names, rebuild string (vectorized)
      ids_mat <- do.call(rbind, strsplit(keys, "\\|"))
      lab_mat <- matrix(lab_levels[as.integer(ids_mat)], ncol = 3)
      keys_lab <- apply(lab_mat, 1, function(v) paste(v, collapse="|"))
      counts <- as.integer(out$counts)
      # accumulate
      tris_all_keys <- union(tris_all_keys, keys_lab)
      tris_counts_list[[s]] <- structure(counts, names = keys_lab)
      tri_totals[s] <- as.numeric(out$tri_total)
    } else {
      tris_counts_list[[s]] <- integer(0); tri_totals[s] <- 0
    }
  }
  if (length(tris_all_keys)){
    Y3 <- Matrix::Matrix(0, nrow = length(tris_all_keys), ncol = length(samples),
                         sparse = TRUE, dimnames = list(tris_all_keys, samples))
    for (s in samples){
      tab <- tris_counts_list[[s]]
      if (length(tab)) Y3[names(tab), s] <- as.integer(tab)
    }
  } else {
    Y3 <- Matrix::Matrix(0, nrow = 0, ncol = length(samples), sparse = TRUE, dimnames = list(character(), samples))
  }

  # exposures (totals) per sample
  # total edges = length of edge list; total triangles from C++
  n_edges <- vapply(per_sample, function(ps) nrow(ps$edges), 0L)
  names(n_edges) <- samples
  n_tris  <- tri_totals

  if (verbose){
    message("Counts ready: |labels|=", K,
            ", singles=", nrow(Y1), ", pairs=", nrow(Y2), ", triangles=", nrow(Y3))
  }

  list(
    samples = samples,
    label_levels = lab_levels,
    counts = list(size1 = Y1, size2 = Y2, size3 = Y3),
    exposure = list(edges = n_edges, triangles = n_tris, cells = Matrix::colSums(Y1)),
    meta = list(max_edge_len = max_edge_len)
  )
}


# Depends: edgeR, Matrix, data.table
motif_edger <- function(
  motif_obj,            # output list from count_motifs_graphs()
  sample_df,            # data.frame; rownames = samples; columns for phenotype/covariates
  design_formula,       # e.g. ~ phenotype + batch
  coef = NULL,          # target coefficient (name or index)
  pseudo = 0.5,         # avoids log(0) in offsets
  alpha = 0.05,         # FDR level for DAGGER (if used)
  verbose = TRUE
){
  requireNamespace("edgeR"); requireNamespace("Matrix"); requireNamespace("data.table")

  samples <- motif_obj$samples
  if (is.null(rownames(sample_df))) stop("sample_df must have rownames = sample IDs")
  sample_df <- as.data.frame(sample_df[samples, , drop = FALSE])
  design <- model.matrix(as.formula(design_formula), data = sample_df)
  if (is.null(coef)){
    coef <- which(colnames(design)!="(Intercept)")[1]
    if (length(coef)==0) stop("No non-intercept term found; set `coef`.")
  } else if (is.character(coef)){
    coef <- match(coef, colnames(design)); if (is.na(coef)) stop("coef name not found.")
  }

  Y1 <- motif_obj$counts$size1
  Y2 <- motif_obj$counts$size2
  Y3 <- motif_obj$counts$size3

  # exposures (vectorized)
  log_cells <- matrix(log(pmax(as.numeric(motif_obj$exposure$cells), 1)), nrow = max(1,nrow(Y1)),
                      ncol = length(samples), byrow = TRUE, dimnames = list(rownames(Y1), samples))
  log_edges <- log(pmax(as.numeric(motif_obj$exposure$edges), 1))
  names(log_edges) <- samples
  log_tris  <- log(pmax(as.numeric(motif_obj$exposure$triangles), 1))
  names(log_tris) <- samples

  # --------- Vectorized OFFSETS
  # Pairs offsets: log N(A) + log N(B) + log(exposure_edges)
  build_pair_offsets <- function(Y2, Y1, log_edges, pseudo){
    if (nrow(Y2)==0) return(Matrix::Matrix(0, nrow = 0, ncol = ncol(Y2),
                                           sparse = TRUE, dimnames = list(character(), colnames(Y2))))
    # parse rownames once
    ab <- do.call(rbind, strsplit(rownames(Y2), "\\|"))
    a <- ab[,1]; b <- ab[,2]
    # pull the two single rows as dense matrices (fast because rows are sparse)
    NA_ <- ifelse(a %in% rownames(Y1), as.matrix(Y1[a, , drop = FALSE]), 0)
    NB_ <- ifelse(b %in% rownames(Y1), as.matrix(Y1[b, , drop = FALSE]), 0)
    # log with pseudo, then add exposures (recycle by row)
    off <- log(pmax(NA_, pseudo)) + log(pmax(NB_, pseudo))
    off <- sweep(off, 2, log_edges[colnames(Y2)], FUN = "+")
    dimnames(off) <- dimnames(Y2)
    off
  }

  # Triangles offsets: log E(AB)+log E(AC)+log E(BC) + log(exposure_tris)
  build_tri_offsets <- function(Y3, Y2, log_tris, pseudo){
    if (nrow(Y3)==0) return(Matrix::Matrix(0, nrow = 0, ncol = ncol(Y3),
                                           sparse = TRUE, dimnames = list(character(), colnames(Y3))))
    abc <- do.call(rbind, strsplit(rownames(Y3), "\\|"))
    a <- abc[,1]; b <- abc[,2]; c <- abc[,3]
    pair_key <- function(u,v) if (u<=v) paste(u,v,sep="|") else paste(v,u,sep="|")
    ABk <- pair_key(a,b); ACk <- pair_key(a,c); BCk <- pair_key(b,c)
    get_pair_mat <- function(keys){
      idx <- match(keys, rownames(Y2))
      # fast extraction: missing rows become zeros
      out <- matrix(0, nrow = length(keys), ncol = ncol(Y2), dimnames = list(keys, colnames(Y2)))
      sel <- which(!is.na(idx))
      if (length(sel)) out[sel, ] <- as.matrix(Y2[idx[sel], , drop = FALSE])
      out
    }
    AB <- get_pair_mat(ABk); AC <- get_pair_mat(ACk); BC <- get_pair_mat(BCk)
    off <- log(pmax(AB, pseudo)) + log(pmax(AC, pseudo)) + log(pmax(BC, pseudo))
    off <- sweep(off, 2, log_tris[colnames(Y3)], FUN = "+")
    dimnames(off) <- dimnames(Y3)
    off
  }

  if (verbose) message("Building offsets (vectorized)…")
  off2 <- build_pair_offsets(Y2, Y1, log_edges, pseudo)
  off3 <- build_tri_offsets(Y3, Y2, log_tris,  pseudo)

  # --------- edgeR fits (per layer; no loops over rows)
  fit_layer <- function(Y, off){
    if (nrow(Y)==0) return(NULL)
    dge <- edgeR::DGEList(counts = Y)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    fit <- edgeR::glmQLFit(dge, design = design, offset = off)
    tst <- edgeR::glmQLFTest(fit, coef = coef)
    tt  <- edgeR::topTags(tst, n = Inf, sort.by = "none")$table
    data.frame(motif = rownames(Y), logFC = tt$logFC, PValue = tt$PValue,
               FDR_BH = p.adjust(tt$PValue, method = "BH"), row.names = NULL)
  }

  if (verbose) message("Fitting edgeR (QL)…")
  res1 <- fit_layer(Y1, log_cells)
  res2 <- fit_layer(Y2, off2)
  res3 <- fit_layer(Y3, off3)

  # --------- DAG FDR (DAGGER if available; otherwise hierarchical BH fallback)
  if (verbose) message("Applying DAG FDR…")
  size1 <- if (!is.null(res1)) res1$motif else character()
  size2 <- if (!is.null(res2)) res2$motif else character()
  size3 <- if (!is.null(res3)) res3$motif else character()
  node_order <- c(size1, size2, size3)
  node_id <- setNames(seq_along(node_order), node_order)

  # build edges (parents -> children)
  edges <- list()
  if (length(size1) && length(size2)){
    ab <- do.call(rbind, strsplit(size2, "\\|"))
    for (i in seq_along(size2)){
      a <- ab[i,1]; b <- ab[i,2]
      if (a %in% size1) edges[[length(edges)+1]] <- c(node_id[[a]], node_id[[size2[i]]])
      if (b %in% size1) edges[[length(edges)+1]] <- c(node_id[[b]], node_id[[size2[i]]])
    }
  }
  if (length(size2) && length(size3)){
    abc <- do.call(rbind, strsplit(size3, "\\|"))
    pair_key <- function(u,v) if (u<=v) paste(u,v,sep="|") else paste(v,u,sep="|")
    for (i in seq_along(size3)){
      A <- abc[i,1]; B <- abc[i,2]; C <- abc[i,3]
      for (pp in c(pair_key(A,B), pair_key(A,C), pair_key(B,C))){
        if (pp %in% size2) edges[[length(edges)+1]] <- c(node_id[[pp]], node_id[[size3[i]]])
      }
    }
  }
  E_dag <- if (length(edges)) do.call(rbind, edges) else matrix(numeric(0), ncol = 2)

  # collect p-values
  p_all <- rep(NA_real_, length(node_order))
  if (length(size1)) p_all[node_id[size1]] <- res1$PValue
  if (length(size2)) p_all[node_id[size2]] <- res2$PValue
  if (length(size3)) p_all[node_id[size3]] <- res3$PValue

  # Conservative hierarchical BH on the DAG
  q_dag <- rep(NA_real_, length(node_order))
  if (verbose) message("Applying hierarchical BH on DAG.")
  q1 <- if (length(size1)) p.adjust(res1$PValue, "BH") else numeric(0)
  keep1 <- if (length(q1)) (q1 <= alpha) else logical(0)
  # size2: require both singletons significant
  q2 <- rep(NA_real_, length(size2))
  if (length(size2)){
    ab <- do.call(rbind, strsplit(size2, "\\|"))
    ok <- logical(length(size2))
    for (i in seq_along(size2)){
      ok[i] <- (ab[i,1] %in% size1[keep1]) && (ab[i,2] %in% size1[keep1])
    }
    pv2 <- res2$PValue; pv2[!ok] <- NA
    if (any(ok)) q2[ok] <- p.adjust(pv2[ok], "BH")
  }
  # size3: require all three parent pairs significant
  q3 <- rep(NA_real_, length(size3))
  if (length(size3)){
    abc <- do.call(rbind, strsplit(size3, "\\|"))
    pair_key <- function(u,v) if (u<=v) paste(u,v,sep="|") else paste(v,u,sep="|")
    ok3 <- logical(length(size3))
    keep2_set <- if (length(size2)) size2[which(q2 <= alpha)] else character(0)
    for (i in seq_along(size3)){
      prs <- c(pair_key(abc[i,1],abc[i,2]), pair_key(abc[i,1],abc[i,3]), pair_key(abc[i,2],abc[i,3]))
      ok3[i] <- all(prs %in% keep2_set)
    }
    pv3 <- res3$PValue; pv3[!ok3] <- NA
    if (any(ok3)) q3[ok3] <- p.adjust(pv3[ok3], "BH")
  }
  if (length(size1)) q_dag[node_id[size1]] <- q1
  if (length(size2)) q_dag[node_id[size2]] <- q2
  if (length(size3)) q_dag[node_id[size3]] <- q3

  add_q <- function(res, layer_keys){
    if (is.null(res) || !length(layer_keys)) return(res)
    res$FDR_DAG <- q_dag[node_id[layer_keys]]
    res
  }
  res1 <- add_q(res1, size1)
  res2 <- add_q(res2, size2)
  res3 <- add_q(res3, size3)

  if (verbose) {
    message("Done. DAG FDR: hierarchical BH.")
  }

  list(
    design = design,
    results = list(size1 = res1, size2 = res2, size3 = res3),
    dag = list(nodes = node_order, edges = E_dag),
    alpha = alpha
  )
}
