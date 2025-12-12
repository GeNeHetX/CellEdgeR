#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// Build adjacency (0-based), symmetric
static void build_adj(int n, const IntegerVector& ei, const IntegerVector& ej,
                      std::vector< std::vector<int> >& adj) {
  adj.assign(n, std::vector<int>());
  int m = ei.size();
  for (int e = 0; e < m; ++e) {
    int u = ei[e] - 1;
    int v = ej[e] - 1;
    if (u == v) continue;
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  for (int i = 0; i < n; ++i) {
    auto& nb = adj[i];
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
  }
}

// Intersect sorted neighbor lists (only neighbors > j to avoid duplicates)
static int intersect_counting(const std::vector<int>& A, const std::vector<int>& B,
                              int j, std::vector<int>& out) {
  out.clear();
  size_t i = 0, k = 0;
  while (i < A.size() && k < B.size()) {
    if (A[i] == B[k]) {
      if (A[i] > j) out.push_back(A[i]); // only keep > j to enforce i<j<k
      ++i;
      ++k;
    } else if (A[i] < B[k]) {
      ++i;
    } else {
      ++k;
    }
  }
  return static_cast<int>(out.size());
}

//' Count triangles by unordered label triplets on an undirected graph.
//'
//' @param n_nodes Number of nodes.
//' @param ei,ej Edge lists (1-based indices; same length).
//' @param labels Integer labels (1..K) for each node.
//' @param label_ids Sequence 1..K (kept for compatibility).
//'
//' @return A list with `keys` (label triplets as strings), `counts`, and `tri_total`
//'   (total number of triangles).
//'
//' @keywords internal
// [[Rcpp::export(name = "count_triangle_labels_cpp")]]
List count_triangle_labels_cpp(
    const IntegerVector& n_nodes,
    const IntegerVector& ei,
    const IntegerVector& ej,
    const IntegerVector& labels,
    const IntegerVector& label_ids,
    const bool count_wedges
) {
  int n = n_nodes[0];
  if (n <= 2) {
    return List::create(
      _["keys"] = CharacterVector(0),
      _["counts"] = IntegerVector(0),
      _["tri_total"] = 0
    );
  }

  std::vector< std::vector<int> > adj;
  build_adj(n, ei, ej, adj);

  std::unordered_map< std::string, int > counter;
  std::unordered_map< std::string, int > wedge_counter;
  counter.reserve(1024);
  wedge_counter.reserve(1024);

  std::vector<int> tmp;
  long long tri_total = 0;
  long long wedge_total = 0;

  for (int i = 0; i < n; ++i) {
    const auto& Ni = adj[i];
    for (size_t jj = 0; jj < Ni.size(); ++jj) {
      int j = Ni[jj];
      if (j <= i) continue; // enforce i<j
      const auto& Nj = adj[j];
      std::vector<int> cand;
      intersect_counting(Ni, Nj, j, cand); // neighbors k of both i and j, with k>j
      tri_total += cand.size();
      if (cand.empty()) continue;
      int li = labels[i]; // 1..K
      int lj = labels[j];
      for (int kk : cand) {
        int lk = labels[kk];
        int a = li, b = lj, c = lk;
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        std::string key = std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(c);
        auto it = counter.find(key);
        if (it == counter.end()) counter.emplace(key, 1);
        else (++(it->second));
      }
    }

    if (count_wedges) {
      size_t m = Ni.size();
      if (m >= 2) {
        for (size_t ii = 0; ii < m - 1; ++ii) {
          int j = Ni[ii];
          const auto& Nj = adj[j];
          for (size_t jj = ii + 1; jj < m; ++jj) {
            int k = Ni[jj];
            if (std::binary_search(Nj.begin(), Nj.end(), k)) continue;
            ++wedge_total;
            int center = labels[i];
            int leaf1 = labels[j];
            int leaf2 = labels[k];
            int small_leaf = std::min(leaf1, leaf2);
            int large_leaf = std::max(leaf1, leaf2);
            std::string wedge_key = std::to_string(center) + "_" + std::to_string(small_leaf) + "_" + std::to_string(large_leaf);
            auto itw = wedge_counter.find(wedge_key);
            if (itw == wedge_counter.end()) wedge_counter.emplace(wedge_key, 1);
            else ++(itw->second);
          }
        }
      }
    }
  }

  CharacterVector keys(counter.size());
  IntegerVector vals(counter.size());
  int idx = 0;
  for (auto& kv : counter) {
    keys[idx] = kv.first;
    vals[idx] = kv.second;
    ++idx;
  }

  CharacterVector wedge_keys(wedge_counter.size());
  IntegerVector wedge_vals(wedge_counter.size());
  int wedge_idx = 0;
  for (auto& kv : wedge_counter) {
    wedge_keys[wedge_idx] = kv.first;
    wedge_vals[wedge_idx] = kv.second;
    ++wedge_idx;
  }

  return List::create(
    _["keys"] = keys,
    _["counts"] = vals,
    _["tri_total"] = static_cast<double>(tri_total)
    ,
    _["wedge_keys"] = wedge_keys,
    _["wedge_counts"] = wedge_vals,
    _["wedge_total"] = static_cast<double>(wedge_total)
  );
}
