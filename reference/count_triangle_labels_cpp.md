# Count triangles by unordered label triplets on an undirected graph.

Count triangles by unordered label triplets on an undirected graph.

## Usage

``` r
count_triangle_labels_cpp(n_nodes, ei, ej, labels, label_ids, count_wedges)
```

## Arguments

- n_nodes:

  Number of nodes.

- ei, ej:

  Edge lists (1-based indices; same length).

- labels:

  Integer labels (1..K) for each node.

- label_ids:

  Sequence 1..K (kept for compatibility).

## Value

A list with `keys` (label triplets as strings), `counts`, and
`tri_total` (total number of triangles).
