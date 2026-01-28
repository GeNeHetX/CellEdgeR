# Create a synthetic dataset for quick testing

Generates a named list of sample data frames with random coordinates and
labels, suitable for exercising the CellEdgeR workflow without real
data. When `boost_motif` is `TRUE`, the second half of samples receives
tightly clustered label triplets to increase a triangle motif (default
labels: a/b/c or repeats when fewer labels exist). The number of
triplets scales with `n\_cells`.

## Usage

``` r
make_demo_samples(n_cells = 300, n_labels = 3, n_samples = 24,
  seed = NULL, boost_motif = TRUE)
```

## Arguments

- n_cells:

  Number of cells per sample.

- n_labels:

  Number of distinct labels (letters starting at "a").

- n_samples:

  Number of samples to generate.

- seed:

  Optional integer seed for reproducibility; if `NULL`, the current RNG
  state is used.

- boost_motif:

  Logical; if `TRUE`, add clustered triplets to enrich a triangle motif
  in the second half of samples.

## Value

A named list of data frames (one per sample) with columns `x`, `y`, and
`label`.

## Examples

``` r
demo_samples <- make_demo_samples(seed = 42)
str(demo_samples[[1]])
#> 'data.frame':    10000 obs. of  3 variables:
#>  $ x    : num  0.915 0.937 0.286 0.83 0.642 ...
#>  $ y    : num  0.528 0.646 0.834 0.346 0.622 ...
#>  $ label: chr  "c" "a" "b" "a" ...
```
