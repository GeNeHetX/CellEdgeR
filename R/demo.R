#' Create a synthetic dataset for quick testing
#'
#' Generates a named list of sample data frames with random coordinates and labels,
#' suitable for exercising the CellEdgeR workflow without real data. When
#' \code{boost_motif} is \code{TRUE}, the second half of samples receives tightly
#' clustered label triplets to increase a triangle motif (default labels: a/b/c or repeats
#' when fewer labels exist). The number of triplets scales with \code{n_cells}.
#'
#' @param n_cells Number of cells per sample.
#' @param n_labels Number of distinct labels (letters starting at "a").
#' @param n_samples Number of samples to generate.
#' @param seed Optional integer seed for reproducibility; if \code{NULL}, the current RNG state is used.
#' @param boost_motif Logical; if \code{TRUE}, add clustered triplets to enrich a triangle motif
#'   in the second half of samples.
#'
#' @return A named list of data frames (one per sample) with columns \code{x}, \code{y}, and \code{label}.
#' @examples
#' demo_samples <- make_demo_samples(seed = 42)
#' str(demo_samples[[1]])
#' @export
make_demo_samples <- function(
  n_cells = 10000,
  n_labels = 3,
  n_samples = 24,
  seed = NULL,
  boost_motif = TRUE
) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  if (!is.numeric(n_cells) || length(n_cells) != 1L || n_cells <= 0) {
    stop("n_cells must be a single positive number.")
  }
  if (!is.numeric(n_labels) || length(n_labels) != 1L || n_labels <= 0) {
    stop("n_labels must be a single positive number.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || n_samples <= 0) {
    stop("n_samples must be a single positive number.")
  }
  if (!is.logical(boost_motif) || length(boost_motif) != 1L || is.na(boost_motif)) {
    stop("boost_motif must be TRUE or FALSE.")
  }

  labels_a <- letters[seq_len(n_labels)]
  boost_labels <- if (n_labels >= 3) {
    labels_a[1:3]
  } else if (n_labels == 2) {
    c(labels_a[1], labels_a[2], labels_a[1])
  } else {
    rep(labels_a[1], 3)
  }
  boost_size <- length(boost_labels)
  boost_frac <- 0.6
  boost_jitter <- 0.002

  stats::setNames(
    lapply(seq_len(n_samples), function(i) {
      pool <- labels_a
      df <- data.frame(
        x = stats::runif(n_cells),
        y = stats::runif(n_cells),
        label = sample(pool, n_cells, replace = TRUE),
        stringsAsFactors = FALSE
      )
      if (boost_motif && i > n_samples / 2) {
        n_triplets <- min(floor(n_cells * boost_frac / boost_size), floor(n_cells / boost_size))
        if (n_triplets > 0) {
          idx <- sample(seq_len(n_cells), n_triplets * boost_size)
          for (t in seq_len(n_triplets)) {
            base_x <- stats::runif(1)
            base_y <- stats::runif(1)
            sel <- idx[((t - 1) * boost_size + 1):(t * boost_size)]
            df$x[sel] <- pmin(pmax(base_x + stats::rnorm(boost_size, sd = boost_jitter), 0), 1)
            df$y[sel] <- pmin(pmax(base_y + stats::rnorm(boost_size, sd = boost_jitter), 0), 1)
            df$label[sel] <- boost_labels
          }
        }
      }
      df
    }),
    paste0("s", seq_len(n_samples))
  )
}
