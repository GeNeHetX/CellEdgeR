#' Create a synthetic dataset for quick testing
#'
#' Generates a named list of sample data frames with random coordinates and labels,
#' suitable for exercising the CellEdgeR workflow without real data.
#'
#' @param n_cells Number of cells per sample.
#' @param n_labels Number of distinct labels (letters starting at "a").
#' @param n_samples Number of samples to generate.
#' @param seed Optional integer seed for reproducibility; if \code{NULL}, the current RNG state is used.
#'
#' @return A named list of data frames (one per sample) with columns \code{x}, \code{y}, and \code{label}.
#' @examples
#' demo_samples <- make_demo_samples(seed = 42)
#' str(demo_samples[[1]])
#' @export
make_demo_samples <- function(n_cells = 100, n_labels = 3, n_samples = 24, seed = NULL) {
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

  labels_a <- letters[seq_len(n_labels)]
  labels_b <- c(labels_a[1], labels_a[1], labels_a)

  stats::setNames(
    lapply(seq_len(n_samples), function(i) {
      pool <- if (i > n_samples / 2) labels_b else labels_a
      data.frame(
        x = stats::runif(n_cells),
        y = stats::runif(n_cells),
        label = sample(pool, n_cells, replace = TRUE),
        stringsAsFactors = FALSE
      )
    }),
    paste0("s", seq_len(n_samples))
  )
}
