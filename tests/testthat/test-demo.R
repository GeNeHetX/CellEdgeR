test_that("make_demo_samples returns expected shapes and names", {
  res <- make_demo_samples(n_cells = 10, n_labels = 2, n_samples = 4, seed = 123)
  expect_equal(length(res), 4)
  expect_equal(names(res), paste0("s", 1:4))
  expect_true(all(vapply(res, nrow, integer(1)) == 10))
  expect_true(all(vapply(res, function(df) all(c("x", "y", "label") %in% names(df)), logical(1))))
})
