# tests/testthat/test-set_pg.R

test_that("set_pg returns correct matrix dimensions", {

  intv <- sample(10, size = 1)
  s <- sample(10, size = 1)
  warmup <- s * intv
  n_sim <- 2 * warmup
  pg <- par.control(pg = "sync", intv = intv, min = 0, max = 1, seed = 123)

  out <- set_pg(warmup, n_sim, s, pg)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n_sim, s))
})

test_that("set_pg errors if pg is not a list", {
  expect_error(
    set_pg(warmup = 5, n_sim = 6, s = 3, pg = "sync"),
    "'pg' must be a list"
  )
})

test_that("set_pg sync pattern fills matrix correctly", {

  intv <- sample(10, size = 1)
  s <- sample(10, size = 1)
  warmup <- s * intv
  n_sim <- 2 * warmup

  pg <- par.control(pg = "sync", intv = intv, min = 0, max = 1, seed = 42)

  out <- set_pg(warmup, n_sim, s, pg)

  idx <- seq(1, warmup, by = pg$intv)

  # Check that all rows in idx have non-zero values
  expect_true(all(rowSums(out[idx, , drop = FALSE]) > 0))
})

test_that("set_pg ord pattern fills matrix correctly", {

  intv <- sample(10, size = 1)
  s <- sample(10, size = 1)
  warmup <- s * intv
  n_sim <- 2 * warmup
  pg <- par.control(pg = "ord", intv = intv, min = 0, max = 1, seed = 123)

  out <- set_pg(warmup, n_sim, s, pg)

  # ord: at least some cells should be non-zero
  expect_true(any(out != 0))

  # ord: all species introduced
  expect_true(all(colSums(out > 0) == 1))

  # ord: correct interval
  idx <- seq(1, warmup, by = intv)
  expect_true(all(rowSums(out > 0)[idx] == 1))
})

test_that("set_pg ord warns if introduction exceeds warmup", {
  n_sim <- 6
  s <- 6
  warmup <- 3
  pg <- par.control(pg = "ord", intv = 1, min = 0, max = 1, seed = 123)

  expect_warning(
    set_pg(warmup, n_sim, s, pg),
    "Introduction times exceed 'warmup'"
  )
})

test_that("set_pg reproducibility with seed", {
  n_sim <- 6
  s <- 3
  warmup <- 5
  pg <- par.control(pg = "sync", intv = 1, min = 0, max = 1, seed = 123)

  out1 <- set_pg(warmup, n_sim, s, pg)
  out2 <- set_pg(warmup, n_sim, s, pg)

  expect_equal(out1, out2)
})

