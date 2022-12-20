context("Unweighted AUC")

### Parameters
set.seed(0409)
n <- 200

### Data
p <- runif(n)
y <- rbinom(n, 1, p)
w <- rep(1, n)

### Tests
test_that("The unweighted exact AUC methods exactly equal the AUC", {
  proc_est      <- suppressMessages(unlist(pROC::auc(y, p)[[1]]))
  auc_est       <- AUC(y, p)$est
  auc_f_est     <- AUC(y, p, AUC_method = "factorial")$est
  auc_d_est     <- AUC(y, p, AUC_method = "default")$est
  factorial_est <- AUC_factorial(y, p)
  default_est   <- AUC_default(y,p)
  
  expect_identical(proc_est, auc_est)
  expect_identical(auc_est, auc_f_est)
  expect_identical(auc_f_est, auc_d_est)
  expect_identical(auc_d_est, factorial_est)
  expect_identical(factorial_est, default_est)
})

test_that("For unity weights, the exact wAUC methods exactly equal the default AUC", {
  proc_est <- suppressMessages(unlist(pROC::auc(y, p)[[1]]))
  weighted_exact_auc_est <- wAUC_exact(y, p, w)$est
  weighted_e_auc_est <- wAUC(y, p, w, method = "exact")$est
  
  expect_identical(proc_est, weighted_exact_auc_est)
  expect_identical(weighted_exact_auc_est, weighted_e_auc_est)
})

test_that("For unity weights, the resampling methods nearly equal default c-statistic and its ci", {
  proc_fit <- suppressMessages(pROC::auc(y, p))
  proc_ci  <- suppressMessages(pROC::ci(proc_fit)[1:3])
  
  resampling_fit <- wAUC(y, p, w, method = "re", I = 5000)
  resampling_ci  <- unlist(resampling_fit[c("ci.lb", "estimate", "ci.ub")])
  
  diff <- proc_ci - resampling_ci
  expect_lt(diff[[1]], .005)
  expect_lt(diff[[2]], .005)
  expect_lt(diff[[3]], .005)
})

# Note that the difference may actually be ~1e-16 for some samples with a different seed.
test_that("For unity weights, the exact methods yield identical point estimates", {
  w_auc_est_exact      <- wAUC(y, p, w, method = "exact")
  w_auc_est_exact_slow <- wAUC_exact_slow(y, p, w)
  
  est_exact      <- w_auc_est_exact$estimate
  est_exact_slow <- w_auc_est_exact_slow$estimate
  
  expect_identical(est_exact, est_exact_slow)
})

test_that("wAUC can be printed", {
  expect_true(inherits(print_re <- capture.output(print(wAUC(y, p, w, method = "re", I = 5))), "character"))
  expect_true(inherits(print_ex <- capture.output(print(wAUC(y, p, w, method = "exact", I = 5))), "character"))
  expect_true(inherits(print__r <- capture.output(print(wAUC_resample(y, p, w, I = 5))), "character"))
})


# library("microbenchmark")
# 
# n <- 20
# p <- runif(n)
# y <- rbinom(n, 1, p)
# 
# microbenchmark(wAUC:::AUC_default(y, p), 
#                wAUC:::AUC(y, p),
#                wAUC:::AUC_factorial(y,p), 
#                suppressMessages(pROC::auc(y, p)))
# 
# n <- 200
# p <- runif(n)
# y <- rbinom(n, 1, p)
# 
# microbenchmark(wAUC:::AUC_default(y, p), 
#                wAUC:::AUC(y, p),
#                wAUC:::AUC_factorial(y,p),
#                suppressMessages(pROC::auc(y, p)))
# AUC_default is always faster

