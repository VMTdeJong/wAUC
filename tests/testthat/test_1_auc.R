context("Unweighted AUC")

### Parameters
set.seed(0409)
n <- 200

### Data
p <- runif(n)
y <- rbinom(n, 1, p)
w <- rep(1, n)

### Tests
test_that("The exact AUC methods exactly equal the AUC", {
  proc_est      <- unlist(pROC::auc(y, p)[[1]])
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

test_that("For uniform weights, the exact wAUC methods exactly equal the default AUC", {
  proc_est <- unlist(pROC::auc(y, p)[[1]])
  weighted_exact_auc_est <- wAUC_exact(y, p, w)$est
  weighted_e_auc_est <- wAUC(y, p, w, method = "exact")$est
  
  expect_identical(proc_est, weighted_exact_auc_est)
  expect_identical(weighted_exact_auc_est, weighted_e_auc_est)
})

test_that("For uniform weights, the resampling methods nearly equal default c-statistic and its ci", {
  proc_fit <- pROC::auc(y, p)
  proc_ci  <- pROC::ci(proc_fit)[1:3]
  
  resampling_fit <- wAUC(y, p, w, method = "resampling", replace = TRUE, I = 5000)
  resampling_ci  <- unlist(resampling_fit[c("ci.lb", "estimate", "ci.ub")])
  
  diff <- proc_ci - resampling_ci
  expect_lt(diff[[1]], .005)
  expect_lt(diff[[2]], .005)
  expect_lt(diff[[3]], .005)
})




# library("microbenchmark")
# 
# n <- 20
# p <- runif(n)
# y <- rbinom(n, 1, p)
# 
# microbenchmark(AUC_default(y, p), AUC_factorial(y,p))
# 
# n <- 200
# p <- runif(n)
# y <- rbinom(n, 1, p)
# 
# microbenchmark(AUC_default(y, p), AUC_factorial(y,p))
# 
# n <- 2000
# p <- runif(n)
# y <- rbinom(n, 1, p)
# 
# microbenchmark(AUC_default(y, p), AUC_factorial(y,p))