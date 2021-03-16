context("Weighted AUC")

### Parameters
set.seed(0409)
n <- 200
w_min <- .1
w_max <- .4

### Data
p <- runif(n)
y <- rbinom(n, 1, p)
w <- runif(n, min = w_min, max = w_max)

### Tests
test_that("For random weights, the exact and resampling methods yield nearly equal point estimates", {
  w_auc_est_exact   <- wAUC(y, p, w, method = "exact")
  w_auc_est_replace <- wAUC(y, p, w, method = "resampling", replace = TRUE)

  est_replace <- w_auc_est_replace$estimate
  est_exact <- w_auc_est_exact$estimate

  expect_lt(est_exact - est_replace, .005)
})

test_that("For random weights, the two exact methods yield identical point estimates", {
  w_auc_est_exact      <- wAUC(y, p, w, method = "exact")
  w_auc_est_exact_slow <- wAUC_exact_slow(y, p, w)
  
  est_exact      <- w_auc_est_exact$estimate
  est_exact_slow <- w_auc_est_exact_slow$estimate
  
  expect_identical(est_exact, est_exact_slow)
})


# library("microbenchmark")
# 
# n <- 100
# p <- runif(n)
# y <- rbinom(n, 1, p)
# w <- runif(n, w_min, w_max)
# I <- 250
# 
# microbenchmark(wAUC:::wAUC_exact(y, p, w),
#                wAUC:::wAUC(y, p, w, method = "exact"),
#                wAUC:::wAUC_exact_slow(y, p, w),
#                wAUC:::wAUC(y, p, w, method = "resampling", I = I),
#                wAUC:::wAUC(y, p, w, method = "resampling", I = I, AUC_method = "AUC_factorial")
# )
