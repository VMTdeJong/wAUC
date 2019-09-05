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
test_that("For random weights, the exact and resampling method yield nearly equal point estimates", {
  w_auc_est_exact     <- wAUC(y, p, w, method = "exact")
  w_auc_est_replace   <- wAUC(y, p, w, method = "resampling", replace = TRUE)
  w_auc_est_noreplace <- wAUC(y, p, w, method = "resampling", replace = FALSE)
  
  est_replace   <- w_auc_est_replace$estimate
  est_noreplace <- w_auc_est_noreplace$estimate
  est_exact     <- w_auc_est_replace$estimate
 
  expect_lt(est_replace - est_noreplace, .005)
  expect_lt(est_noreplace - est_exact  , .005)
  expect_lt(est_exact - est_replace    , .005)
})


# library("microbenchmark")
# 
# n <- 100
# p <- runif(n)
# y <- rbinom(n, 1, p)
# w <- runif(n, w_min, w_max)
# I <- 250
# 
# microbenchmark(AUC_default(y, p), 
#                AUC_factorial(y, p),
#                wAUC_exact(y, p, w),
#                wAUC(y, p, w, method = "exact"))
# 
# microbenchmark(wAUC(y, p, w, method = "exact"),
#                wAUC(y, p, w, method = "resampling", I = I, replace = TRUE),
#                wAUC(y, p, w, method = "resampling", I = I, replace = FALSE),
#                wAUC(y, p, w, method = "resampling", I = I, replace = TRUE,  AUC_method = "AUC_factorial"),
#                wAUC(y, p, w, method = "resampling", I = I, replace = FALSE, AUC_method = "AUC_factorial"))
