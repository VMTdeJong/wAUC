#' Weighted c-statistic
#'
#' Calculate a weighted AUC or c-statistic, exact or by resampling
#'
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#'
#' @aliases wAUC wAUC_resampling wAUC_exact
#'
#' @usage wAUC(y, p, w, na.rm = TRUE, method = "resampling", ...)
#' wAUC_resampling(y, p, w, na.rm = TRUE, I = 1000, level = .95, nonunity = FALSE,
#'                 nonzero = FALSE, ret.resamples = FALSE, ...)
#' wAUC_exact(y, p, w, na.rm = TRUE, nonunity = FALSE, nonzero = FALSE,  ...)
#'
#' @param y numeric outcome vector: each value 1 or 0
#' @param p numeric predicted probabilities vector: 0 <= p <= 1
#' @param w numeric weights. w >= 0 for exact, 0 <= w <= 1 for resampling. Set to 
#' \code{rep(1,length(y))} to obtain the unweighted AUC / c-statistic.
#' @param na.rm Should NA's be removed?
#' @param method Use the "resampling" or "exact" method for calculation?
#' @param I number of bootstrap resamples
#' @param level Confidence level or quantile for the propensity weighted bootstrap.
#' @param ret.resamples Should all the resamples be returned?
#' @param nonunity Should a correction be applied when the estimate is exactly one?
#' @param nonzero Should a correction be applied when the estimate is exactly zero?
#' @param ... Passed to \code{wAUC_resampling} and \code{wAUC_exact}.
#'
#' @return A list of class wAUC. The exact method contains only a scalar: the AUC / c-statistic, as 
#' well as a list of the options. The resampling method contains a list with the statistics and 
#' possibly the resamples, as well as the options. With the resampling method, the point estimate
#' and 95\% CI are given by the .025, .975 and .5 quantiles of the bootstrap, respectively.
#'
#' @details Perfect (c = 1) and perfectly wrong (c = 0) values may lead to computational issues when
#' (meta-)analyzed. Set nonunity and/or nonzero to TRUE to avoid this, which then respectively substracts or
#' adds a value equal to one tie when this is the case.
#'
#' Note that unique_events and unique_nonevents may be switched. No check is performed for this.
#'
#' @export
wAUC <- function(y, p, w, na.rm = TRUE, method = "resampling", ...) {
  if (pmatch(method, "resampling", nomatch = 0)) {
    out <- wAUC_resampling(y, p, w, na.rm, ...)
    
    out$estimate <- out$statistics$median
    out$ci.lb <- out$statistics$pct.lb
    out$ci.ub <- out$statistics$pct.ub
  }
  else
    out <- wAUC_exact(y, p, w, na.rm, ...)
  out$call <- match.call()
  out
}

#' @export
wAUC_resampling <- function(y, p, w, na.rm = TRUE, I = 1000, level = .95, 
                            nonunity = FALSE, nonzero = FALSE,
                            ret.resamples = FALSE, ...) {
  resamples <- rep(NA, I)
  n_unique_obs <- matrix(nrow = I, ncol = 2)
  
  for (i in seq_len(I)) {
    indices <- sample(x = length(p), size = length(p), replace = T, prob = w)
    n_unique_obs[i, ] <- by(indices, y, function(x) length(unique(x)))[c(1,2)]
    resamples[i] <- AUC(y[indices], p[indices], na.rm = na.rm, nonunity = nonunity, 
                        nonzero = nonzero, ...)$estimate
  }
  colnames(n_unique_obs) <- names(by(indices, y, function(x) length(unique(x)))[c(1,2)])
  
  out <- list(statistics = data.frame(median = stats::median(resamples, na.rm = na.rm),
                                      mean = mean(resamples, na.rm = na.rm),
                                      pct.lb = stats::quantile(resamples, probs = (1 - level)/2),
                                      pct.ub = stats::quantile(resamples, probs = 1 - (1 - level)/2),
                                      se = stats::sd(resamples, na.rm = na.rm),
                                      iqr = stats::IQR(resamples, na.rm = na.rm),
                                      mad = stats::mad(resamples),
                                      unique_events = mean(n_unique_obs[ , 2]),
                                      unique_nonevents = mean(n_unique_obs[ , 1]),
                                      row.names = "Estimate"),
              options = list(na.rm = na.rm,
                             I = I,
                             level = level,
                             nonunity = nonunity,
                             nonzero = nonzero,
                             ret.resamples = ret.resamples)
  )
  if (ret.resamples) {
    out$resamples <- resamples
    out$n_unique_obs <- n_unique_obs
  }
  class(out) <- c("wAUC_resampling", "wAUC", class(out))
  out
}

#' @export
wAUC_exact <- function(y, p, w, na.rm = TRUE, nonunity = FALSE, nonzero = FALSE,  ...){
  out <- list(options = list(na.rm = na.rm,
                              nonunity = nonunity,
                              nonzero = nonzero))
  class(out) <- c("wAUC_exact", "wAUC", class(out))

  score_order <- order(p, decreasing=TRUE)
  y <- as.logical(y[score_order])
  p <- p[score_order]
  w <- w[score_order]

  pos_p <- p[y]
  neg_p <- p[!y]
  pos_w <- w[y]
  neg_w <- w[!y]

  # + 1 converts the -1 or +1 from sign() to 0 or 2, 
  # then divide final estimate by 2 to obtain 0 or 1.
  s <- (outer(1:sum(y), 1:sum(!y), function(i, j) (sign(pos_p[i] - neg_p[j]))) + 1)
  W <- outer(pos_w, neg_w)

  if (nonunity && all(s == 0))
    s <- c(s, 1)
  if (nonzero && all(s == 2))
    s <- c(s, 1)

  out$estimate <- (mean(s*W, na.rm = na.rm) / mean(W, na.rm = na.rm))/2
  out
}

# This function is always slower than the one above, and should always produce
# the exact same value. The only reason that it is included is because this uses
# the same notation as we did in the manuscript.
wAUC_exact_slow <- function(y, p, w, na.rm = TRUE, nonunity = FALSE, nonzero = FALSE,  ...){
  out <- list(options = list(na.rm = na.rm,
                             nonunity = nonunity,
                             nonzero = nonzero))
  class(out) <- c("wAUC_exact", "wAUC", class(out))
  
  score_order <- order(p, decreasing=TRUE)
  y <- as.logical(y[score_order])
  p <- p[score_order]
  w <- w[score_order]
  
  pos_p <- p[y]
  neg_p <- p[!y]
  pos_w <- w[y]
  neg_w <- w[!y]
  
  cs <- matrix(nrow = length(pos_p), ncol = length(neg_p))
  W <- outer(pos_w, neg_w)

  for (i in seq_along(pos_p))
    for (q in seq_along(neg_p))
      cs[i, q] <- (pos_p[i] > neg_p[q]) * W[i, q]
  
  if (any(ties <- pos_p %in% neg_p))  # save ties for pos_p, if any ties, assess which
    for (tie in which(ties))          # then only visit the pos_p that are tied (with neg_p)
      for (q in seq_along(neg_p))     # But visit all neg_p, because we don't know which neg_p were tied yet.
        if (pos_p[tie] == neg_p[q])
          cs[tie, q] <- W[tie, q] / 2 # weight * 1/2 for a tie
  
  out$estimate <- sum(cs, na.rm = na.rm) / sum(W)
  
  out
}

#' @author Valentijn de Jong
#' @method print wAUC
#' @export
print.wAUC <- function(x, digits = 3, ...) {
  cat("Weighted AUC: ", round(x$estimate, digits = digits), ".\n", sep = "")
}
#' @author Valentijn de Jong
#' @method print wAUC_resampling
#' @export
print.wAUC_resampling <- function(x, digits = 3, ...) {
    print(round(with(x, data.frame(ci.lb = ci.lb, estimate = estimate, ci.ub = ci.ub, row.names = "Weighted AUC")), digits = 3))
}
