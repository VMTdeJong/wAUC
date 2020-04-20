#' Weighted c-statistic
#'
#' Calculate a weighted AUC or c-statistic, exactly or by resampling
#'
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#'
#' @aliases wAUC wAUC_resampling wAUC_exact
#'
#' @usage wAUC(y, p, w, na.rm = TRUE, rescale.w = FALSE, method = "resampling", ...)
#' wAUC_resampling(y, p, w, na.rm = TRUE, rescale.w = FALSE, replace = TRUE,
#'                 size = length(p), I = 1000, level = .95, nonunity = FALSE,
#'                 nonzero = FALSE, ret.resamples = FALSE,
#'                 AUC_method = "default", ...)
#' wAUC_exact(y, p, w, na.rm = TRUE, rescale.w = FALSE, nonunity = FALSE,
#'            nonzero = FALSE,  ...)
#'
#' @param y numeric outcome vector: each value 1 or 0
#' @param p numeric predicted probabilities vector: 0 <= p <= 1
#' @param w numeric weights. w >= for exact, 0 <= w <= 1 for resampling
#' @param na.rm Should NA's be removed?
#' @param rescale.w Should w be rescaled to 0 <= w <= 1?
#' @param method Use the "resampling" or "exact" method for calculation?
#' @param replace Should sample be performed with replacement?
#' @param size If sampling is performed with replacement, how many observations should be sampled per
#' propensity weighted bootstrap?
#' @param I number of resamples
#' @param level Confidence level or percentile for the propensity weighted bootstrap percentiles with
#' replacement.
#' @param ret.resamples Should all the resamples be returned?
#' @param AUC_method default or factorial. default is always faster.
#' @param nonunity Should a correction be applied when the estimate is exactly one?
#' @param nonzero Should a correction be applied when the estimate is exactly zero?
#' @param ... Passed to \code{wAUC_resampling} and \code{wAUC_exact}.
#'
#' @return A list of class wAUC. The exact method contains only a scalar: the AUC / c-statistic. The
#' resampling version contains a list with the statistics and possibly the resamples.
#'
#' @details Perfect (c = 1) and perfectly wrong (c = 0) values may lead to computational issues when
#' (meta-)analyzed. Set nonunity to nonzero to TRUE to avoid this, which then respectively substracts or
#' adds a value equal to one tie when this is the case.
#'
#' For the resampling methods, setting replace = TRUE essentially gives a weighted bootstrap
#' estimate. Hence, a confidence interval is returned. Setting replace = FALSE, yields smaller
#' samples for any w < 1, but returns the exact same same for all w == 1. In the latter case
#' the percentiles have the same value as the point estimate and cannot be interpreted as
#' confidence intervals. Hence, no ci is returned when replace = FALSE.
#'
#' Note that unique_events and unique_nonevents may be switched. No check is performed for this.
#'
#' @export
wAUC <- function(y, p, w, na.rm = TRUE, rescale.w = FALSE, method = "resampling", ...) {
  if (pmatch(method, "resampling", nomatch = 0)) {
    out <- wAUC_resampling(y, p, w, na.rm, ...)
    out$estimate <- out$statistics$median

    if (out$options$replace) {
      out$ci.lb <- out$statistics$pct.lb
      out$ci.ub <- out$statistics$pct.ub
    }
  }
  else
    out <- wAUC_exact(y, p, w, na.rm, rescale.w, ...)
  out$call <- match.call()
  out
}

#' @export
wAUC_resampling <- function(y, p, w, na.rm = TRUE, rescale.w = FALSE, replace = TRUE, size = length(p),
                            I = 1000, level = .95, nonunity = FALSE, nonzero = FALSE,
                            ret.resamples = FALSE, AUC_method = "default", ...) {
  if (rescale.w) {
    w <- w - min(w)
    w <- w/max(w)
  }
  resamples <- rep(NA, I)
  n_unique_obs <- matrix(nrow = I, ncol = 2)

  if (replace) {
    for (i in seq_len(I)) {
      indices <- sample(length(p), size, replace = T, prob = w)
      n_unique_obs[i, ] <- by(indices, y, function(x) length(unique(x)))[c(1,2)]
      resamples[i] <- AUC(y[indices], p[indices], na.rm = na.rm, nonunity = nonunity, nonzero = nonzero,
                          AUC_method = AUC_method)$estimate
    }
    colnames(n_unique_obs) <- names(by(indices, y, function(x) length(unique(x)))[c(1,2)])
  } else {
    for (i in seq_len(I)) {
      indices <- which(as.logical(stats::rbinom(length(y), 1, prob = w)))
      n_unique_obs[i, ] <- by(indices, y[indices], function(x) length(unique(x)))[c(1,2)]
      resamples[i] <- AUC(y[indices], p[indices], na.rm = na.rm, nonunity = nonunity, nonzero = nonzero,
                          AUC_method = AUC_method)$estimate
    }
    colnames(n_unique_obs) <- names(by(indices, y[indices], function(x) length(unique(x)))[c(1,2)])
  }

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
                             replace = replace,
                             size = size,
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
wAUC_exact <- function(y, p, w, na.rm = TRUE, rescale.w = FALSE, nonunity = FALSE, nonzero = FALSE,  ...){
  if (rescale.w) {
    w <- w - min(w)
    w <- w/max(w)
  }
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

  s <- (outer(sum(y):1, 1:sum(!y), function(i, j) (sign(pos_p[i] - neg_p[j]))) + 1)
  W <- outer(pos_w, neg_w)

  if (nonunity && all(s == 0))
    s <- c(s, 1)
  if (nonzero && all(s == 2))
    s <- c(s, 1)

  out$estimate <- (mean(s*W, na.rm = na.rm) / mean(W, na.rm = na.rm))/2
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
#' @S3
print.wAUC_resampling <- function(x, digits = 3, ...) {
  if (x$options$replace)
    print(round(with(x, data.frame(ci.lb = ci.lb, estimate = estimate, ci.ub = ci.ub, row.names = "Weighted AUC")), digits = 3))
  else
    print.wAUC(x, digits, ...)
}
