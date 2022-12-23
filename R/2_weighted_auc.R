#' Weighted c-statistic
#'
#' Calculate a weighted AUC or c-statistic, exact or by resampling
#'
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#'
#' @aliases wAUC wAUC_resample wAUC_exact
#'
#' @usage wAUC(y, p, w, na.rm = TRUE, method = "resample", ...)
#' wAUC_resample(y, p, w, na.rm = TRUE, I = 1000, level = .95, 
#' ret.resamples = FALSE, ...)
#' wAUC_exact(y, p, w, na.rm = TRUE, ...)
#'
#' @param y numeric outcome vector: each value 1 or 0
#' @param p numeric predicted probabilities vector: 0 <= p <= 1
#' @param w numeric weights. w >= 0 for exact, 0 <= w <= 1 for resample. Set to 
#' \code{NULL} or \code{rep(1,length(y))} to obtain the unweighted AUC / 
#' c-statistic.
#' @param na.rm Should NA's be removed?
#' @param method Use the "resample" or "exact" method for calculation?
#' @param I number of bootstrap resamples
#' @param level Confidence level or quantile for the propensity weighted 
#' bootstrap.
#' @param ret.resamples Should all the resamples be returned?
#' @param ... Passed to \code{wAUC_resample} and \code{wAUC_exact}.
#'
#' @return A list of class wAUC. The exact method contains a scalar: the 
#' AUC / c-statistic, as well as a list of the options. The resample method 
#' contains a list with the statistics and possibly the resamples, as well as 
#' the options. With the resample method, the point estimate and 95\% CI are 
#' given by the .025, .975 and .5 quantiles of the bootstrap, respectively.
#'
#' @details The confidence intervals produced by the resampling method assume 
#' that a predefined model is validated/tested. If model selection is performed, 
#' then the model selection should be performed within the bootstrap resampling 
#' procedure to account for the uncertainty by the model selection.  
#'
#' Note that unique_events and unique_nonevents may be switched. No check is 
#' performed for this.
#'
#' @export
#' @examples
#' # Fictional data
#' y <- rep(c(rep(0, 10), rep(1, 10)),5)
#' p <- rep(c(0:9, 1:10)/10, 5)
#' w <- (1:100)/100
#' 
#' # Unweighted AUC / concordance statistic
#' wAUC(y, p, NULL)
#' wAUC(y, p, w = rep(1,length(y)))
#' 
#' # Weighted AUC / concordance statistic
#' # With these 'random' weights the result is similar to the unweighted version
#' wAUC(y, p, w)
#' wAUC(y, p, w, method = "exact")
#' wAUC(y, p, w, method = "resample") # Currently the default
#' 
#' # With these obviously non-random weights, we get a very different result
#' w2 <- (p + y)/2
#' wAUC(y, p, w2, method = "exact")
wAUC <- function(y, p, w, na.rm = TRUE, method = "resample", ...) {
  if (pmatch(method, "resample", nomatch = 0))
    out <- wAUC_resample(y, p, w, na.rm, ...)
  else
    out <- wAUC_exact(y, p, w, na.rm, AUC_method = method, ...)
  out$call <- match.call()
  out
}

#' @export
wAUC_resample <- function(y, p, w, na.rm = TRUE, I = 1000, level = .95, 
                            ret.resamples = FALSE, ...) {
  if (is.null(w)) {
    weighted <- FALSE
    w <- rep(1, length(y))
  } else weighted <- TRUE
  
  resamples <- rep(NA, I)
  n_unique_obs <- matrix(nrow = I, ncol = 2)
  
  for (i in seq_len(I)) {
    indices <- sample(x = length(p), size = length(p), replace = T, prob = w)
    n_unique_obs[i, ] <- by(indices, y, function(x) length(unique(x)))[c(1,2)]
    resamples[i] <- AUC(y[indices], p[indices], na.rm = na.rm, warn_separation = FALSE, ...)$estimate
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
                             weighted = weighted,
                             ret.resamples = ret.resamples)
  )
  out$estimate <- out$statistics$median
  out$ci.lb <- out$statistics$pct.lb
  out$ci.ub <- out$statistics$pct.ub
  
  if (ret.resamples) {
    out$resamples <- resamples
    out$n_unique_obs <- n_unique_obs
  }
  class(out) <- c("wAUC_resample", "wAUC", class(out))
  out$separation <- test_separation(out$statistics$median, ...)
  out
}

#' @export
wAUC_exact <- function(y, p, w, na.rm = TRUE, ...){
  if (is.null(w)) {
    weighted <- FALSE
    w <- rep(1, length(y))
  } else weighted <- TRUE
  
  out <- list(options = list(na.rm = na.rm,
                             weighted = weighted))
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

  out$estimate <- (mean(s*W, na.rm = na.rm) / mean(W, na.rm = na.rm))/2
  out$separation <- test_separation(out$estimate, ...)
  out
}

# This function is always slower than the one above, and should always produce
# the (numerically) exact same value. The only reason that it is included is 
# because this uses the same notation as we did in the manuscript.
wAUC_exact_slow <- function(y, p, w, na.rm = TRUE, ...){
  if (is.null(w)) {
    weighted <- FALSE
    w <- rep(1, length(y))
  } else weighted <- TRUE
  
  out <- list(options = list(na.rm = na.rm,
                             weighted = weighted))
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
  
  if (any(ties <- pos_p %in% neg_p))  # Save ties for pos_p, if any ties, assess which
    for (tie in which(ties))          # then only visit the pos_p that are tied (with neg_p).
      for (q in seq_along(neg_p))     # But visit all neg_p, because we don't know which neg_p were tied yet.
        if (pos_p[tie] == neg_p[q])
          cs[tie, q] <- W[tie, q] / 2 # weight
  
  out$estimate <- sum(cs, na.rm = na.rm) / sum(W)
  out$separation <- test_separation(out$estimate, ...)
  out
}

#' @author Valentijn de Jong
#' @method print wAUC
#' @export
print.wAUC <- function(x, digits = 3, ...) {
  name <- if (x$options$weighted) "Weighted AUC" else "AUC"
  cat(name, ": ", round(x$estimate, digits = digits), "\n", sep = "")
  print_separation(x$estimate)
  invisible(x)
}
#' @author Valentijn de Jong
#' @method print wAUC_resample
#' @export
print.wAUC_resample <- function(x, digits = 3, ...) {
  name <- if (x$options$weighted) "Weighted AUC" else "AUC"
  print(round(with(x, data.frame(ci.lb = ci.lb, estimate = estimate, ci.ub = ci.ub, row.names = name)), digits = 3))
  print_separation(x$estimate)  
  invisible(x)
}

test_separation <- function(x, warn_separation = TRUE, ...) {
  if (identical(x, 1)) {
    if (warn_separation)
      warning("The estimated (weighted) AUC was exactly equal to 1. The point estimate and the confidence interval can be highly misleading.")
    return(TRUE)
  }
  if (identical(x, 0)) {
    if (warn_separation)
      warning("The estimated (weighted) AUC was exactly equal to 0. The point estimate and the confidence interval can be highly misleading.")
    return(TRUE)
  }  
  return(FALSE)
}

print_separation <- function(x, ...) {
  if (identical(x, 1))
    cat("The estimated (weighted) AUC was exactly equal to 1. The point estimate and the confidence interval can be highly misleading.\n")
  if (identical(x, 0))
    cat("The estimated (weighted) AUC was exactly equal to 0. The point estimate and the confidence interval can be highly misleading.\n")
}
