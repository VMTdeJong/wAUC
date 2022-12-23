# #' AUC
# #'
# #' Calculate an AUC / c-statistic
# #'
# #' @param y numeric outcome vector: each value 1 or 0
# #' @param p numeric predicted probabilities vector: 0 <= p <= 1
# #' @param na.rm Should NA's be removed?
# #' @param ... passed on.

AUC <- function(y, p, na.rm = TRUE, ...) {
  out <- list(call = match.call(),
              AUC_method = list(...)$AUC_method)
  class(out) <- c("AUC", class(out))
  
  if (!is.null(out$AUC_method) && pmatch(out$AUC_method, "factorial", nomatch = 0))
    out$estimate <- AUC_factorial(y, p, na.rm, ...)
  else
    out$estimate <- AUC_default(y, p, na.rm, ...)

  out$separation <- test_separation(out$estimate, ...)
  out
}

AUC_default <- function(y, p, ...) {
  cats <- sort(unique(y))
  n_cat <- length(cats)
  if (n_cat > 2)
    stop(paste("Number of y categories must be 2. y had ", n_cat, " categories.", sep = ""))
  if (n_cat == 1) {
    warning("AUC is undefined when no pairs can be made due to a lack of a comparison category.")
    return(NaN)
  }

  n0 <- sum(y == cats[2])
  n1 <- length(y) - n0
  r  <- rank(p)
  S0 <- sum(as.numeric(r[y == cats[2]]))

  return((S0 - n0 * (n0 + 1)/2)/(as.numeric(n0) * as.numeric(n1)))
}

AUC_factorial <- function(y, p, na.rm = TRUE, ...){
  cats <- sort(unique(y))
  n_cat <- length(cats)
  if (n_cat > 2)
    stop(paste("Number of y categories must be 2. y had ", n_cat, " categories.", sep = ""))
  if (n_cat == 1) {
    warning("AUC is undefined when no pairs can be made due to a lack of a comparison category.")
    return(NaN)
  }
  score_order <- order(p, decreasing = TRUE)
  y <- as.logical(y[score_order])
  p <- p[score_order]

  pos_p <- p[y]
  neg_p <- p[!y]

  s <- outer(sum(y):1, 1:sum(!y), function(i, j) (1 + sign(pos_p[i] - neg_p[j])))

  mean(s, na.rm = na.rm)/2
}

print.AUC <- function(x, digits = 3, ...) {
  cat("AUC: ", round(x$estimate, digits = digits), ".\n", sep = "")
  print_separation(x$estimate)
}
