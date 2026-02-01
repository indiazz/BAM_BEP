MDM_on_scores_chain2 <- function(S_matrix, q = floor(ncol(S)^(1/3)),
                                statistic = "S") {
  S <- S_matrix
  stopifnot(is.matrix(S), is.numeric(S))
  K <- nrow(S)
  T_n <- ncol(S)
  stopifnot(K >= 2L, T_n >= 10L)
  d <- S[1:(K-1), , drop = FALSE] - S[2:K, , drop = FALSE]
  dbar <- rowMeans(d)
  
  G <- function(d, h) {
    SCM <- matrix(0, nrow(d), nrow(d))
    for (t in (h + 1):ncol(d)) {
      SCM <- SCM + (d[, t, drop=FALSE] - dbar) %*% t(d[, t - h, drop=FALSE] - dbar)
    }
    SCM / ncol(d)
  }
  Ohm <- G(d, 0)
  if (q > 0) {
    for (h in 1:q) {
      TEMP <- G(d, h)
      Ohm <- Ohm + TEMP + t(TEMP)
    }
  }
  Waldstat <- as.numeric(T_n * t(dbar) %*% solve(Ohm) %*% dbar)
  
  # finite sample correction
  if (identical(statistic, "Sc")) {
    c_T <- 1 - (1 + 2*q)/T_n + q*(q + 1)/(T_n^2)
    Waldstat <- c_T * Waldstat
  }
  pval <- pchisq(Waldstat, df = nrow(d), lower.tail = FALSE)
  structure(list(
    statistic = Waldstat,
    p.value   = pval,
    parameter = c(df = nrow(d), q = q, T = T_n)
  ), class = "htest")
}


