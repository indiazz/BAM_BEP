selection <- function (realized = NULL, evaluated = NULL, q, alpha, statistic = "Sc", S_matrix = NULL) 
{
  .in.MDM.test <- getFromNamespace(".in.MDM.test", "multDM")
  p <- 0
  S <- S_matrix
  models <- (1:nrow(S))
  if (is.null(rownames(S))) {
    rownames(S) <- seq(1:nrow(S))
  }
  numb_models <- nrow(S)
  while (p < alpha && length(models) > 1) {
    d <- S
    d <- d_t(d)
    mdm <- .in.MDM.test(d = d, q = q, statistic = statistic)
    p <- mdm[[2]]$p.value
    if (p < alpha) {
      z <- which.max(abs(mdm[[1]]))
      if (mdm[[1]][z] > 0) {
        models <- models[-z]
        z.drop <- z
      }
      else {
        models <- models[-(z + 1)]
        z.drop <- z + 1
      }
      S <- S[-z.drop, , drop = FALSE]
      models <- (1:nrow(S))
    }
  }
  numb_models <- numb_models - nrow(S)
  result <- matrix(NA, ncol = 3, nrow = nrow(S))
  rownames(result) <- rownames(S)
  result[-nrow(result), 2] <- mdm[[1]]
  result[-nrow(result), 1] <- rank(abs(result[-nrow(result), 2]))
  result[, 3] <- rowMeans(S)
  colnames(result) <- c("rank", "s", "mean loss")
  result <- list(result, as.numeric(p), alpha, numb_models)
  names(result) <- c("outcomes", "p.value", "alpha", "eliminated")
  class(result) <- "MDM"
  return(result)
}
