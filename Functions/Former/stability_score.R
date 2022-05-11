#' Stability in clustering score
#'
#' Computes the stability in clustering score from co-membership proportions of
#' models with a given parameter controlling the number of clusters and for
#' different thresholds in co-membership proportions. The score measures how
#' unlikely it is that the clustering procedure is uniform (i.e. uninformative)
#' for a given combination of parameters.
#'
#' @inheritParams StabilityMetrics
#'
#' @details XXX
#'
#' @return A vector of stability scores obtained with the different thresholds
#'   in selection proportions.
#'
#' @family stability metric functions
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#' @examples
#' # Simulating set of selection proportions
#' selprop <- round(runif(n = 20), digits = 2)
#'
#' # Computing stability scores for different thresholds
#' score <- StabilityScore(selprop, pi_list = c(0.6, 0.7, 0.8), K = 100)
#' @export
StabilityScoreClustering <- function(selprop, nc, pi_list = seq(0.6, 0.9, by = 0.01), K) {
  # Preparing objects
  if (is.matrix(selprop)) {
    M <- selprop
    selprop <- selprop[upper.tri(selprop)]
  } else {
    M <- matrix(0, nrow = n, ncol = n)
    M[upper.tri(M)] <- selprop
    M <- M + t(M)
  }

  # Clustering on the consensus matrix
  sh_clust <- hclust(as.dist(1 - M), method = "complete")

  # Identifying stable clusters
  theta <- cutree(sh_clust, k = nc)
  print(table(theta))
  pi_list <- sh_clust$height[length(sh_clust$height) - (nc - 2)]
  print(pi_list)
  cat("\n")

  # Probability that items i and j belong to the same cluster
  p_unif <- sum(table(theta) * (table(theta) - 1)) / (nrow(M) * (nrow(M) - 1)) # P(i) * P(j | i)
  print(p_unif)
  N <- nrow(M) * (nrow(M) - 1) / 2

  loglik <- 0
  for (i in 1:(nrow(M) - 1)) {
    for (j in (i + 1):nrow(M)) {
      loglik <- loglik + dbinom(round(M[i, j] * K), size = K, prob = p_unif, log = TRUE)
    }
  }

  n_1 <- p_unif * N
  n_3 <- (1 - p_unif) * N
  p_1 <- dbinom(K, size = K, prob = p_unif, log = TRUE)
  p_3 <- dbinom(0, size = K, prob = p_unif, log = TRUE)
  best_score <- -(n_1 * p_1 + n_3 * p_3)

  score <- -loglik / best_score

  return(score)
}
