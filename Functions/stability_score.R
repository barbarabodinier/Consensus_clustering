#' Consensus score
#'
#' Computes the consensus score from the consensus matrix obtained with a given
#' parameter controlling the number of clusters. The score is defined as a ratio
#' of likelihoods under the assumption that all pairs of items have the same
#' probability of belonging to the same cluster for the observed and most stable
#' (binary) consensus matrices.
#'
#' @inheritParams StabilityScore
#' @param coprop consensus matrix.
#' @param nc number of clusters used to construct the consensus matrix
#'   \code{coprop}.
#' @param K number of subsampling iterations used to construct the consensus
#'   matrix \code{coprop}.
#' @param linkage character string indicating the type of linkage used in
#'   hierarchical clustering to define the stable clusters. Possible values
#'   include \code{"complete"}, \code{"single"} and \code{"average"} (see
#'   argument \code{"method"} in \code{\link[stats]{hclust}} for a full list).
#'
#' @return A vector of stability scores obtained with the different thresholds
#'   in selection proportions.
#'
#' @family stability metric functions
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#' @examples
#' # Simulating random co-membership proportions
#' set.seed(1)
#' n <- 10
#' coprop <- round(runif(n = n * (n - 1) / 2), digits = 2)
#' M <- matrix(0, n, n)
#' M[upper.tri(M)] <- coprop
#' M <- M + t(M)
#'
#' # Computing consensus score
#' ConsensusScore(M, nc = 3, K = 100)
#' @export
# # Native
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       tmp=tmp+1
#       loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   n_1=p_unif*N
#   n_3=(1-p_unif)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # Equiproba
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   p_unif=1/(nc^2)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       tmp=tmp+1
#       loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   prop=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=prop*N
#   n_3=(1-prop)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # Equiproba within clusters
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   p_unif=1/(nc^2)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       if (theta[i]==theta[j]){
#         tmp=tmp+1
#         loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#       }
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   prop=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=prop*N
#   n_3=(1-prop)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   # best_score=(n_1*p_1+n_3*p_3)
#   best_score=(n_1*p_1)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # Within clusters
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       if (theta[i]==theta[j]){
#         tmp=tmp+1
#         loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#       }
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   prop=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=prop*N
#   n_3=(1-prop)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   # best_score=(n_1*p_1+n_3*p_3)
#   best_score=(n_1*p_1)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # At least as stable within clusters
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       if (theta[i]==theta[j]){
#         tmp=tmp+1
#         loglik=loglik+pbinom(round(coprop[i,j]*K)-1, size=K, prob=p_unif, lower.tail = FALSE, log = TRUE)
#       }
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   prop=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=prop*N
#   n_3=(1-prop)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   # best_score=(n_1*p_1+n_3*p_3)
#   best_score=(n_1*p_1)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # At least as stable within and between clusters
# ConsensusScore <- function(coprop, nc, K=K, linkage="complete") {
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   # print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   # print(p_unif)
#   # print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   # print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       if (theta[i]==theta[j]){
#         loglik=loglik+pbinom(round(coprop[i,j]*K)-1, size=K, prob=p_unif, lower.tail = FALSE, log = TRUE)
#       } else {
#         loglik=loglik+pbinom(round(coprop[i,j]*K), size=K, prob=p_unif, lower.tail = TRUE, log = TRUE)
#       }
#     }
#   }
#   # print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   n_1=p_unif*N
#   n_3=(1-p_unif)*N
#   # print(n_1)
#   # print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#   # best_score=(n_1*p_1)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# At least as stable within and between clusters and normalised
ConsensusScore <- function(coprop, nc, K = K, linkage = "complete") {
  # Clustering on the consensus matrix
  sh_clust <- hclust(as.dist(1 - coprop), method = linkage)

  # Identifying stable clusters
  theta <- cutree(sh_clust, k = nc)
  # print(table(theta))

  # Probability that items i and j belong to the same cluster
  N <- length(theta) * (length(theta) - 1) / 2
  p_unif <- sum(table(theta) * (table(theta) - 1)) / (2 * N) # P(i) * P(j | i)
  # print(p_unif)
  # print(sum(coprop[upper.tri(coprop)])/N)
  # p_unif=sum(coprop[upper.tri(coprop)])/N
  # print(p_unif)
  # print(N)

  # Calculating log-likelihood for observed consensus matrix
  loglik <- 0
  for (i in 1:(nrow(coprop) - 1)) {
    for (j in (i + 1):nrow(coprop)) {
      if (theta[i] == theta[j]) {
        loglik <- loglik + pbinom(round(coprop[i, j] * K) - 1, size = K, prob = p_unif, lower.tail = FALSE, log = TRUE)
      } else {
        loglik <- loglik + pbinom(round(coprop[i, j] * K), size = K, prob = p_unif, lower.tail = TRUE, log = TRUE)
      }
    }
  }

  # Calculating numbers of within and between cluster pairs
  n_1 <- p_unif * N
  n_3 <- (1 - p_unif) * N

  # Calculating log-likelihood for least stable (uniform) consensus matrix
  p_1 <- pbinom(round(p_unif * K) - 1, size = K, prob = p_unif, lower.tail = FALSE, log = TRUE)
  p_3 <- pbinom(round(p_unif * K), size = K, prob = p_unif, lower.tail = TRUE, log = TRUE)
  worst_score <- (n_1 * p_1 + n_3 * p_3)

  # Calculating log-likelihood for most stable (binary) consensus matrix
  p_1 <- dbinom(K, size = K, prob = p_unif, log = TRUE)
  p_3 <- dbinom(0, size = K, prob = p_unif, log = TRUE)
  best_score <- (n_1 * p_1 + n_3 * p_3)

  # Calculating consensus score as likelihood ratio
  score <- (-loglik + worst_score) / (-best_score + worst_score) # ( x - min ) / ( max - min )

  return(score)
}


# # At least as stable within and between clusters equiproba
# ConsensusScore <- function(coprop, nc, coprobnull=NULL, linkage="complete") {
#   K=length(coprobnull)
#
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   # p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   p_unif=1/(nc^2)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   # p_unif=sum(coprop[upper.tri(coprop)])/N
#   # print(p_unif)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix
#   loglik=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       if (theta[i]==theta[j]){
#         loglik=loglik+pbinom(round(coprop[i,j]*K)-1, size=K, prob=p_unif, lower.tail = FALSE, log = TRUE)
#       } else {
#         loglik=loglik+pbinom(round(coprop[i,j]*K), size=K, prob=p_unif, lower.tail = TRUE, log = TRUE)
#       }
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   # p_unif=min(p_unif, sum(coprop[upper.tri(coprop)])/N)
#   prop=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=prop*N
#   n_3=(1-prop)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#   # best_score=(n_1*p_1)
#   print(loglik)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # Using poisson binomial (problem is not all items are sampled)
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   K=length(coprobnull)
#
#   # # Probability that items i and j belong to the same cluster
#   # N=length(theta)*(length(theta)-1)/2
#   # p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   # print(p_unif)
#   # print(sum(coprop[upper.tri(coprop)])/N)
#   # print(N)
#
#   # Calculating log-likelihood for observed consensus matrix (under expected p under the null)
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       tmp=tmp+1
#       # loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#       loglik=loglik+PoissonBinomial::dpbinom(x=round(coprop[i,j]*K), probs = coprobnull, log = TRUE)
#     }
#   }
#   print(loglik)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix (under the expected p under the null)
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   n_1=p_unif*N
#   n_3=(1-p_unif)*N
#   print(n_1)
#   print(n_3)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#   print(best_score)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# # Using poisson binomial for both (not working)
# ConsensusScore <- function(coprop, nc, coprobnull, linkage="complete") {
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   K=length(coprobnull)
#
#   # # Probability that items i and j belong to the same cluster
#   # N=length(theta)*(length(theta)-1)/2
#   # p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   # print(p_unif)
#   # print(sum(coprop[upper.tri(coprop)])/N)
#   # print(N)
#
#   # Calculating log-likelihood for observed consensus matrix (under expected p under the null)
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       tmp=tmp+1
#       # loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#       loglik=loglik+PoissonBinomial::dpbinom(x=round(coprop[i,j]*K), probs = coprobnull, log = TRUE)
#     }
#   }
#   print(loglik)
#
#   # # Calculating log-likelihood for most stable (binary) consensus matrix (under the expected p under the null)
#   # N=length(theta)*(length(theta)-1)/2
#   # p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   # n_1=p_unif*N
#   # n_3=(1-p_unif)*N
#   # print(n_1)
#   # print(n_3)
#   # p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   # p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   # best_score=(n_1*p_1+n_3*p_3)
#
#   # Expected counts under the most stable scenario
#   coprobnull=sort(coprobnull)
#   n_exp=rep(0, K+1)
#   names(n_exp)=as.character(K:0)
#   N=sum(upper.tri(coprop))
#
#   n_exp[1]=coprobnull[1]*N
#   for (i in 1:(K-1)){
#     n_exp[i+1]=(coprobnull[i+1]-coprobnull[i])*N
#   }
#   n_exp[K+1]=N-sum(n_exp)
#
#   best_score=0
#   for (i in 1:length(n_exp)){
#     best_score=best_score+n_exp[i]*PoissonBinomial::dpbinom(x=as.numeric(names(n_exp[i])), probs = coprobnull, log = TRUE)
#   }
#   print(best_score)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }


# ConsensusScore <- function(coprop, nc, K, linkage="complete") {
#   # Clustering on the consensus matrix
#   sh_clust=hclust(as.dist(1-coprop), method = linkage)
#
#   # Identifying stable clusters
#   theta=cutree(sh_clust, k=nc)
#   print(table(theta))
#
#   # Probability that items i and j belong to the same cluster
#   N=length(theta)*(length(theta)-1)/2
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
#   print(p_unif)
#   print(sum(coprop[upper.tri(coprop)])/N)
#   print(N)
#
#   # Calculating log-likelihood for observed consensus matrix (under expected p under the null)
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # based on expected under the null: P(i) * P(j | i)
#   loglik=0
#   tmp=0
#   for (i in 1:(nrow(coprop)-1)){
#     for (j in (i+1):nrow(coprop)){
#       tmp=tmp+1
#       loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#     }
#   } # XXX maybe a problem with different cluster sizes (ratio above 1)
#
#   # Calculating log-likelihood for most stable (binary) consensus matrix
#   p_unif=sum(coprop[upper.tri(coprop)])/N # based on observed p
#   n_1=p_unif*N
#   n_3=(1-p_unif)*N
#   print(n_1)
#   print(n_3)
#
#   p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # based on expected under the null: P(i) * P(j | i)
#   p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
#   p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
#   best_score=(n_1*p_1+n_3*p_3)
#
#   # Calculating consensus score as likelihood ratio
#   score=loglik/best_score
#
#   return(score)
# }
