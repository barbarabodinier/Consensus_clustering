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
#' n=10
#' coprop <- round(runif(n = n*(n-1)/2), digits = 2)
#' M=matrix(0, n, n)
#' M[upper.tri(M)]=coprop
#' M=M+t(M)
#' 
#' # Computing consensus score
#' ConsensusScore(M, nc=3, K=100)
#'
#' @export
# # Native
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


# Using min/max for the ratio
ConsensusScore <- function(coprop, nc, K, linkage="complete") {
  # Clustering on the consensus matrix
  sh_clust=hclust(as.dist(1-coprop), method = linkage)

  # Identifying stable clusters
  theta=cutree(sh_clust, k=nc)
  print(table(theta))

  # Probability that items i and j belong to the same cluster
  N=length(theta)*(length(theta)-1)/2
  p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
  print(p_unif)
  print(sum(coprop[upper.tri(coprop)])/N)
  print(N)

  # Calculating log-likelihood for observed consensus matrix (under expected p under the null)
  p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
  loglik=0
  tmp=0
  for (i in 1:(nrow(coprop)-1)){
    for (j in (i+1):nrow(coprop)){
      tmp=tmp+1
      loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
    }
  } # XXX maybe a problem with different cluster sizes (ratio above 1)
  loglik1=loglik
  print("loglik1")
  print(loglik1)

  # Calculating log-likelihood for observed consensus matrix (under observed p)
  p_unif=sum(coprop[upper.tri(coprop)])/N
  loglik=0
  tmp=0
  for (i in 1:(nrow(coprop)-1)){
    for (j in (i+1):nrow(coprop)){
      tmp=tmp+1
      loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
    }
  } # XXX maybe a problem with different cluster sizes (ratio above 1)
  loglik2=loglik
  print("loglik2")
  print(loglik2)

  # Calculating log-likelihood for most stable (binary) consensus matrix (under the expected p under the null)
  p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
  n_1=p_unif*N
  n_3=(1-p_unif)*N
  print(n_1)
  print(n_3)
  p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
  p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
  best_score=(n_1*p_1+n_3*p_3)
  best_score1=best_score
  print("best_score1")
  print(best_score1)

  # Calculating log-likelihood for most stable (binary) consensus matrix (under the observed p)
  p_unif=sum(coprop[upper.tri(coprop)])/N
  n_1=p_unif*N
  n_3=(1-p_unif)*N
  print(n_1)
  print(n_3)
  p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
  p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
  best_score=(n_1*p_1+n_3*p_3)
  best_score2=best_score
  print("best_score2")
  print(best_score2)

  # Calculating consensus score as likelihood ratio
  score=loglik1/min(best_score1, best_score2)

  return(score)
}


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


StabilityRatio <- function(selprop, K, group = NULL) {
  # Preparing objects
  if (is.matrix(selprop)) {
    selprop <- selprop[upper.tri(selprop)]
  }
  
  # Using group penalisation (extracting one per group)
  if (!is.null(group)) {
    selprop <- selprop[cumsum(group)]
  }
  
  # Computing the number of features (edges/variables)
  N <- length(selprop)
  
  # Computing the average number of selected features
  q <- round(sum(selprop, na.rm = TRUE))
  print(q)
  print(N)
  p_unif=q/N
  
  # Calculating log-likelihood for observed selection proportions
  loglik=0
  for (i in 1:length(selprop)){
    loglik=loglik+dbinom(round(selprop[i]*K), size=K, prob=p_unif, log = TRUE)
  }
  
  # Calculating log-likelihood for most stable (binary) selection proportions
  n_1=p_unif*N
  n_3=(1-p_unif)*N
  p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
  p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
  best_score=(n_1*p_1+n_3*p_3)
  
  # Calculating consensus score as likelihood ratio
  score=loglik/best_score
  
  return(score)
}


