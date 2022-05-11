#' Clustering performance
#'
#' Computes different metrics of clustering performance by comparing true and
#' predicted pairwise co-membership. This function can only be used in
#' simulation studies (i.e. when the true cluster membership is known).
#'
#' @param theta vector of group membership, binary and symmetric matrix of
#'   co-membership, or output of \code{\link{GraphicalModel}}.
#' @param theta_star vector of group membership, or binary and symmetric matrix
#'   of co-membership.
#' @param ... additional arguments to be passed to \code{\link{Clusters}}.
#'
#' @return A matrix of selection metrics including:
#'
#'   \item{TP}{number of True Positives (TP)} \item{FN}{number of False
#'   Negatives (TN)} \item{FP}{number of False Positives (FP)} \item{TN}{number
#'   of True Negatives (TN)} \item{sensitivity}{sensitivity, i.e. TP/(TP+FN)}
#'   \item{specificity}{specificity, i.e. TN/(TN+FP)} \item{accuracy}{accuracy,
#'   i.e. (TP+TN)/(TP+TN+FP+FN)} \item{precision}{precision (p), i.e.
#'   TP/(TP+FP)} \item{recall}{recall (r), i.e. TP/(TP+FN)}
#'   \item{F1_score}{F1-score, i.e. 2*p*r/(p+r)} \item{rand}{Rand Index, i.e.
#'   (TP+TN)/(TP+FP+TN+FN)} \item{ari}{Adjusted Rand Index (ARI), i.e.
#'   2*(TP*TN-FP*FN)/((TP+FP)*(TN+FP)+(TP+FN)*(TN+FN))} \item{jaccard}{Jaccard
#'   index, i.e. TP/(TP+FP+FN)} \item{ami}{Adjusted Mutual Information (see
#'   \code{\link[aricode]{AMI}})}
#'
#' @family functions for model performance
#'
#' @examples
#' \dontrun{
#' # Simulation with 5 groups of correlated variables
#' set.seed(1)
#' pk <- sample(1:10, size = 5, replace = TRUE)
#' print(pk)
#' simul <- SimulateGraphical(
#'   n = 100, pk = pk,
#'   nu_within = 0.6,
#'   nu_between = 0.1,
#'   v_within = c(0.1, 1),
#'   v_between = c(0, 0.3),
#'   v_sign = -1
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(simul$data),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # Stability grouping
#' stab <- GraphicalModel(
#'   xdata = simul$data,
#'   Lambda = 1:ncol(simul$data),
#'   implementation = HierarchicalClustering
#' )
#'
#' # Clustering performance
#' ClusteringPerformance(
#'   theta = stab,
#'   theta_star = Clusters(BlockDiagonal(pk))
#' )
#' ClusteringPerformance(
#'   theta = Clusters(stab),
#'   theta_star = Clusters(BlockDiagonal(pk))
#' ) # alternative formulation with membership
#' ClusteringPerformance(
#'   theta = CoMembership(Clusters(stab)),
#'   theta_star = Clusters(BlockDiagonal(pk))
#' ) # alternative formulation with co-membership
#' }
#'
#' @export
ClusteringPerformance <- function(theta, theta_star, ...) {
  # Re-formatting input theta
  if (any(class(theta) %in% c("graphical_model", "clustering"))) {
    theta <- Clusters(theta, ...)
  }

  # Re-formatting input theta_star
  if (any(class(theta_star) %in% c("simulation_clustering"))) {
    theta_star <- theta_star$theta
  }

  # Initialising unused parameters
  cor <- NULL
  thr <- 0.5

  # Computing co-membership matrices
  if (!is.matrix(theta)) {
    theta <- CoMembership(theta)
  }
  if (!is.matrix(theta_star)) {
    theta_star <- CoMembership(theta_star)
  }

  # Storing similarities/differences between estimated and true sets
  Asum <- theta + 2 * theta_star

  tmp <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
  rand <- (tmp$TP + tmp$TN) / (tmp$TP + tmp$FP + tmp$TN + tmp$FN)
  ari <- 2 * ((tmp$TP * tmp$TN) - (tmp$FP * tmp$FN)) / ((tmp$TP + tmp$FP) * (tmp$TN + tmp$FP) + (tmp$TP + tmp$FN) * (tmp$TN + tmp$FN))
  jaccard <- (tmp$TP) / (tmp$TP + tmp$FP + tmp$FN)
  ami <- aricode::AMI(c1 = as.vector(theta), c2 = as.vector(theta_star))
  tmp <- cbind(tmp, rand = rand, ari = ari, jaccard = jaccard, ami = ami)
  return(tmp)
}
