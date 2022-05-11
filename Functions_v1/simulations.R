#' Simulation of data with underlying clusters
#'
#' Simulates mixture multivariate Normal data with clusters of observations
#' (rows) sharing similar profiles along (a subset of) variables (columns). The
#' conditional independence structure between the variables can be simulated or
#' provided in argument \code{adjacency}. The same covariance is used across all
#' clusters. Independent variables are simulated by default
#' (\code{nu_within=0}).
#'
#' @inheritParams SimulateGraphical
#' @param n vector of the number of observations per cluster in the simulated
#'   data. The number of observations in the simulated data is \code{sum(n)}.
#' @param pk vector of the number of variables in the simulated data.
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional independence structure between variables.
#' @param theta_xc optional binary vector encoding which variables (columns)
#'   contribute to the clustering structure between observations (rows). If
#'   \code{theta_xc=NULL}, variables contributing to the clustering are sampled
#'   with probability \code{nu_x/sum(pk)}.
#' @param nu_xc expected proportion of variables contributing to the clustering
#'   over the total number of variables. This argument is only used if
#'   \code{theta_xc} is not provided.
#' @param ev_xc vector of marginal expected proportion of explained for each
#'   variable contributing to the clustering. This parameter is only used for
#'   variables with a nonzero entry in \code{theta_xc}.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis applied on the
#'   predictors. This is the largest eigenvalue of the correlation (if
#'   \code{scale=TRUE}) or covariance (if \code{scale=FALSE}) matrix divided by
#'   the sum of eigenvalues. If \code{ev_xx=NULL} (the default), the constant u
#'   is chosen by maximising the contrast of the correlation matrix.
#'
#' @seealso \code{\link{MakePositiveDefinite}}, \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{sum(n)}
#'   observation and \code{sum(pk)} variables} \item{theta}{simulated (true)
#'   cluster membership.} \item{theta}{adjacency matrix of the graph encoding
#'   the conditional independence structure between variables.}
#'   \item{theta_xc}{binary vector encoding variables contributing to the
#'   clustering structure.} \item{ev}{vector of marginal expected proportions of
#'   explained variance for each variable.}
#'
#' @examples
#' \dontrun{
#' ## Example with 3 clusters
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15), ev_xc = 0.95
#' )
#' print(simul)
#' plot(simul)
#'
#'
#' ## Example with 2 variables contributing to clustering
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8))
#' )
#' print(simul)
#'
#' # Visualisation of the data
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = simul$data,
#'   colours = c("navy", "white", "red")
#' )
#' simul$ev # marginal proportions of explained variance
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with more distinct clusters
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.9, 0.8, rep(0, 8))
#' )
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with correlated contributors
#'
#' # Data simulation
#' pk <- 10
#' adjacency <- matrix(0, pk, pk)
#' adjacency[1, 2] <- adjacency[2, 1] <- 1
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = pk,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.9, 0.8, rep(0, 8)),
#'   adjacency = adjacency,
#'   pd_strategy = "min_eigenvalue",
#'   v_within = 0.6, v_sign = -1
#' )
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#' # Checking marginal proportions of explained variance
#' mymodel <- lm(simul$data[, 1] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#' mymodel <- lm(simul$data[, 2] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#' }
#' @export
SimulateClustering <- function(n = c(10, 10), pk = 10, adjacency = NULL,
                               theta_xc = NULL, nu_xc = 0.1, ev_xc = NULL,
                               implementation = HugeAdjacency, topology = "random",
                               nu_within = 0, nu_between = NULL,
                               v_within = c(0.5, 1), v_between = c(0, 0.1),
                               v_sign = c(-1, 1), continuous = TRUE,
                               pd_strategy = "diagonally_dominant", ev_xx = NULL, scale = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                               output_matrices = FALSE) {
  # Checking the inputs
  if (!is.null(theta_xc)){
    if (sum(pk)!=length(theta_xc)){
      warning("Arguments 'pk' and 'theta_xc' are not compatible. Argument 'pk' has been set to length('theta_xc').")
      pk=length(theta_xc)
    }
  }
  if (!is.null(adjacency)){
    if (ncol(adjacency)!=nrow(adjacency)){
      stop("Invalid input for argument 'adjacency'. It must be a square matrix (same number of rows and columns).")
    }
    if (sum(pk)!=ncol(adjacency)){
      warning("Arguments 'pk' and 'theta_xc' are not compatible. Argument 'pk' has been set to ncol('adjacency').")
      pk=ncol(adjacency)
    }
  }
  
  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = sum(n), pk = pk, theta = adjacency,
    implementation = implementation,
    topology = topology,
    nu_within = nu_within,
    nu_between = nu_between,
    output_matrices = output_matrices,
    v_within = v_within,
    v_between = v_between,
    continuous = continuous,
    pd_strategy = pd_strategy, ev_xx = ev_xx, scale = scale,
    u_list = u_list, tol = tol
  )
  
  # Defining number of clusters
  nc <- length(n)
  
  # Defining variables contributing to the clustering
  if (is.null(theta_xc)) {
    theta_xc <- SamplePredictors(pk = sum(pk), q = 1, nu = nu_xc, orthogonal = TRUE)[, 1]
  }
  
  # Simulating marginal proportions of explained variance
  if (is.null(ev_xc)) {
    ev_xc <- stats::runif(n = sum(pk))
  } else {
    if (length(ev_xc) == 1) {
      ev_xc <- rep(ev_xc, sum(pk))
    }
  }
  
  # Re-naming the outputs
  out$adjacency <- out$theta
  
  # Definition of membership
  theta <- NULL
  for (i in 1:length(n)) {
    theta <- c(theta, rep(i, each = n[i]))
  }
  names(theta) <- rownames(out$data)
  out$theta <- theta
  
  # # Building binary cluster membership for each feature
  # V <- stats::model.matrix(~ as.factor(theta) - 1)
  
  # Simulating the cluster-specific means
  if (length(n)>1){
    mu_mat <- matrix(NA, nrow = sum(n), ncol = sum(pk))
    for (k in 1:ncol(mu_mat)) {
      # Defining variance to reach expected proportion of e.v.
      var_mu <- ev_xc[k] * 1 / (1 - ev_xc[k])
      
      # Sampling initial values for cluster-specific means
      mu <- stats::rnorm(n = nc, mean = 0, sd = 1)
      for (i in 1:nrow(mu_mat)) {
        mu_mat[i, k] <- mu[theta[i]]
      }
      
      # Scaling to ensure mean of zero and defined variance
      mu_mat[, k] <- scale(mu_mat[, k])
      mu_mat[, k] <- mu_mat[, k] * sqrt(var_mu)
      # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-mean(mu_mat[,k]))/sqrt(1/(nrow(mu_mat)-1)*sum((mu_mat[, k]-mean(mu_mat[, k]))^2))
      # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-1/nrow(mu_mat)*sum(table(theta)*mu))/sqrt(1/(nrow(mu_mat)-1)*sum((mu_mat[, k]-mean(mu_mat[, k]))^2))
      # Equivalent to: sqrt(var_mu)*(mu_mat[, k]-1/nrow(mu_mat)*sum(table(theta)*mu))/sqrt(1/(nrow(mu_mat)-1)*sum(table(theta)*(mu-1/nrow(mu_mat)*sum(table(theta)*mu))^2))
    }
    
    # Using cluster-specific mean for contributing variables
    for (k in 1:ncol(mu_mat)) {
      if (theta_xc[k] == 1) {
        out$data[, k] <- scale(out$data[, k]) + mu_mat[, k]
      }
    }
  }
  
  # Scaling the output data
  out$data=scale(out$data)
  
  # Definition of contributing variables
  names(theta_xc) <- colnames(out$data)
  out$theta_xc <- theta_xc
  out$ev <- ev_xc * theta_xc
  
  # Defining the class
  class(out) <- "simulation_clustering"
  
  return(out)
}
