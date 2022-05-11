#' (Sparse) clustering algorithm
#'
#' Runs the (sparse) clustering algorithm specified in the argument
#' \code{implementation} and returns matrices of selected variables, variable
#' weights, and the co-membership structure. This function is not using
#' stability.
#'
#' @inheritParams Clustering
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc}
#'   is not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param Lambda vector of penalty parameters.
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{weight}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'   \item{comembership}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of 15 observations belonging to 3 groups
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100
#' )
#'
#' # Running hierarchical clustering
#' myclust <- ClusteringAlgo(
#'   xdata = simul$data, nc = 2:5,
#'   implementation = HierarchicalClustering
#' )
#'
#' # Running sparse hierarchical clustering
#' myclust <- ClusteringAlgo(
#'   xdata = simul$data,
#'   Lambda = c(1.5, 2), nc = 2:5,
#'   implementation = SparseHierarchicalClustering
#' )
#' }
#'
#' @export
ClusteringAlgo <- function(xdata,
                           Lambda = NULL, nc,
                           implementation = HierarchicalClustering, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in 1:ncol(xdata)) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Applying user-defined function for variable selection
  out <- do.call(implementation, args = list(xdata = xdata, nc = nc, Lambda = Lambda, ...))

  if ("weight" %in% names(out)) {
    beta_full <- out$weight

    # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
    if (any(mysd == 0)) {
      selected[, which(mysd == 0)] <- 0
      if (length(dim(beta_full)) == 2) {
        beta_full[, which(mysd == 0)] <- 0
      }
      if (length(dim(beta_full)) == 3) {
        beta_full[, which(mysd == 0), ] <- 0
      }
    }
  } else {
    beta_full <- matrix(1, nrow = length(nc), ncol = ncol(xdata))
    rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))
    colnames(beta_full) <- colnames(xdata)
  }

  return(list(comembership = out$comembership, weight = beta_full))
}
