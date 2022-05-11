#' Calibrated consensus matrix
#'
#' Extracts the (calibrated) consensus matrix.
#'
#' @param stability output of \code{\link{Clustering}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated model is used.
#'
#' @return A binary and symmetric adjacency matrix encoding an undirected graph
#'   with no self-loops.
#'
#' @family calibration functions
#' @seealso \code{\link{Clustering}}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of data with clusters
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15),
#'   pk = 10, nu_xc = 1,
#'   ev_xc = 0.95
#' )
#'
#' # Consensus clustering
#' stab <- Clustering(xdata = x, implementation = HierarchicalClustering)
#' ConsensusMatrix(stab)
#' }
#'
#' @export
ConsensusMatrix <- function(stability, argmax_id = NULL) {
  if (class(stability) != "clustering") {
    stop("Invalid input for argument 'stability'. Only applicable to an object of class 'clustering', i.e. to the utput of Clustering().")
  }

  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability = stability)
  }
  mat <- stability$coprop[, , argmax_id[1]]

  return(mat)
}


#' Stable cluster membership
#'
#' Extracts (calibrated) stable clusters. These correspond to connected
#' components of the graph defined from stable co-membership.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{Clustering}}.
#'
#' @return A vector of cluster memberships.
#'
#' @family calibration functions
#' @seealso \code{\link{BiSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 50, n = 10)
#' }
#' @export
Clusters <- function(stability, argmax_id = NULL) {
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability = stability)
  }

  # Calibrated consensus matrix
  coprop <- ConsensusMatrix(stability = stability, argmax_id = argmax_id[1])

  # Extracting stable clusters from hierarchical clustering
  shclust <- stats::hclust(as.dist(1 - coprop), method = stability$methods$linkage)
  mymembership <- stats::cutree(shclust, k = stability$nc[argmax_id[1], 1])

  return(mymembership)
}


WeightBoxplot <- function(stability, at = NULL, argmax_id = NULL,
                          col = NULL, boxwex = 0.3,
                          xaxt = "s", xlab="", ylab = "Weight", cex.lab = 1.5,
                          las = 3, frame = "F", add = FALSE) {
  # Defining default colours
  if (is.null(col)) {
    col <- "navy"
  }

  # Extracting ID of calibrated parameters
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
  }

  # Extracting weights
  y <- stability$Beta[argmax_id[1], , ]

  # Removing zero weights (for sparse methods)
  y[which(y == 0)] <- NA

  if (is.null(at)) {
    x <- 1:nrow(y)
  } else {
    x <- at
    xaxt <- "n"
  }

  # Showing the distribution over nonzero weights
  boxplot(t(y),
    at = x, xlim = range(x),
    col = col, boxcol = col, whiskcol = col,
    staplecol = col, medcol = darken(col, amount = 0.5),
    whisklty = 1, range = 0, las = las, add = add,
    xlab=xlab, ylab = ylab, cex.lab = cex.lab, frame = frame,
    boxwex = boxwex, xaxt = xaxt
  )
  if (!is.null(at)) {
    axis(side = 1, at = axTicks(side = 1))
  }
}
