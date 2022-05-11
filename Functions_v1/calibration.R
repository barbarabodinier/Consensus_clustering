#' Calibrated parameter indices
#'
#' Extracts the indices of calibrated parameters with respect to the grids
#' provided in \code{Lambda} and \code{pi_list} in \code{stability}.
#'
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}. If \code{stability=NULL}, \code{S} must be
#'   provided.
#' @param clustering logical indicating whether indices of calibrated clustering
#'   parameters (\code{clustering=TRUE}) or variable selection
#'   (\code{clustering=FALSE}) should be extracted. This argument is only used
#'   if \code{stability} is the output from \code{\link{Clustering}}.
#' @param S matrix of stability scores obtained with different combinations of
#'   parameters where rows correspond to different values of the parameter
#'   controlling the level of sparsity in the underlying feature selection
#'   algorithm and columns correspond to different values of the threshold in
#'   selection proportions. If \code{S=NULL}, argument \code{stability} must be
#'   provided.
#'
#' @return A matrix of parameter indices. For multi-block graphical models, rows
#'   correspond to different blocks.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Extracting IDs of calibrated parameters
#' ids <- ArgmaxId(stab)
#' stab$Lambda[ids[1], 1]
#' stab$params$pi_list[ids[2]]
#'
#' # Alternative formulation
#' ids2 <- ArgmaxId(S = stab$S_2d)
#'
#' # Link with Argmax() function
#' args <- Argmax(stab)
#' }
#'
#' @export
ArgmaxId <- function(stability = NULL, clustering=FALSE, S = NULL) {
  if ((is.null(stability)) & (is.null(S))) {
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }
  if (clustering) {
    if (!is.null(S)) {
      stop("Invalid input. Argument 'stability' needs to be supplied with clustering = TRUE.")
    }
  }
  if (is.null(S)) {
    if (clustering) {
      # if (length(unique(stability$Lambda)) > 1) {
      #   # Identifying best number of contributing variables
      #   lambda_hat <- stability$Lambda[which.max(stability$S), 1]
      #   ids <- which(as.character(stability$Lambda) == lambda_hat)
      # } else {
      #   ids <- 1:nrow(stability$Sc)
      # }
      Sc <- stability$Sc
      # Sc_2d <- stability$Sc_2d[ids, , drop = FALSE]
      
      # Identifying best number of clusters
      # argmax_id <- matrix(NA, nrow = 1, ncol = 2)
      argmax_id <- as.matrix(which.max(Sc), 1, 1)
      # argmax_id[, 1] <- ids[id]
      # tmpSc <- Sc_2d[id, ]
      # argmax_id[, 2] <- which.max(tmpSc)
    } else {
      argmax_id <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
      if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
        id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
        argmax_id[, 1] <- rep(id, nrow(argmax_id))
        for (block_id in 1:ncol(stability$Lambda)) {
          if (!is.na(stability$P[id, block_id])) {
            argmax_id[block_id, 2] <- which(stability$params$pi_list == stability$P[id, block_id])
          }
        }
      } else {
        for (block_id in 1:ncol(stability$Lambda)) {
          if (ncol(stability$Lambda) == 1) {
            myS <- stability$S
          } else {
            myS <- stability$S[, block_id, drop = FALSE]
          }
          myS[is.na(myS)] <- 0
          myid <- which.max(myS[, 1])
          argmax_id[block_id, ] <- c(myid, which(stability$params$pi_list == stability$P[myid, block_id]))
        }
      }
    }
  } else {
    argmax_id <- matrix(NA, nrow = 1, ncol = 2)
    myS <- apply(S, 1, max, na.rm = TRUE)
    myS[is.na(myS)] <- 0
    myid <- which.max(myS)
    argmax_id[1, ] <- c(myid, max(which(S[myid, ] == myS[myid])))
  }
  if (clustering){
    colnames(argmax_id)="nc_id"
  } else {
    colnames(argmax_id) <- c("lambda_id", "pi_id")
  }
  return(argmax_id)
}


#' Calibrated parameters
#'
#' Extracts calibrated parameter values in stability selection.
#'
#' @inheritParams ArgmaxId
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}}, or \code{\link{BiSelection}}.
#'
#' @return A matrix of parameter values. If applied to the output of
#'   \code{\link{VariableSelection}} or \code{\link{GraphicalModel}}, the first
#'   column (\code{lambda}) denotes the calibrated hyper-parameter of the
#'   underlying algorithm. The second column (\code{pi}) is the calibrated
#'   threshold in selection/co-membership proportions. For multi-block graphical
#'   models, rows correspond to different blocks. If applied to the output of
#'   \code{\link{BiSelection}}, all columns are named as in object
#'   \code{summary}.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}
#'
#' @examples
#' \dontrun{
#'
#' ## Graphical modelling
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Extracting calibrated parameters
#' Argmax(stab)
#' }
#'
#' @export
Argmax <- function(stability, clustering=FALSE) {
  if (class(stability) == "bi_selection") {
    argmax <- stability$summary
    argmax <- argmax[, colnames(argmax) != "S", drop = FALSE]
  } else {
    argmax <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 1)
    if (clustering) {
      id <- ArgmaxId(stability = stability, clustering=clustering)
      argmax[, 1] <- stability$nc[id[1], 1]
      # argmax[, 2] <- stability$params$pi_list[id[2]]
    } else {
      if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
        id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
        argmax[, 1] <- stability$Lambda[id, ]
        argmax[, 2] <- stability$P[id, ]
      } else {
        for (block_id in 1:ncol(stability$Lambda)) {
          if (ncol(stability$Lambda) == 1) {
            myS <- stability$S
          } else {
            myS <- stability$S[, block_id, drop = FALSE]
          }
          myS[is.na(myS)] <- 0
          myid <- which.max(myS[, 1])
          argmax[block_id, ] <- c(stability$Lambda[myid, block_id], stability$P[myid, block_id])
        }
      }
    }
    if (clustering){
      colnames(argmax)="nc"
    } else {
      colnames(argmax) <- c("lambda", "pi")
    }
  }
  
  return(argmax)
}


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
#' simul=SimulateClustering(n=c(10,30,15),
#'                          pk=10, nu_xc=1,
#'                          ev_xc=0.95)
#'
#' # Consensus clustering
#' stab=Clustering(xdata=x, implementation=HierarchicalClustering)
#' ConsensusMatrix(stab)
#' }
#'
#' @export
ConsensusMatrix <- function(stability, argmax_id = NULL) {
  if (class(stability)!="clustering"){
    stop("Invalid input for argument 'stability'. Only applicable to an object of class 'clustering', i.e. to the utput of Clustering().")
  }
  
  if (is.null(argmax_id)){
    argmax_id=ArgmaxId(stability=stability, clustering=TRUE)
  }
  mat=stability$coprop[,,argmax_id[1]]
  
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
#'
#' }
#' @export
Clusters <- function(stability, argmax_id = NULL) {
  if (is.null(argmax_id)){
    argmax_id=ArgmaxId(stability=stability, clustering=TRUE)
  }
  
  # Calibrated consensus matrix
  coprop=ConsensusMatrix(stability = stability, argmax_id = argmax_id[1])
  
  # Extracting stable clusters from hierarchical clustering
  shclust=stats::hclust(as.dist(1-coprop), method = stability$methods$linkage)
  mymembership=stats::cutree(shclust, k=stability$nc[argmax_id[1], 1])
  
  return(mymembership)
}
