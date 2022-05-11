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
ArgmaxId <- function(stability = NULL, clustering = FALSE, S = NULL) {
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
      if (length(unique(stability$Lambda)) > 1) {
        # Identifying best number of contributing variables
        lambda_hat <- stability$Lambda[which.max(stability$S), 1]
        ids <- which(as.character(stability$Lambda) == lambda_hat)
      } else {
        ids <- 1:nrow(stability$Sc)
      }
      Sc <- stability$Sc[ids, 1]
      Sc_2d <- stability$Sc_2d[ids, , drop = FALSE]

      # Identifying best number of clusters
      argmax_id <- matrix(NA, nrow = 1, ncol = 2)
      id <- which.max(Sc)
      argmax_id[, 1] <- ids[id]
      tmpSc <- Sc_2d[id, ]
      argmax_id[, 2] <- which.max(tmpSc)
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
  colnames(argmax_id) <- c("lambda_id", "pi_id")
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
Argmax <- function(stability, clustering = FALSE) {
  if (class(stability) == "bi_selection") {
    argmax <- stability$summary
    argmax <- argmax[, colnames(argmax) != "S", drop = FALSE]
  } else {
    argmax <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
    if (clustering) {
      id <- ArgmaxId(stability = stability, clustering = clustering)
      argmax[, 1] <- stability$nc[id[1], 1]
      argmax[, 2] <- stability$params$pi_list[id[2]]
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
    colnames(argmax) <- c("lambda", "pi")
  }

  return(argmax)
}


#' Calibrated adjacency matrix
#'
#' Extracts the adjacency matrix of the (calibrated) stability selection
#' graphical model.
#'
#' @param stability output of \code{\link{GraphicalModel}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated model is used.
#'
#' @return A binary and symmetric adjacency matrix encoding an undirected graph
#'   with no self-loops.
#'
#' @family calibration functions
#' @seealso \code{\link{GraphicalModel}}
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
#' # Calibrated adjacency matrix
#' A <- Adjacency(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(20, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' A <- Adjacency(stab, argmax_id = myids)
#' }
#'
#' @export
Adjacency <- function(stability, argmax_id = NULL) {
  if (class(stability) == "bi_selection") {
    if ("selectedY" %in% names(stability)) {
      A <- Square(t(rbind(stability$selectedX, stability$selectedY)))
    } else {
      A <- Square(stability$selectedX)
    }
  } else {
    if (class(stability) == "graphical_model") {
      A <- matrix(0, ncol = ncol(stability$selprop), nrow = nrow(stability$selprop))
    } else {
      A <- matrix(0, ncol = ncol(stability$coprop), nrow = nrow(stability$coprop))
    }
    bigblocks <- BlockMatrix(stability$params$pk)
    if (is.null(argmax_id)) {
      if (class(stability) == "graphical_model") {
        argmax_id <- ArgmaxId(stability)
        argmax <- Argmax(stability)
      } else {
        argmax_id <- ArgmaxId(stability, clustering = TRUE)
        argmax <- Argmax(stability, clustering = TRUE)
      }
    } else {
      argmax <- NULL
      for (block_id in 1:ncol(stability$Lambda)) {
        argmax <- rbind(argmax, c(
          stability$Lambda[argmax_id[block_id, 1], block_id],
          stability$params$pi_list[argmax_id[block_id, 2]]
        ))
      }
    }
    for (block_id in 1:ncol(stability$Lambda)) {
      if (class(stability) == "graphical_model") {
        A_block <- ifelse(stability$selprop[, , argmax_id[block_id, 1]] >= argmax[block_id, 2], 1, 0)
      } else {
        A_block <- ifelse(stability$coprop[, , argmax_id[block_id, 1]] >= argmax[block_id, 2], 1, 0)
      }
      # A_block[lower.tri(A_block)] <- 0
      # A_block <- A_block + t(A_block) # for symmetry
      if (length(stability$params$pk) > 1) {
        A_block[bigblocks != block_id] <- 0
      }
      A <- A + A_block
    }
  }
  A[is.na(A)] <- 0
  return(A)
}


#' Calibrated adjacency matrix
#'
#' Extracts the adjacency matrix of the (calibrated) stability selection
#' graphical model.
#'
#' @param stability output of \code{\link{GraphicalModel}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated model is used.
#'
#' @return A binary and symmetric adjacency matrix encoding an undirected graph
#'   with no self-loops.
#'
#' @family calibration functions
#' @seealso \code{\link{GraphicalModel}}
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
    argmax_id <- ArgmaxId(stability = stability, clustering = TRUE)
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
#' @param adjacency adjacency matrix or output of \code{\link{Clustering}}.
#' @param implementation function to use for community detection. This function
#'   must use \code{graph} as argument. Recommended functions include the
#'   following from \code{\link[igraph:igraph-package]{igraph}}:
#'   \code{\link[igraph]{components}} (connected components),
#'   \code{\link[igraph]{cluster_louvain}} (Louvain method),
#'   \code{\link[igraph]{cluster_walktrap}} (using random walks).
#' @param weighted logical indicating if edges should be weighted by
#'   co-membership proportions for the community detection (if applicable). Only
#'   used if \code{adjacency} is the output of \code{Clustering}.
#' @param ... additional parameters passed to the functions provided in
#'   \code{implementation}.
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
#' # Consensus clustering
#' stab <- GraphicalModel(
#'   xdata = simul$data,
#'   Lambda = seq(2, ncol(simul$data)),
#'   implementation = HierarchicalClustering
#' )
#'
#' # Stable cluster membership
#' groups <- Clusters(stab)
#'
#' # Network representation of stable co-membership
#' set.seed(1)
#' plot(Graph(CoMembership(groups),
#'   satellites = TRUE,
#'   node_colour = groups
#' ))
#' }
#' @export
Clusters <- function(adjacency = NULL, argmax_id = NULL,
                     weighted = FALSE,
                     implementation = igraph::components, ...) {
  if (!is.matrix(adjacency)) {
    # Computing stable co-membership matrix
    consensus <- ConsensusMatrix(adjacency, argmax_id = argmax_id)
    adjacency <- Adjacency(stability = adjacency, argmax_id = argmax_id)
  } else {
    if (weighted) {
      warning("The weights provided in 'adjacency' are used.")
      consensus <- adjacency
    }
  }

  # Extracting stable connected components
  # mymembership <- igraph::components(Graph(adjacency, satellites = TRUE))$membership
  if (weighted & ("weights" %in% names(formals(implementation)))) {
    mycommunities <- do.call(implementation, args = list(
      graph = Graph(adjacency, satellites = TRUE),
      weights = consensus,
      ...
    ))
  } else {
    mycommunities <- do.call(implementation, args = list(
      graph = Graph(adjacency, satellites = TRUE),
      ...
    ))
  }
  mymembership <- igraph::membership(mycommunities)

  return(mymembership)
}
