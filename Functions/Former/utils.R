#' Pairwise co-membership
#'
#' Generates a symmetric and binary matrix indicating if a pair of features
#' belongs to the same cluster.
#'
#' @param groups vector of group membership.
#'
#' @return A symmetric and binary matrix.
#'
#' @examples
#' \dontrun{
#' # Simulated grouping structure
#' mygroups <- c(rep(1, 3), rep(2, 5), rep(3, 2))
#'
#' # Co-membership matrix
#' CoMembership(mygroups)
#' }
#'
#' @export
CoMembership <- function(groups) {
  if (length(unique(groups)) > 1) {
    # Building binary cluster membership for each feature
    V <- stats::model.matrix(~ as.factor(groups) - 1)

    # Building cluster co-membership
    comembership <- V %*% t(V)
  } else {
    comembership <- matrix(1, nrow = length(groups), ncol = length(groups))
  }

  # Re-formatting co-membership matrix
  diag(comembership) <- 0
  rownames(comembership) <- colnames(comembership) <- names(groups)

  return(comembership)
}
