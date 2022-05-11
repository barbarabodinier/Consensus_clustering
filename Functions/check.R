#' Checking input parameters (graphical model)
#'
#' Checks if input parameters are valid. For invalid parameters, this function
#' (i) stops the run and generates an error message, or (ii) sets the invalid
#' parameter to its default value and reports it in a warning message.
#'
#' @inheritParams Clustering
#'
#' @keywords internal
CheckInputClustering <- function(xdata, Lambda = NULL, 
                                pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3,
                                implementation = HierarchicalClustering, scale = TRUE,
                                resampling = "subsampling", 
                                verbose = TRUE) {
  # List of arguments
  myargs <- c(
    "xdata", "Lambda", 
    "pi_list", "K", "tau", "seed", "n_cat", 
    "scale",
    "verbose"
  )
  
  # Checking the inputs (xdata)
  xdata <- as.matrix(xdata)
  if (sum(is.na(xdata)) > 0) {
    stop("Invalid input for argument 'xdata'. Missing values are not allowed in 'xdata'.")
  }
  if ((nrow(xdata) < 10) | (ncol(xdata) < 1)) {
    stop("Invalid input for argument 'xdata'. Not enough xdata.")
  }
  
  # Checking the inputs (pi_list)
  pi_list <- sort(pi_list)
  if (n_cat == 3) {
    if (any(pi_list > 0.5) & any(pi_list < 1)) {
      if ((min(pi_list) < 0.5) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0.5) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list > 0) & any(pi_list < 1)) {
      if ((min(pi_list) < 0) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }
  
  # Checking the inputs (K)
  K <- as.numeric(K)
  if ((length(K) != 1) | is.na(K)) {
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K <- 100
  }
  
  # Checking the inputs (tau)
  tau <- as.numeric(tau)
  if ((length(tau) != 1) | is.na(tau) | (tau >= 1) | (tau <= 0)) {
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau <- 0.5
  }
  
  # Checking the inputs (seed)
  seed <- as.numeric(seed)
  if ((length(seed) != 1) | is.na(seed)) {
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed <- 1
  }

  # Checking the inputs (implementation)
  if (!is.function(implementation)) {
    stop("Invalid input for argument 'implementation'. This argument must be a function to use for graphical modelling.")
  }
  
  # Checking the inputs (scale)
  scale <- as.logical(scale)
  if ((length(scale) != 1) | is.na(scale)) {
    stop("Invalid input for argument 'scale'. The argument 'scale' must be logical (TRUE or FALSE).")
  }
  
  # Checking the inputs (resampling)
  if ((!is.function(resampling)) & (!is.character(resampling))) {
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }
  
  # Checking the inputs (verbose)
  verbose <- as.logical(verbose)
  if ((length(verbose) != 1) | is.na(verbose)) {
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose <- TRUE
  }
  
  # Checking the inputs (Lambda)
  nblocks=1
  if (!is.null(Lambda)) {
    if (is.matrix(Lambda)) {
      if ((ncol(Lambda) != nblocks) & (ncol(Lambda) != 1)) {
        stop(paste0("Invalid input for argument 'Lambda'. The argument 'Lambda' must be a matrix as many columns as blocks (N=", nblocks, ")."))
      }
      if (ncol(Lambda) == 1) {
        Lambda <- as.numeric(as.vector(Lambda))
      } else {
        Lambda_copy <- Lambda
        Lambda <- NULL
        for (k in 1:ncol(Lambda_copy)) {
          Lambda <- cbind(Lambda, as.numeric(Lambda_copy[, k]))
        }
      }
    } else {
      Lambda <- as.numeric(Lambda)
    }
    if (any(is.na(Lambda))) {
      if (all(is.na(Lambda))) {
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda <- as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
  }
  
  # Assigning checked values to the parent function
  for (i in 1:length(myargs)) {
    assign(myargs[i], get(myargs[i]), envir = parent.frame(n = 1))
  }
}

