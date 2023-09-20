InternalValidationIndex <- function(xdata, grouping, index = "silhouette") {
  if (index == "silhouette") {
    mysilhouette <- cluster::silhouette(
      x = grouping,
      dist = dist(x = xdata, method = "euclidian")
    )
    out <- mean(mysilhouette[, 3])
  }
  if (index == "ch") {
    out <- fpc::calinhara(
      x = xdata,
      clustering = grouping,
      cn = length(unique(grouping))
    )
  }
  if (index == "db") {
    out <- clusterSim::index.DB(x = xdata, cl = grouping, p = 2, q = 1)$DB
    # Note: DB to be minimised
  }
  return(out)
}


InternalCalibration <- function(xdata, stability = NULL,
                                nc_max = 20, method = "hclust", linkage = "complete",
                                index = "silhouette") {
  if (is.null(stability)) {
    dist <- dist(xdata, method = "euclidian")
    score <- rep(NA, nc_max)

    if (method == "hclust") {
      hclust <- hclust(d = dist, method = linkage)
      for (k in 2:nc_max) {
        myclusters <- cutree(hclust, k = k)
        score[k] <- InternalValidationIndex(
          xdata = xdata,
          grouping = myclusters,
          index = index
        )
      }
    }

    if (method == "pam") {
      for (k in 2:nc_max) {
        set.seed(1)
        myclusters <- cluster::pam(x = dist, k = k, diss = TRUE, cluster.only = TRUE)
        score[k] <- InternalValidationIndex(
          xdata = xdata,
          grouping = myclusters,
          index = index
        )
      }
    }

    if (method == "kmeans") {
      for (k in 2:nc_max) {
        set.seed(1)
        mykmeans <- stats::kmeans(x = xdata, centers = k)
        myclusters <- mykmeans$cluster
        score[k] <- InternalValidationIndex(
          xdata = xdata,
          grouping = myclusters,
          index = index
        )
      }
    }

    if (method == "gmm") {
      for (k in 2:nc_max) {
        set.seed(1)
        myclust <- mclust::Mclust(data = xdata, G = k, verbose = FALSE)
        myclusters <- myclust$classification
        score[k] <- InternalValidationIndex(
          xdata = xdata,
          grouping = myclusters,
          index = index
        )
      }
    }
  } else {
    nc_max <- max(stability$nc)
    score <- rep(NA, nc_max)
    nc_list <- stability$nc
    nc_list <- nc_list[which(nc_list != 1)]
    for (k in nc_list) {
      myclusters <- Clusters(stability = stability, argmax_id = k)
      score[k] <- InternalValidationIndex(
        xdata = xdata,
        grouping = myclusters,
        index = index
      )
    }
  }

  return(score)
}


InformationTheoryCalibration <- function(xdata, nc_max = 20, index = "bic") {
  score <- rep(NA, nc_max)

  for (k in 2:nc_max) {
    set.seed(1)
    myclust <- mclust::Mclust(data = xdata, G = k, verbose = FALSE)
    score[k] <- myclust[[index]]
  }

  return(score)
}


GMM0Calibration <- function(xdata, nc_max = 20, index = "bic") {
  score <- rep(NA, nc_max)

  for (k in 2:nc_max) {
    set.seed(1)
    myclust <- mclust::Mclust(data = xdata, G = k, verbose = FALSE)
    score[k] <- myclust[[index]]
    if (which.max(score) == k) {
      myclusters <- myclust$classification
    }
  }

  return(list(scores = score, clusters = myclusters))
}


GMM1Calibration <- function(xdata, nc_max = 20, search = "headlong", index = "bic") {
  score <- rep(NA, nc_max)
  for (k in 2:nc_max) {
    myclustvarsel <- clustvarsel(xdata,
      search = search,
      G = k,
      verbose = FALSE
    )
    score[k] <- myclustvarsel$model[[index]]
    if (which.max(score) == k) {
      myclusters <- myclustvarsel$model$classification
      selected <- rep(0, ncol(xdata))
      names(selected) <- colnames(xdata)
      selected[names(myclustvarsel$subset)] <- 1
    }
  }
  return(list(scores = score, clusters = myclusters, selected = selected))
}


GMM2Calibration <- function(xdata, nc_max = 20, index = "micl") {
  score <- rep(NA, nc_max)
  for (k in 2:nc_max) {
    print(k)
    myvarsellcm <- VarSelCluster(
      x = xdata,
      gvals = k,
      vbleSelec = TRUE,
      crit.varsel = toupper(index)
    )
    if (index == "micl") {
      score[k] <- MICL(myvarsellcm)
    }
    if (index == "bic") {
      score[k] <- BIC(myvarsellcm)
    }
    if (index == "aic") {
      score[k] <- AIC(myvarsellcm)
    }
    if (which.max(score) == k) {
      myclusters <- fitted(myvarsellcm)
      selected <- rep(0, ncol(xdata))
      names(selected) <- colnames(xdata)
      selected[myvarsellcm@model@names.relevant] <- 1
    }
  }
  return(list(scores = score, clusters = myclusters, selected = selected))
}


SilhouetteScore <- function(x, nc_max = 20, method = "hclust", linkage = "complete") {
  dist <- dist(x, method = "euclidian")
  score <- rep(NA, nc_max)

  if (method == "hclust") {
    hclust <- hclust(d = dist, method = linkage)
    for (k in 2:nc_max) {
      myclusters <- cutree(hclust, k = k)
      mysilhouette <- silhouette(x = myclusters, dist = dist)
      score[k] <- mean(mysilhouette[, 3])
    }
  }

  if (method == "kmeans") {
    for (k in 2:nc_max) {
      set.seed(1)
      mykmeans <- stats::kmeans(x = x, centers = k)
      myclusters <- mykmeans$cluster
      mysilhouette <- silhouette(x = myclusters, dist = dist)
      score[k] <- mean(mysilhouette[, 3])
    }
  }

  if (method == "pam") {
    for (k in 2:nc_max) {
      set.seed(1)
      myclusters <- cluster::pam(x = dist, k = k, diss = TRUE, cluster.only = TRUE)
      mysilhouette <- silhouette(x = myclusters, dist = dist)
      score[k] <- mean(mysilhouette[, 3])
    }
  }

  return(score)
}


hclusCut <- function(x, k, d.meth = "euclidean", linkage = "complete", ...) {
  return(list(cluster = cutree(hclust(dist(x, method = d.meth), method = linkage, ...), k = k)))
}


kmeansCut <- function(x, k, d.meth = "euclidean", ...) {
  set.seed(1)
  return(list(cluster = stats::kmeans(x, centers = k)$cluster))
}


gmmCut <- function(x, k, d.meth = "euclidean", ...) {
  set.seed(1)
  return(list(cluster = mclust::Mclust(data = simul$data, G = id, verbose = FALSE)$classification))
}


pamCut <- function(x, k, d.meth = "euclidean", ...) {
  set.seed(1)
  return(list(cluster = cluster::pam(x = dist(x, method = d.meth), k = k, diss = TRUE, cluster.only = TRUE)))
}


GapStatistic <- function(xdata, nc_max = 20, iters = 25, method = "hclust", linkage = "complete") {
  if (method == "hclust") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = hclusCut,
      linkage = linkage,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  if (method == "kmeans") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = kmeansCut,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  if (method == "gmm") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = gmmCut,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  if (method == "pam") {
    out <- clusGap(
      x = as.matrix(xdata),
      FUNcluster = pamCut,
      K.max = nc_max, B = iters,
      verbose = FALSE
    )
  }
  return(as.data.frame(out$Tab))
}


AreaUnderCDF <- function(M, thr_list = seq(0.01, 1, by = 0.01)) {
  n <- nrow(M)
  x <- M[upper.tri(M)]
  CDF <- rep(NA, length(thr_list))
  for (k in 1:length(thr_list)) {
    CDF[k] <- sum(x <= thr_list[k]) / (n * (n - 1) / 2)
  }

  area <- 0
  for (k in 2:length(CDF)) {
    area <- area + (thr_list[k] - thr_list[k - 1]) * CDF[k]
  }

  return(list(CDF = CDF, area = area))
}


DeltaArea <- function(areas) {
  delta <- areas[1]
  for (k in 2:(length(areas))) {
    delta <- c(delta, (areas[k] - areas[k - 1]) / areas[k - 1])
  }
  names(delta) <- names(areas)
  return(delta)
}


DeltaAreaCDF <- function(stability, thr_list = seq(0.01, 1, by = 0.01)) {
  areas <- rep(NA, dim(stability$coprop)[3])
  for (k in 2:dim(stability$coprop)[3]) {
    areas[k] <- AreaUnderCDF(M = stability$coprop[, , k])$area
  }
  delta <- DeltaArea(areas[-1])
  delta <- c(NA, delta)
  names(delta) <- stability$nc[, 1]
  return(delta)
}


PAC <- function(stab, x1 = 0.1, x2 = 0.9) {
  maxK <- max(stab$nc)
  Kvec <- 2:maxK
  score <- rep(NA, length(Kvec))
  for (i in Kvec) {
    M <- stab$coprop[, , i]
    Fn <- ecdf(M[lower.tri(M)])
    score[i - 1] <- Fn(x2) - Fn(x1)
  }
  score <- c(NA, score)
  names(score) <- stab$nc
  return(score)
}


PINSDiscrepancy <- function(x, stab, method = "hclust", linkage = "complete") {
  # Clustering on original (full) data
  if (method == "hclust") {
    myclustering <- HierarchicalClustering(xdata = x, nc = stab$nc, linkage = linkage)
  }
  if (method == "pam") {
    myclustering <- PAMClustering(xdata = x, nc = stab$nc)
  }
  if (method == "kmeans") {
    myclustering <- KMeansClustering(xdata = x, nc = stab$nc)
  }
  if (method == "gmm") {
    myclustering <- GMMClustering(xdata = x, nc = stab$nc)
  }
  origS <- list()
  for (k in stab$nc) {
    origS <- c(origS, list(myclustering$comembership[, , k]))
  }
  origS[[1]] <- list(NULL)

  # Clustering on perturbed data (subsamples)
  pertS <- list()
  for (k in 1:dim(stab$coprop)[3]) {
    pertS <- c(pertS, list(stab$coprop[, , k]))
  }
  pertS[[1]] <- list(NULL)

  # Calculating PINS discrepancy
  discrepancy <- PINSPlus:::CalcPerturbedDiscrepancy(
    origS = origS,
    pertS = pertS,
    clusRange = 2:max(stab$nc)
  )$AUC

  return(discrepancy)
}


MonteCarloScore <- function(x, stab, iters = 25, objective = "entropy", method = "hclust", linkage = "complete") {
  if (method == "hclust") {
    clusteralg <- "hc"
  }
  if (method == "kmeans") {
    clusteralg <- "km"
  }
  if (method == "pam") {
    clusteralg <- "pam"
  }

  # Running M3C for reference distribution
  out <- M3C(
    mydata = t(x),
    iters = iters,
    clusteralg = clusteralg,
    innerLinkage = linkage,
    finalLinkage = "complete",
    maxK = max(stab$nc),
    pItem = stab$params$tau,
    repsref = stab$params$K,
    repsreal = 2,
    seed = 1,
    objective = objective,
    removeplots = TRUE,
    silent = TRUE
  )

  # Calculation of PAC scores
  if (objective == "PAC") {
    real <- data.frame(K = stab$nc, PAC_REAL = PAC(stab))
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Calculation of entropy scores
  if (objective == "entropy") {
    entropies <- rep(NA, dim(stab$coprop)[3])
    for (i in 2:dim(stab$coprop)[3]) {
      entropies[i] <- M3C:::entropy(stab$coprop[, , i])
    }
    real <- data.frame(K = stab$nc, PAC_REAL = entropies)
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Storing the reference
  ls <- out$refpacscores

  # Calculating reference mean
  real$PAC_REF <- colMeans(ls)

  # Checking PAC values
  ptemp <- real$PAC_REAL
  ptemp[ptemp == 0] <- 0.0001
  pacreal <- ptemp

  # Calculating RCSI score
  diffM <- sweep(log(ls), 2, log(pacreal))
  real$RCSI <- colMeans(diffM)
  real$RCSI_SE <- (apply(diffM, 2, sd)) / sqrt(nrow(ls))

  # Calculating p-value
  pvals <- vapply(seq_len(ncol(ls)), function(i) {
    distribution <- as.numeric(ls[, i])
    ((length(distribution[distribution < real$PAC_REAL[i]])) + 1) / (iters + 1) # (b+1)/(m+1)=pval
  }, numeric(1))
  real$MONTECARLO_P <- pvals

  if (objective == "PAC") {
    variance <- apply(ls, 2, var)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      params2 <- M3C:::estBetaParams(mu = mean, var = var)
      pbeta(realpac, params2[[1]], params2[[2]])
    }, numeric(1))
    real$BETA_P <- pvals2
    real$P_SCORE <- -log10(real$BETA_P)
  } else if (objective == "entropy") {
    variance <- apply(ls, 2, sd)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      pnorm(realpac, mean = mean, sd = var)
    }, numeric(1))
    real$NORM_P <- pvals2
    real$P_SCORE <- -log10(real$NORM_P)
    colnames(real)[2:3] <- c("ENTROPY_REAL", "ENTROPY_REF")
  }

  real <- rbind(c(1, rep(NA, ncol(real) - 1)), real)
  return(real)
}


AllPerf <- function(stab) {
  perf <- NULL
  for (i in 1:dim(stab$coprop)[3]) {
    perf <- rbind(
      perf,
      ClusteringPerformance(
        theta = Clusters(stab, argmax_id = i),
        theta_star = simul
      )
    )
  }
  perf <- cbind(nc = stab$nc, perf)
  return(perf)
}


ScatterPerf <- function(x, perf, xaxt = "s", xlab = "", ylab = "ARI", col = "navy") {
  id <- ManualArgmaxId(x)
  mycolours <- rep(col, nrow(perf))
  mycolours[id] <- "darkred"
  plot(x, perf$ari,
    panel.first = c(
      abline(h = perf$ari[id], col = "darkred", lty = 3),
      abline(v = x[id], col = "darkred", lty = 3)
    ),
    # panel.first=c(abline(h=axisTicks(range(perf$ari), log=FALSE), lty=3, col="grey"),
    #               abline(v=axisTicks(range(delta, na.rm=TRUE), log=FALSE), lty=3, col="grey")),
    pch = 19, cex = 3,
    las = 1,
    col = colorspace::lighten(mycolours, amount = 0.8),
    xaxt = xaxt,
    cex.lab = 1.5,
    xlab = xlab,
    ylab = ylab
  )
  text(x, perf$ari, labels = perf$nc, col = colorspace::darken(mycolours, amount = 0.2))
}


ManualCalibPlot <- function(y, x = NULL, xaxt = "s", xlab = "Number of clusters", ylab = "", col = "navy") {
  if (is.null(x)) {
    x <- 1:length(y)
  }
  id <- ManualArgmaxId(y)
  mycolours <- rep(col, length(y))
  mycolours[id] <- "darkred"
  plot(x, y,
    panel.first = c(
      abline(h = y[id], col = "darkred", lty = 3),
      abline(v = x[id], col = "darkred", lty = 3)
    ),
    # panel.first=c(abline(h=axisTicks(range(perf$ari), log=FALSE), lty=3, col="grey"),
    #               abline(v=axisTicks(range(delta, na.rm=TRUE), log=FALSE), lty=3, col="grey")),
    pch = 19, cex = 3,
    las = 1,
    col = colorspace::lighten(mycolours, amount = 0.8),
    xaxt = xaxt,
    cex.lab = 1.5,
    xlab = xlab,
    ylab = ylab
  )
  text(x, y, labels = x, col = colorspace::darken(mycolours, amount = 0.2))
}


ManualArgmaxId <- function(x, digits = 10) {
  return(max(which(round(x, digits = digits) == max(round(x, digits = digits), na.rm = TRUE))))
}


# CalibrationCurve <- function(stability,
#                              bty = "o", xlab = "", ylab = "",
#                              col = c("navy", "forestgreen", "tomato"),
#                              legend = TRUE, ncol = 1) {
#   y <- stability$Sc
#   x <- stability$nc
#   z <- round(stability$Lambda, digits = 5)
#
#   mycolours <- colorRampPalette(col)(length(unique(z)))
#   names(mycolours) <- unique(z)
#   plot(NA,
#     xlim = c(0, max(stability$nc)), ylim = c(0, 1),
#     xlab = xlab, ylab = ylab,
#     las = 1, cex.lab = 1.5, bty = bty
#   )
#   for (lambda in unique(z)) {
#     ids <- which(z == lambda)
#     points(x[ids], y[ids], pch = 18, col = mycolours[as.character(lambda)])
#     lines(x[ids], y[ids], lty = 1, lwd = 0.5, col = mycolours[as.character(lambda)])
#   }
#   abline(v = Argmax(stability)[1], lty = 2, col = "darkred")
#   if (legend) {
#     if (length(unique(stability$Q)) == 1) {
#       legend("topright",
#         legend = unique(formatC(stability$Lambda, format = "f", digits = 2)),
#         pch = 15, col = mycolours, bty = "n", title = expression(lambda), ncol = ncol
#       )
#     } else {
#       legend("topright",
#         legend = paste0(
#           unique(formatC(stability$Lambda[, 1], format = "f", digits = 2)),
#           " (", unique(stability$Q[, 1]), ")"
#         ),
#         pch = 15, col = mycolours, bty = "n", title = expression(lambda), ncol = ncol
#       )
#     }
#   }
# }


SparseHierarchicalClustering <- function(xdata, nc = NULL, Lambda,
                                         use_permutations = FALSE,
                                         linkage = "complete",
                                         scale = TRUE, rows = TRUE, ...) {
  # Checking sparcl package is installed
  if (!requireNamespace("sparcl")) {
    stop("This function requires the 'sparcl' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::hclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "wbound", "silent", "method")]

  # Calibrating Lambda
  if (use_permutations) {
    # Running automated calibration from sparcl
    sink(nullfile())
    myperm <- do.call(sparcl::HierarchicalSparseCluster.permute, args = c(
      list(x = xdata, wbound = Lambda[, 1]),
      tmp_extra_args
    ))
    sink()
    Lambda <- cbind(myperm$bestw)
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * nrow(Lambda)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Iterating over the pair of parameters
  id <- 0
  for (i in 1:nrow(Lambda)) {
    # Running sparse hierarchical clustering
    myclust <- do.call(sparcl::HierarchicalSparseCluster, args = c(
      list(x = xdata, wbound = Lambda[i, 1], silent = TRUE, method = linkage),
      tmp_extra_args
    ))

    # Defining clusters
    mygroups <- do.call(stats::cutree, args = list(tree = myclust$hc, k = nc))
    if (is.null(dim(mygroups))) {
      mygroups <- cbind(mygroups)
    }
    for (j in 1:nrow(nc)) {
      adjacency[, , id + j] <- CoMembership(groups = mygroups[, j])
      weight[id + j, ] <- myclust$ws[, 1]
    }
    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight, Lambda = Lambda))
}


SparseKMeansClustering <- function(xdata, nc = NULL, Lambda,
                                   use_permutations = FALSE,
                                   scale = TRUE, rows = TRUE, ...) {
  # Checking sparcl package is installed
  if (!requireNamespace("sparcl")) {
    stop("This function requires the 'sparcl' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::hclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "wbound", "silent", "K")]

  # Defining the number of iterations
  if (use_permutations) {
    n_iter_lambda <- 1
  } else {
    n_iter_lambda <- nrow(Lambda)
    Lambda_iter=Lambda
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * n_iter_lambda))
  weight <- matrix(NA, nrow = nrow(nc) * n_iter_lambda, ncol = ncol(xdata))

  # Iterating over the pair of parameters
  id <- 0
  for (i in 1:n_iter_lambda) {
    for (j in 1:nrow(nc)) {
#      print(nc[j])
      # Calibrating Lambda
      if (use_permutations) {
        # Running automated calibration from sparcl
        sink(nullfile())
        myperm <- tryCatch(do.call(sparcl::KMeansSparseCluster.permute, args = c(
          list(x = xdata, K = nc[j], wbound = Lambda[, 1]),
          tmp_extra_args
        )),
        error = function(e) {
          message("Not run.")
        }
        )
        sink()
        if (!is.null(myperm)){
          Lambda_iter <- cbind(myperm$bestw)
          
          # Running sparse K means clustering
          myclust <- tryCatch(
            do.call(sparcl::KMeansSparseCluster, args = c(
              list(x = xdata, K = nc[j], wbound = Lambda_iter, silent = TRUE),
              tmp_extra_args
            )),
            error = function(e) {
              message("Not run.")
            }
          )
        } else {
          myclust=NULL
        }
      } else {
        # Running sparse K means clustering
        myclust <- tryCatch(
          do.call(sparcl::KMeansSparseCluster, args = c(
            list(x = xdata, K = nc[j], wbound = Lambda[i, 1], silent = TRUE),
            tmp_extra_args
          )),
          error = function(e) {
            message("Not run.")
          }
        )
      }

      # Defining clusters
      if (!is.null(myclust)) {
        mygroups <- myclust[[1]]$Cs
        adjacency[, , id + j] <- CoMembership(groups = mygroups)
        weight[id + j, ] <- myclust[[1]]$ws
      } else {
        # Single cluster and unique weights in case of error
        mygroups <- rep(1, nrow(xdata))
        adjacency[, , id + j] <- CoMembership(groups = mygroups)
        weight[id + j, ] <- rep(1, ncol(xdata))
      }
    }
    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(
    comembership = adjacency,
    weight = weight,
    Lambda = Lambda_iter
  ))
}


# Code from clustvarsel package
# Fixing oneBIC$modelName error
clvarselhlfwd <- function(X, G = 1:9,
                          emModels1 = c("E", "V"),
                          emModels2 = mclust.options("emModelNames"),
                          samp = FALSE, sampsize = 2000,
                          hcModel = "VVV",
                          allow.EEE = TRUE, forcetwo = TRUE,
                          BIC.upper = 0, BIC.lower = -10,
                          itermax = 100,
                          verbose = interactive()) {
  X <- as.matrix(X)
  n <- nrow(X) # number of rows=number of observations
  d <- ncol(X) # number of columns=number of variables
  G <- setdiff(G, 1)

  # If needed, sample the subset of observations for hierarchical clustering
  if (samp) {
    sub <- sample(1:n, min(sampsize, n), replace = FALSE)
  } else {
    sub <- seq.int(1, n)
  }

  # First Step - selecting single variable and ordering
  if (verbose) cat(paste("iter 1\n+ adding step\n"))
  maxBIC <- BICdiff <- oneBIC <- rep(NA, d)
  ModelG <- vector(mode = "list", length = d)
  for (i in 1:d)
  {
    xBIC <- NULL
    # Fit the single variable cluster models
    try(
      xBIC <- Mclust(X[, i],
        G = G, modelNames = emModels1,
        initialization = list(subset = sub),
        verbose = FALSE
      ),
      silent = TRUE
    )
    if (is.null(xBIC)) {
      try(xBIC <- Mclust(X[, i], G = G, modelNames = emModels1),
        silent = TRUE
      )
    }
    # If we get all NA's from "V" starting hierarchical values use "E"
    if ((allow.EEE) & sum(is.finite(xBIC$BIC)) == 0) {
      try(
        xBIC <- Mclust(X[, i],
          G = G, modelNames = emModels1,
          initialization = list(
            hcPairs = hcE(X[sub, i]),
            subset = sub
          ),
          verbose = FALSE
        ),
        silent = TRUE
      )
    }
    # maxBIC is the maximum BIC over all clustering models fit
    if (sum(is.finite(xBIC$BIC)) == 0) {
      maxBIC[i] <- NA
    } else {
      maxBIC[i] <- max(xBIC$BIC[is.finite(xBIC$BIC)])
    }
    # Fit and get BIC for a single component no-cluster normal model
    try(
      oneBIC[i] <- Mclust(X[, i],
        G = 1, modelNames = emModels1,
        initialization = list(subset = sub),
        verbose = FALSE
      )$BIC[1],
      silent = TRUE
    )
    # Difference between maximum BIC for clustering and BIC for no clustering
    BICdiff[i] <- c(maxBIC[i] - oneBIC[i])
    ModelG[[i]] <- c(xBIC$modelName, xBIC$G)
  }

  # Find the single variable with the biggest difference between
  # clustering and no clustering
  m <- max(BICdiff[is.finite(BICdiff)])
  arg <- which(BICdiff == m, arr.ind = TRUE)[1]
  # This is our first selected variable/S is the matrix of currently selected
  # clustering variables
  S <- X[, arg, drop = FALSE]
  # BICS is the BIC value for the clustering model with the variable(s) in S
  BICS <- maxBIC[arg]
  temp <- order(BICdiff[-arg], decreasing = TRUE)
  # NS is the matrix of currently not selected variables
  NS <- as.matrix(X[, -arg])
  # This orders NS in terms of strongest evidence of univariate clustering
  # versus no clustering
  NS <- NS[, temp, drop = FALSE]
  # info records the proposed variable, BIC for the S matrix and difference
  # in BIC for clustering versus no clustering on S, whether it was an
  # addition step and if it was accepted
  info <- data.frame(
    Var = colnames(S),
    BIC = BICS, BICdiff = BICdiff[arg],
    Step = "Add", Decision = "Accepted",
    Model = ModelG[[arg]][1],
    G = ModelG[[arg]][2],
    stringsAsFactors = FALSE
  )
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)

  # Second Step - selecting second variable
  if (verbose) {
    print(info[, c(1, 3:5), drop = FALSE])
    cat(paste("iter 2\n+ adding step\n"))
  }

  depBIC <- cindepBIC <- cdiff <- rep(NA, ncol(NS))
  ModelG <- vector(mode = "list", length = ncol(NS))
  crit <- -Inf
  i <- 0
  # We only run until we find a variable whose difference in BIC between
  # being included in the clustering variables versus conditionally
  # independent of the clustering is greater than BIC.upper
  while ((crit <= BIC.upper) & (i < ncol(NS))) {
    i <- i + 1

    # Calculate the BIC for the regression of the proposed variable
    # on the variable in S
    regBIC <- BICreg(y = NS[, i], x = S)

    sBIC <- NULL
    # Fit the cluster model on the two variables
    try(
      sBIC <- Mclust(cbind(S, NS[, i]),
        G = G, modelNames = emModels2,
        initialization = list(
          hcPairs = hc(hcModel,
            data = cbind(S, NS[, i])[sub, ]
          ),
          subset = sub
        ),
        verbose = FALSE
      ),
      silent = TRUE
    )
    # If we get all NA's from "VVV" starting hierarchical values use "EEE"
    if ((allow.EEE) & sum(is.finite(sBIC$BIC)) == 0) {
      try(
        sBIC <- Mclust(cbind(S, NS[, i]),
          G = G, modelNames = emModels2,
          initialization = list(
            hcPairs = hc("EEE",
              data = cbind(S, NS[, i])[sub, ]
            ),
            subset = sub
          ),
          verbose = FALSE
        ),
        silent = TRUE
      )
    }
    # depBIC is the BIC for the clustering model with both variables
    if (sum(is.finite(sBIC$BIC)) > 0) {
      depBIC[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
    }
    # cindepBIC is the BIC for the clustering model on S and the regression
    # model of the new variable on S
    cindepBIC[i] <- regBIC + BICS
    cdiff[i] <- depBIC[i] - cindepBIC[i]
    if (!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
    crit <- cdiff[i]
    ModelG[[i]] <- c(sBIC$modelName, sBIC$G)
  }
  depBIC <- depBIC[1:i]
  cindepBIC <- cindepBIC[1:i]
  cdiff <- cdiff[1:i]

  # i is the index of those variables not selected but whose evidence of
  # clustering BIC did not fall below BIC.lower or those not looked at yet
  if (cdiff[i] > BIC.upper) {
    # i.e. evidence is stronger for including variable in S
    k <- c(colnames(S), colnames(NS)[i])
    S <- cbind(S, NS[, i])
    colnames(S) <- k
    BICS <- depBIC[i]
    info <- rbind(info, c(
      colnames(NS)[i], BICS, cdiff[i],
      "Add", "Accepted", ModelG[[i]]
    ))
    ns <- s <- NULL
    if (i < ncol(NS)) {
      ns <- seq(i + 1, ncol(NS))
    }
    if (i > 1) {
      s <- seq(i - 1)[which(cdiff[-i] > BIC.lower)]
    }
    ind <- c(s, ns)
    if (!is.null(ind)) {
      nks <- c(colnames(NS)[ind])
      # NS is the not selected clustering variables whose recently calculated
      # evidence of clustering BIC was higher than BIC.lower or variables not
      # yet looked at
      NS <- as.matrix(NS[, ind])
      colnames(NS) <- nks
    } else {
      NS <- NULL
    }
  } else {
    if ((cdiff[i] < BIC.upper) & (forcetwo)) {
      # if the evidence is weaker but we're forcing choice of second variable
      m <- max(cdiff[is.finite(cdiff)])
      i <- which(cdiff == m, arr.ind = TRUE)[1]
      k <- c(colnames(S), colnames(NS)[i])
      S <- cbind(S, NS[, i])
      colnames(S) <- k
      BICS <- depBIC[i]
      info <- rbind(info, c(
        colnames(NS)[i], BICS, cdiff[i],
        "Add", "Accepted", ModelG[[i]]
      ))
      nks <- c(colnames(NS)[-i])
      # NS is the not selected clustering variables whose recently calculated
      # evidence of clustering BIC was higher than BIC.lower
      NS <- as.matrix(NS[, -i])
      temp <- cdiff[-i]
      if (sum(temp > BIC.lower) != 0) {
        NS <- as.matrix(NS[, c(which(temp > BIC.lower))])
        colnames(NS) <- nks[c(which(temp > BIC.lower))]
      } else {
        NS <- NULL
      }
    } else {
      m <- max(cdiff[is.finite(cdiff)])
      i <- which(cdiff == m, arr.ind = TRUE)[1]
      info <- rbind(info, c(
        colnames(NS)[i], BICS, cdiff[i],
        "Add", "Rejected", ModelG[[i]]
      ))
    }
  }
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)

  if (verbose) print(info[2, c(1, 3:5), drop = FALSE])

  criterion <- 1
  iter <- 0
  while ((criterion == 1) & (iter < itermax)) {
    iter <- iter + 1
    check1 <- colnames(S)

    if (verbose) cat(paste("iter", iter + 2, "\n"))

    # Addition step
    if (verbose) cat("+ adding step\n")
    # For the special case where we have removed all the clustering variables/S
    # is empty [LS: is really needed??]
    if ((NCOL(NS) != 0 & !is.null(ncol(NS))) &
      (ncol(S) == 0) || (is.null(ncol(S)))) { # cat("\n really needed ??\n")
      depBIC <- 0
      DepBIC <- NULL
      crit <- -10
      cdiff <- 0
      Cdiff <- NULL
      oneBIC <- rep(NA, d)
      ModelG <- vector(mode = "list", length = d)
      i <- 0
      crit <- -10
      while ((crit <= BIC.upper) & (i < ncol(NS))) {
        xBIC <- NULL
        i <- i + 1
        # Fit the cluster models
        try(
          xBIC <- Mclust(X[, i],
            G = G, modelNames = emModels1,
            initialization = list(subset = sub),
            verbose = FALSE
          ),
          silent = TRUE
        )
        # If we get all NA's from "V" starting hierarchical values use "E"
        if ((allow.EEE) & sum(is.finite(xBIC$BIC)) == 0) {
          try(
            xBIC <- Mclust(X[, i],
              G = G, modelNames = emModels1,
              initialization = list(
                hcPairs = hcE(X[sub, i]),
                subset = sub
              ),
              verbose = FALSE
            ),
            silent = TRUE
          )
        }
        # depBIC is the maximum BIC over all clustering models fit
        if (sum(is.finite(xBIC$BIC)) == 0) {
          depBIC <- NA
        } else {
          depBIC <- max(xBIC$BIC[is.finite(xBIC$BIC)])
        }
        DepBIC <- c(DepBIC, depBIC)
        # Fit and get BIC for a single component no-cluster normal model
        try(oneBIC <- Mclust(X[, i], G = 1, modelNames = "V", verbose = FALSE)$BIC[1],
          silent = TRUE
        )
        # Difference between maximum BIC for clustering and BIC for no clustering
        cdiff <- c(depBIC - oneBIC)
        if (!is.finite(cdiff)) cdiff <- BIC.upper
        Cdiff <- c(Cdiff, cdiff)
        crit <- cdiff
        ModelG[[i]] <- c(xBIC$modelName, xBIC$G)
      }
      if (cdiff > BIC.upper) { # ie. evidence is stronger for including variable in S
        k <- c(colnames(NS)[i])
        S <- as.matrix(NS[, i])
        colnames(S) <- k
        BICS <- depBIC
        info <- rbind(info, c(
          colnames(NS)[i], BICS, cdiff,
          "Add", "Accepted", ModelG[[i]]
        ))
        ns <- s <- NULL
        # i is the index of those variables not selected but whose evidence of
        # clustering BIC did not fall below BIC.lower or those not looked at yet
        if (i < ncol(NS)) ns <- seq(i + 1, ncol(NS))
        if (i > 1) s <- seq(i - 1)[which(Cdiff[-i] > BIC.lower)]
        ind <- c(s, ns)
        if (!is.null(ind)) {
          nks <- c(colnames(NS)[ind])
          # NS is the not selected clustering variables whose recently
          # calculated evidence of clustering BIC was higher than BIC.lower
          # or variables not yet looked at
          NS <- as.matrix(NS[, ind])
          colnames(NS) <- nks
        } else {
          NS <- NULL
        }
      } else {
        m <- max(Cdiff[is.finite(Cdiff)])
        i <- which(Cdiff == m, arr.ind = TRUE)[1]
        info <- rbind(info, c(
          colnames(NS)[i], BICS, Cdiff[i],
          "Add", "Rejected", ModelG[[i]]
        ))
        ind <- seq(ncol(NS))[which(Cdiff > BIC.lower)]
        if (!is.null(ind)) {
          k <- colnames(NS)[ind]
          # Exclude variables in NS whose evidence of clustering in this step
          # was lower than BIC.lower
          NS <- as.matrix(NS[, ind])
          colnames(NS) <- k
        } else {
          NS <- NULL
        }
      }
    }
    # [LS: is really needed?? -- END]
    else {
      # Addition Step in general (for all cases except when S is empty)
      if ((NCOL(NS) != 0) & !is.null(ncol(NS))) {
        depBIC <- cindepBIC <- cdiff <- rep(NA, ncol(NS))
        ModelG <- vector(mode = "list", length = ncol(NS))
        crit <- -Inf
        i <- 0
        # We only run until we find a variable whose difference in BIC
        # between being included in the clustering variables versus
        # conditionally independent of the clustering is greater than BIC.upper
        while (crit <= BIC.upper & i < ncol(NS)) {
          sBIC <- NULL
          i <- i + 1
          # Calculate the BIC for the regression of the proposed variable
          # on the variable(s) in S
          regBIC <- BICreg(y = NS[, i], x = S)

          # Fit the cluster model on the S variables with the proposed variable
          try(
            sBIC <- Mclust(cbind(S, NS[, i]),
              G = G, modelNames = emModels2,
              initialization = list(
                hcPairs = hc(hcModel,
                  data = cbind(S, NS[, i])[sub, ]
                ),
                subset = sub
              ),
              verbose = FALSE
            ),
            silent = TRUE
          )
          # If we get all NA's from "VVV" starting hierarchical values use "EEE"
          if ((allow.EEE) & (sum(is.finite(sBIC$BIC)) == 0)) {
            try(
              sBIC <- Mclust(cbind(S, NS[, i]),
                G = G, modelNames = emModels2,
                initialization = list(
                  hcPairs = hc("EEE",
                    data = cbind(S, NS[, i])[sub, ]
                  ),
                  subset = sub
                ),
                verbose = FALSE
              ),
              silent = TRUE
            )
          }
          # depBIC is the BIC for the clustering model with both S and proposed
          # variable
          if (sum(is.finite(sBIC$BIC)) > 0) {
            depBIC[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
          }
          # cindepBIC is the BIC for the clustering model on S and the
          # regression model of the new variable on S
          cindepBIC[i] <- regBIC + BICS
          cdiff[i] <- depBIC[i] - cindepBIC[i]
          if (!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
          crit <- cdiff[i]
          ModelG[[i]] <- c(sBIC$modelName, sBIC$G)
        }
        depBIC <- depBIC[1:i]
        cindepBIC <- cindepBIC[1:i]
        cdiff <- cdiff[1:i]
        if (cdiff[i] > BIC.upper) {
          # i.e. evidence is stronger for including variable in S
          k <- c(colnames(S), colnames(NS)[i])
          nks <- c(colnames(NS)[-i])
          S <- cbind(S, NS[, i])
          colnames(S) <- k
          BICS <- depBIC[i]
          info <- rbind(info, c(
            colnames(NS)[i], BICS, cdiff[i],
            "Add", "Accepted", ModelG[[i]]
          ))
          ns <- s <- NULL
          # Exclude variables in NS whose evidence of clustering in this step
          # was lower than BIC.lower
          if (i < ncol(NS)) ns <- seq(i + 1, ncol(NS))
          if (i > 1) s <- seq(i - 1)[which(cdiff[-i] > BIC.lower)]
          ind <- c(s, ns)
          if (!is.null(ind)) {
            nks <- colnames(NS)[ind]
            NS <- as.matrix(NS[, ind])
            colnames(NS) <- nks
          } else {
            NS <- NULL
          }
        } else {
          m <- max(cdiff[is.finite(cdiff)])
          i <- which(cdiff == m, arr.ind = TRUE)[1]
          info <- rbind(info, c(
            colnames(NS)[i], depBIC[i], cdiff[i],
            "Add", "Rejected", ModelG[[i]]
          ))
          ind <- seq(1, ncol(NS))[which(cdiff > BIC.lower)]
          if (!is.null(ind)) {
            k <- colnames(NS)[ind]
            # Exclude variables in NS whose evidence of clustering in this
            # step was lower than BIC.lower
            NS <- as.matrix(NS[, ind])
            colnames(NS) <- k
          } else {
            NS <- NULL
          }
        }
      }
    }

    # Removal Step for the special case where S contains only a single variable
    if (verbose) cat("- removing step\n")
    if (ncol(S) == 1) {
      cdiff <- 0
      oneBIC <- NA
      try(
        oneBIC <- Mclust(S,
          G = 1, modelNames = "V",
          initialization = list(subset = sub),
          verbose = FALSE
        )$BIC[1],
        silent = TRUE
      )
      # Difference between maximum BIC for clustering and BIC for no clustering
      cdiff <- BICS - oneBIC
      if (is.na(cdiff)) cdiff <- BIC.upper
      # check if difference is negative
      if (cdiff <= BIC.upper) {
        # if negative remove the variable from S and set the BIC for
        # the model to NA
        BICS <- NA
        info <- rbind(info, c(
          colnames(S), BICS, cdiff,
          "Remove", "Accepted",
          # oneBIC$modelName, oneBIC$G))
          "V", "G"
        )) # BB: removed oneBIC$modelName and oneBIC$G
        # Only return variable to NS if difference is greater than BIC.lower
        if (cdiff > BIC.lower) {
          k <- c(colnames(NS), colnames(S))
          NS <- cbind(NS, S)
          colnames(NS) <- k
          S <- NULL
        } else {
          S <- NULL
        }
      } else {
        info <- rbind(info, c(
          colnames(S), BICS, cdiff,
          "Remove", "Rejected",
          # oneBIC$modelName, oneBIC$G))
          "V", "G"
        )) # BB: removed oneBIC$modelName and oneBIC$G
      }
    } else {
      # Removal step in general (for all cases except when S is a single
      # variable or empty)
      if (ncol(S) >= 2) {
        depBIC <- BICS
        cindepBIC <- rdep <- cdiff <- rep(NA, ncol(S))
        ModelG <- vector(mode = "list", length = ncol(S))
        crit <- Inf
        i <- 0
        # Check if the data is at least 3 dimensional
        name <- if (ncol(S) > 2) emModels2 else emModels1
        # We only run until we find a variable whose difference in BIC between
        # being included in the clustering variables versus conditionally
        # independent of the clustering is lower than BIC.upper
        while (crit > BIC.upper & (i < ncol(S))) {
          i <- i + 1
          # Calculate the BIC for the regression of the proposed variable
          # from S on the other variable(s) in S
          regBIC <- BICreg(y = S[, i], x = S[, -i])
          # Fit the cluster model on the S variables without the
          # proposed variable
          sBIC <- NULL
          try(
            sBIC <- Mclust(S[, -i],
              G = G, modelNames = name,
              initialization =
                list(
                  hcPairs = hc(hcModel,
                    data = S[sub, -i, drop = FALSE]
                  ),
                  subset = sub
                ),
              verbose = FALSE
            ),
            silent = TRUE
          )
          # If we get all NA's from "VVV" starting hierarchical values use "EEE"
          if (allow.EEE & (ncol(S) >= 3) & sum(is.finite(sBIC$BIC)) == 0) {
            try(
              sBIC <- Mclust(S[, -i],
                G = G, modelNames = name,
                initialization = list(
                  hcPairs = hc("EEE",
                    data = S[sub, -i]
                  ),
                  subset = sub
                ),
                verbose = FALSE
              ),
              silent = TRUE
            )
          } else {
            if ((allow.EEE) & (ncol(S) == 2) & sum(is.finite(sBIC$BIC)) == 0) {
              try(
                sBIC <- Mclust(as.matrix(S[, -i]),
                  G = G, modelNames = name,
                  initialization = list(
                    hcPairs = hcE(S[sub, -i, drop = FALSE]),
                    subset = sub
                  ),
                  verbose = FALSE
                ),
                silent = TRUE
              )
            }
          }
          if (sum(is.finite(sBIC$BIC)) > 0) {
            rdep[i] <- max(sBIC$BIC[is.finite(sBIC$BIC)])
          }
          # cindepBIC is the BIC for the clustering model on the other variables
          # in S and the regression model of the proposed variable on the other
          # variables in S
          cindepBIC[i] <- regBIC + rdep[i]
          cdiff[i] <- depBIC - cindepBIC[i]
          if (!is.finite(cdiff[i])) cdiff[i] <- BIC.upper
          crit <- cdiff[i]
          ModelG[[i]] <- c(sBIC$modelName, sBIC$G)
        }
        if ((cdiff[i] < BIC.upper) & (cdiff[i] > BIC.lower)) { # i.e. evidence is stronger for excluding variable from S but
          # still including it in NS
          BICS <- rdep[i]
          info <- rbind(info, c(
            colnames(S)[i], BICS, cdiff[i],
            "Remove", "Accepted", ModelG[[i]]
          ))
          k <- c(colnames(NS), colnames(S)[i])
          nk <- colnames(S)[-i]
          NS <- cbind(NS, S[, i])
          S <- as.matrix(S[, -i])
          colnames(NS) <- k
          colnames(S) <- nk
        } else {
          if (cdiff[i] < BIC.lower) { # exclude variable entirely
            BICS <- rdep[i]
            info <- rbind(info, c(
              colnames(S)[i], BICS, cdiff[i],
              "Remove", "Accepted", ModelG[[i]]
            ))
            nk <- colnames(S)[-i]
            S <- as.matrix(S[, -i])
            colnames(S) <- nk
          } else {
            m <- min(cdiff[is.finite(cdiff)])
            i <- which(cdiff == m, arr.ind = TRUE)[1]
            info <- rbind(info, c(
              colnames(S)[i], rdep[i], cdiff[i],
              "Remove", "Rejected", ModelG[[i]]
            ))
          }
        }
      }
    }
    info$BIC <- as.numeric(info$BIC)
    info$BICdiff <- as.numeric(info$BICdiff)

    if (verbose) {
      print(info[seq(nrow(info) - 1, nrow(info)), c(1, 3:5), drop = FALSE])
    }

    # Check if the variables in S have changed or not
    check2 <- colnames(S)
    if (is.null(check2)) # all variables have been removed
      {
        criterion <- 0
      } else
    # if they have changed (either added one or removed one or changed one)
    # then continue the algorithm (criterion is 1) otherwise stop
    # (criterion is 0)
    {
      if (length(check2) != length(check1)) {
        criterion <- 1
      } else {
        criterion <- if (sum(check1 == check2) != length(check1)) 1 else 0
      }
    }
  }

  if (iter >= itermax) {
    warning("Algorithm stopped because maximum number of iterations was reached")
  }

  # List the selected variables and the matrix of steps' information
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)
  # reorder steps.info
  info <- info[, c(1, 4, 2, 6, 7, 3, 5), drop = FALSE]
  colnames(info) <- c(
    "Variable proposed", "Type of step",
    "BICclust", "Model", "G", "BICdiff", "Decision"
  )
  varnames <- colnames(X)
  subset <- sapply(colnames(S), function(x) which(x == varnames))

  if (length(subset) == 0) {
    subset <- varnames # BB: adding to avoid error
  }

  out <- list(
    variables = varnames,
    subset = subset,
    steps.info = info,
    search = "headlong",
    direction = "forward"
  )

  return(out)
}
