## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(M3C)
# now we have loaded the mydata and desx objects (with the package automatically)
# mydata is the expression data for GBM
# desx is the annotation for this data

## ----fig.width=4,fig.height=4-------------------------------------------------
pca(mydata,legendtextsize = 10,axistextsize = 10,dotsize=2)

## ---- results='hide'----------------------------------------------------------
# for vignette
res <- M3C(mydata, des = desx, removeplots = TRUE, iters=25, 
           objective='PAC', fsize=8, lthick=1, dotsize=1.25)

## -----------------------------------------------------------------------------
res$scores

## -----------------------------------------------------------------------------
for (k in seq(2,10)){
  myresults <- res$realdataresults[[k]]$ordered_annotation
  chifit <- suppressWarnings(chisq.test(table(myresults[c('consensuscluster','class')])))
  print(chifit$p.value)
}

## ----fig.width=5,fig.height=3.5-----------------------------------------------
res$plots[[1]]

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[2]]

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[4]]

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[3]]

## ----fig.show = 'hide'--------------------------------------------------------
data <- res$realdataresults[[4]]$ordered_data 
annon <- res$realdataresults[[4]]$ordered_annotation 
ccmatrix <- res$realdataresults[[4]]$consensus_matrix
head(annon)

## -----------------------------------------------------------------------------
# library(ComplexHeatmap)
# ccl <- list()
# x <- c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen", 
#                 "violetred", "orange", "midnightblue", "grey31", "black")
# names(x) <- as.character(seq(1,11,by=1))
# for (i in seq(2,10)){
#   # get cc matrix and labels
#   ccmatrix <- res$realdataresults[[i]]$consensus_matrix
#   annon <- res$realdataresults[[i]]$ordered_annotation
#   # do heatmap
#   n <- 10
#   seq <- rev(seq(0,255,by=255/(n)))
#   palRGB <- cbind(seq,seq,255)
#   mypal <- rgb(palRGB,maxColorValue=255)
#   ha = HeatmapAnnotation(
#   df= data.frame(Cluster=as.character(annon[,1])), col = list(Cluster=x))
#   ccl[[i]] <-  Heatmap(ccmatrix, name = "Consensus_index", top_annotation = ha, 
#                col=mypal, show_row_dend = FALSE,
#                show_column_dend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
#                show_column_names = FALSE)
# }

## ----fig.width=5,fig.height=4-------------------------------------------------
pca(data,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)

## ----fig.show = 'hide',results='hide'-----------------------------------------
res <- M3C(mydata, method = 2, objective='PAC',fsize=8, lthick=1, dotsize=1.25)

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[3]]

## ----fig.show = 'hide',results='hide'-----------------------------------------
res <- M3C(mydata, method = 2,fsize=8, lthick=1, dotsize=1.25)

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[2]]

## ----fig.width=4,fig.height=3-------------------------------------------------
res$plots[[3]]

## -----------------------------------------------------------------------------
  filtered_results <- featurefilter(mydata, percentile=10, method='MAD', topN=5)

