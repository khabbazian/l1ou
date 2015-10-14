#'
#' Computes bootstrap support for shift positions
#'
#' Takes a given shift configuration previously detected from data along with shift magnitudes
#' and OU parameters, to calculate bootstrap support for shift positions. 
#' The non-parametric bootstrap procedure calculates phylogenetically-uncorrelated standardized residuals,
#' one at each node. These residuals are sampled with replacement, then mapped back onto the tree
#' to create bootstrap replicates. Each replicate is analyzed with the l1ou method and user-specified options.
#'
#'@param tree ultrametric phylogenetic tree of class phylo with branch lengths.
#'@param model an object output by \code{\link{estimate_shift_configuration}}. 
#'@param nItrs number of independent iterations (bootstrap replicates).
#'@param multicore logical. If TRUE, nCores processes are used in parallel. 
#'@param nCores desired number of parallel processes.
#'@return vector of size the number of edges in the tree. Each entry is the proportion of bootstrap replicates for which a shift is detected on the corresponding edge. 
#'
#'@details The results of sequential and parallel runs are not necessarily equal, because different seeds might be used for different bootstrap replicates.
#'         For multiple cores to be used, the \code{parallel} library needs to be installed.
#'         To change options for the analysis of each bootstrap replicate,
#'         like the information criterion or the maximum allowed number of shifts, modify model$opt.
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1] 
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' result <- l1ou_bootstrap_support(lizard.tree, eModel, nItrs=2)
#' # using only 2 replicates in vastly insufficient in general,
#' # but used here to make the illustrative example run faster.
#' nEdges <- length(lizard.tree$edge[,1])
#' ew <- rep(1,nEdges) 
#' ew[eModel$shift.configuration] <- 3
#' lizard.tree$edge.label <- round(result * 100, digits=1)
#' # below: code to avoid annotating edges with support at or below 10%
#' lizard.tree$edge.label <- ifelse(lizard.tree$edge.label>10,
#'                                  paste0(lizard.tree$edge.label,"%"), NA)
#' plot_l1ou(lizard.tree, eModel, edge.ann.cex=0.7, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
l1ou_bootstrap_support <- function(tree, model, nItrs=100, multicore=FALSE, nCores = 2){

    if(multicore)
        multicore = require("parallel")

    if(ncol(model$Y)==1){
        return(bootstrap_support_univariate(tree=tree, model=model, nItrs=nItrs, multicore=multicore, nCores=nCores))
    }
    if(ncol(model$Y)>1){
        return(bootstrap_support_multivariate(tree=tree, model=model, nItrs=nItrs, multicore=multicore, nCores=nCores))
    }
}

bootstrap_support_univariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2){

    RE    = sqrt_OU_covariance(tree, alpha=model$alpha)

    C.IH  = t(RE$sqrtInvSigma)
    C.H   = RE$sqrtSigma

    Y     = model$Y
    YY    = C.IH%*%(Y - model$mu )

    detection.vec = rep(0, nrow(tree$edge))

    if(multicore == FALSE){
        for(itr in 1:nItrs){
            YYstar = sample(YY, replace = TRUE)
            Ystar  = (C.H%*%YYstar) + model$mu 
            eM     = estimate_shift_configuration(tree, Ystar, l1ou.options = model$l1ou.options)
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
        }
        return(detection.vec/nItrs)
    }

    shift.configuration.list = 
        mclapply(X=1:nItrs, FUN=function(itr){

                     set.seed( 101 + itr)
                     YYstar = sample(YY, replace = TRUE)
                     Ystar  = (C.H%*%YYstar) + model$mu  

                     eM  <-  tryCatch({
                         estimate_shift_configuration(tree, Ystar, l1ou.options =model$l1ou.options)
                     }, error = function(e) {
                         print("l1OU error, return NA")
                         return(NA) }  )

                     if(all(is.na(eM))) {return(NA)}
                     return(eM$shift.configuration)
           }, mc.cores = nCores)

    valid.count <- 0
    for( i in 1:length(shift.configuration.list)){
        if( all(is.na( shift.configuration.list[[i]] )) ){
            next
        }
        valid.count <- valid.count + 1
        detection.vec[ shift.configuration.list[[i]] ] = 
            detection.vec[ shift.configuration.list[[i]] ] + 1
    }

    return(detection.vec/valid.count)
}

bootstrap_support_multivariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2){

    Y = as.matrix(model$Y)
    stopifnot( length(model$alpha) == ncol(Y) )

    YY        = Y
    C.Hlist   = list()
    for( idx in 1:ncol(Y) ){
        RE    = sqrt_OU_covariance(tree, alpha = model$alpha[[idx]] ) 
        C.IH  = t(RE$sqrtInvSigma) 
        C.Hlist[[idx]] = RE$sqrtSigma
        YY[, idx]      = C.IH%*%(Y[, idx] - model$mu[ ,idx])
    }

    detection.vec = rep(0, nrow(tree$edge))

    if( multicore == FALSE ){
        for(itr in 1:nItrs){

            Ystar   = YY
            idx.vec = sample(1:nrow(YY), replace = TRUE)
            for( idx in 1:ncol(YY) ){
                YYstar        = YY[idx.vec, idx]
                Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx] 
            }
            eM  <-  tryCatch({
                estimate_shift_configuration(tree, Ystar,  l1ou.options=model$l1ou.options)
            }, error = function(e) {
                print("l1OU error, return NA")
                return(NA) }  )

            if(all(is.na(eM))) {next}
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
        }
    }

    shift.configuration.list = 
        mclapply(X=1:nItrs, FUN=function(itr){
                     Ystar   = YY
                     set.seed( 101 + itr)
                     idx.vec = sample(1:nrow(YY), replace = TRUE)
                     for( idx in 1:ncol(YY) ){
                         YYstar        = YY[idx.vec, idx]
                         Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx] 
                     }
                     eM  <-  tryCatch({
                         estimate_shift_configuration(tree, Ystar, l1ou.options = model$l1ou.options)
                     }, error = function(e) {
                         print("l1OU error, return NA")
                         return(NA) }  )

                     if(all(is.na(eM))) {return(NA)}
                     return(eM$shift.configuration)
                }, mc.cores = nCores)

    valid.count <- 0
    for( i in 1:length(shift.configuration.list)){
        if( all(is.na( shift.configuration.list[[i]] )) ){
            next
        }
        valid.count <- valid.count + 1
        detection.vec[ shift.configuration.list[[i]] ] = 
            detection.vec[ shift.configuration.list[[i]] ] + 1
    }

    return(detection.vec/valid.count)
}
