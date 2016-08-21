#'
#' Computes bootstrap support for shift positions
#'
#' Takes a given shift configuration previously detected from data along with shift magnitudes
#' and OU parameters, to calculate bootstrap support for shift positions. 
#' The non-parametric bootstrap procedure calculates phylogenetically-uncorrelated standardized residuals,
#' one at each node. These residuals are sampled with replacement, then mapped back onto the tree
#' to create bootstrap replicates. Each replicate is analyzed with the l1ou method and user-specified options.
#'
#'@param model an object output by \code{\link{estimate_shift_configuration}}. 
#'@param nItrs number of independent iterations (bootstrap replicates).
#'@param multicore logical. If TRUE, nCores processes are used in parallel. 
#'@param nCores desired number of parallel processes.
#'@param quietly logical. If FALSE, a summary of each iteration will be printed out.
#'
#'
#'@return vector of size the number of edges in the tree. Each entry is the proportion of bootstrap replicates for which a shift is detected on the corresponding edge. 
#'
#'
#'@details The results of sequential and parallel runs are not necessarily equal, because different seeds might be used for different bootstrap replicates.
#'         For multiple cores to be used, the \code{parallel} library needs to be installed.
#'         To change options for the analysis of each bootstrap replicate,
#'         like the information criterion or the maximum allowed number of shifts, modify model$opt.
#'
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1] 
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' result <- l1ou_bootstrap_support(eModel, nItrs=2)
#' # using only 2 replicates in vastly insufficient in general,
#' # but used here to make the illustrative example run faster.
#' nEdges <- Nedge(lizard.tree)
#' e.w <- rep(1,nEdges) 
#' e.w[eModel$shift.configuration] <- 3
#' e.l <- round(result$detection.rate*100, digits=1)
#' # to avoid annotating edges with support at or below 10%
#' e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
#' plot(eModel, edge.label=e.l, edge.ann.cex=0.7, edge.label.ann=TRUE, cex=0.5, label.offset=0.02, edge.width=e.w)
#'
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
l1ou_bootstrap_support <- function(model, nItrs=100, multicore=FALSE, nCores = 2, quietly=TRUE){
 
    if (!inherits(model, "l1ou"))  stop("object \"model\" is not of class \"l1ou\".")

    if(multicore)
        multicore = require("parallel")

    tree = model$tree
    if(ncol(model$Y)==1){
        return(bootstrap_support_univariate(tree=tree, model=model, nItrs=nItrs, 
                                            multicore=multicore, nCores=nCores, quietly=quietly))
    }
    if(ncol(model$Y)>1){
        return(bootstrap_support_multivariate(tree=tree, model=model, nItrs=nItrs, 
                                              multicore=multicore, nCores=nCores, quietly=quietly))
    }
}

bootstrap_support_univariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2, quietly=FALSE){

    RE    = sqrt_OU_covariance(tree, alpha=model$alpha, 
                               root.model = model$l1ou.options$root.model,
                               check.order=F, check.ultrametric=F)

    C.IH  = t(RE$sqrtInvSigma)
    C.H   = RE$sqrtSigma

    Y     = model$Y
    YY    = C.IH%*%(Y - model$mu )

    seed.vec <- sample(.Machine$integer.max, nItrs+1, replace=TRUE)

    detection.vec = rep(0, nrow(tree$edge))
    all.shift.configurations <- list()

    if(quietly==FALSE)
        print(paste0("iteration #:nShifts:shift configuraitons"))

    valid.count <- 0
    if(multicore == FALSE){
        for(itr in 1:nItrs){
            set.seed(seed.vec[[itr]])
            YYstar = sample(YY, replace = TRUE)
            Ystar  = as.matrix( (C.H%*%YYstar) + model$mu )
            rownames(Ystar) <- rownames(Y)

            eM  <-  tryCatch({
                estimate_shift_configuration(tree, Ystar, l1ou.options =model$l1ou.options)
            }, error = function(e) {
                print("l1OU error, return NA")
                return(NA) }  )
            if(all(is.na(eM))) {return(NA)}

            valid.count <- valid.count + 1

            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
            all.shift.configurations[[itr]] <- eM$shift.configuration

            if(quietly==FALSE){
                print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                             paste0(eM$shift.configuration, collapse=" ") ) )
            }

        }
        set.seed(seed.vec[[nItrs+1]])
        return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
    }


    all.shift.configurations = 
        mclapply(X=1:nItrs, FUN=function(itr){

                     set.seed(seed.vec[[itr]])
                     YYstar = sample(YY, replace = TRUE)
                     Ystar  = as.matrix( (C.H%*%YYstar) + model$mu )
                     rownames(Ystar) <- rownames(Y)

                     eM  <-  tryCatch({
                         estimate_shift_configuration(tree, Ystar, l1ou.options =model$l1ou.options)
                     }, error = function(e) {
                         print("l1OU error, return NA")
                         return(NA) }  )

                     if(all(is.na(eM))) {return(NA)}

                     if(quietly==FALSE){
                         print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                                      paste0(eM$shift.configuration, collapse=" ") ) )
                     }
                     return(eM$shift.configuration)
           }, mc.cores = nCores)

    valid.count <- 0
    na.indices <- c()
    for( i in 1:length(all.shift.configurations)){
        if( all(is.na( all.shift.configurations[[i]] )) ){
            na.indices <- c(na.indices, i)
            next
        }
        valid.count <- valid.count + 1
        detection.vec[ all.shift.configurations[[i]] ] = 
            detection.vec[ all.shift.configurations[[i]] ] + 1
    }

    if(length(na.indices)>0)
        all.shift.configurations <- all.shift.configurations[-1*na.indices]

    set.seed(seed.vec[[nItrs+1]])
    return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
}

bootstrap_support_multivariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2, quietly=FALSE){

    Y = as.matrix(model$Y)
    stopifnot( length(model$alpha) == ncol(Y) )

    seed.vec <- sample(.Machine$integer.max, nItrs+1, replace=TRUE)

    YY        = Y
    C.Hlist   = list()
    for( idx in 1:ncol(Y) ){
        RE    = sqrt_OU_covariance(tree, alpha = model$alpha[[idx]], 
                                   root.model = model$l1ou.options$root.model,
                                   check.order=F, check.ultrametric=F ) 
        C.IH  = t(RE$sqrtInvSigma) 
        C.Hlist[[idx]] = RE$sqrtSigma
        YY[, idx]      = C.IH%*%(Y[, idx] - model$mu[ ,idx])

    }
    if(quietly==FALSE)
        print(paste0("iteration #:nShifts:shift configuraitons"))

    detection.vec = rep(0, nrow(tree$edge))
    all.shift.configurations <- list()

    valid.count <- 0
    if( multicore == FALSE ){
        for(itr in 1:nItrs){

            set.seed(seed.vec[[itr]])
            Ystar   = YY
            idx.vec = sample(1:nrow(YY), replace = TRUE)
            for( idx in 1:ncol(YY) ){
                YYstar        = YY[idx.vec, idx]
                Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx] 
            }
            rownames(Ystar) <- rownames(Y)

            eM  <-  tryCatch({
                estimate_shift_configuration(tree, Ystar,  l1ou.options=model$l1ou.options)
            }, error = function(e) {
                print("l1OU error, return NA")
                return(NA) }  )

            if(all(is.na(eM))) {next}

            if(quietly==FALSE){
                print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                             paste0(eM$shift.configuration, collapse=" ") ) )
            }

            valid.count  <- valid.count + 1
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
            all.shift.configurations[[itr]] <- eM$shift.configuration
        }
        stopifnot( valid.count > 0 )
        set.seed(seed.vec[[nItrs+1]])
        return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
    }

    all.shift.configurations = 
        mclapply(X=1:nItrs, FUN=function(itr){
                     Ystar   = YY
                     set.seed(seed.vec[[itr]])
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

                     if(quietly==FALSE){
                         print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                                      paste0(eM$shift.configuration, collapse=" ") ) )
                     }

                     return(eM$shift.configuration)
                }, mc.cores = nCores)

    valid.count <- 0
    na.indices <- c()
    for( i in 1:length(all.shift.configurations )){
        if( all(is.na( all.shift.configurations [[i]] )) ){
            na.indices <- c(na.indices, i)
            next
        }
        valid.count <- valid.count + 1

        stopifnot( valid.count > 0 )
        detection.vec[ all.shift.configurations [[i]] ] = 
            detection.vec[ all.shift.configurations [[i]] ] + 1
    }
    if(length(na.indices)>0)
        all.shift.configurations <- all.shift.configurations[-1*na.indices]

    set.seed(seed.vec[[nItrs+1]]) ## To make sure after both mclapply and for-loop we have same seed for the reproducibility  
    return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
}

