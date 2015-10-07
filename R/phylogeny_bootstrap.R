#'
#' computes the bootstrap support for the detected shift configuration 
#'
#'@param tr an ultrametric phylogenetic tree of type phylo with branch lengths and tip labels.
#'@param model the object generetes by function estimate_shift_configuration. 
#'@param nItrs number of independent iterations (bootstrap independent replicates).
#'@param multicore logical. If TRUE, it runs nCores processes in parallel. See details. 
#'@param nCores desired number of parallel processes.
#'@return vector of size of the number of edges in the tree. Each entry is the proportion of bootstrap replicates for which a shift is detected on the corresponding edge. 
#'
#'@details The results of sequential and parallel runs are not necessarily equal, because different seeds might be used for different bootstrap replicates.
#'         To change options, like the information criterion or maximum allowed number of shifts, modify model$opt.
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree);
#' Y <- lizard.traits[,1]; 
#' eModel <- estimate_shift_configuration(lizard.tree, Y);
#' result <- l1ou_bootstrap_support(lizard.tree, eModel, nItrs=2);
#' print(result);
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
l1ou_bootstrap_support <- function(tr, model, nItrs=100, multicore=FALSE, nCores = 2){

    if(multicore)
        multicore = require("parallel");

    if(ncol(model$Y)==1){
        return(bootstrap_support_univariate(tr=tr, model=model, nItrs=nItrs, multicore=multicore, nCores=nCores));
    }
    if(ncol(model$Y)>1){
        return(bootstrap_support_multivariate(tr=tr, model=model, nItrs=nItrs, multicore=multicore, nCores=nCores));
    }
}

bootstrap_support_univariate <- function(tr, model, nItrs, multicore=FALSE, nCores=2){

    RE    = sqrt_OU_covariance(tr, alpha=model$alpha);

    C.IH  = t(RE$sqrtInvSigma);
    C.H   = RE$sqrtSigma;

    Y     = model$Y;
    YY    = C.IH%*%(Y - model$mu );

    detection.vec = rep(0, nrow(tr$edge));

    if(multicore == FALSE){
        for(itr in 1:nItrs){
            YYstar = sample(YY, replace = TRUE);
            Ystar  = (C.H%*%YYstar) + model$mu; 
            eM     = estimate_shift_configuration(tr, Ystar, l1ou.options = model$l1ou.options);
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1;
        }
        return(detection.vec/nItrs);
    }

    shift.configuration.list = 
        mclapply(X=1:nItrs, FUN=function(itr){

                     set.seed( 101 + itr);
                     YYstar = sample(YY, replace = TRUE);
                     Ystar  = (C.H%*%YYstar) + model$mu ; 

                     eM  <-  tryCatch({
                         estimate_shift_configuration(tr, Ystar, l1ou.options =model$l1ou.options);
                     }, error = function(e) {
                         print("l1OU error, return NA");
                         return(NA); }  );

                     if(all(is.na(eM))) {return(NA);}
                     return(eM$shift.configuration);
           }, mc.cores = nCores);

    valid.count <- 0;
    for( i in 1:length(shift.configuration.list)){
        if( all(is.na( shift.configuration.list[[i]] )) ){
            next;
        }
        valid.count <- valid.count + 1;
        detection.vec[ shift.configuration.list[[i]] ] = 
            detection.vec[ shift.configuration.list[[i]] ] + 1;
    }

    return(detection.vec/valid.count);
}

bootstrap_support_multivariate <- function(tr, model, nItrs, multicore=FALSE, nCores=2){

    Y = as.matrix(model$Y);
    stopifnot( length(model$alpha) == ncol(Y) );

    YY        = Y;
    C.Hlist   = list();
    for( idx in 1:ncol(Y) ){
        RE    = sqrt_OU_covariance(tr, alpha = model$alpha[[idx]] ); 
        C.IH  = t(RE$sqrtInvSigma); 
        C.Hlist[[idx]] = RE$sqrtSigma;
        YY[, idx]      = C.IH%*%(Y[, idx] - model$mu[ ,idx]);
    }

    detection.vec = rep(0, nrow(tr$edge));

    if( multicore == FALSE ){
        for(itr in 1:nItrs){

            Ystar   = YY;
            idx.vec = sample(1:nrow(YY), replace = TRUE);
            for( idx in 1:ncol(YY) ){
                YYstar        = YY[idx.vec, idx];
                Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx]; 
            }
            eM  <-  tryCatch({
                estimate_shift_configuration(tr, Ystar,  l1ou.options=model$l1ou.options);
            }, error = function(e) {
                print("l1OU error, return NA");
                return(NA); }  );

            if(all(is.na(eM))) {next;}
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1;
        }
    }

    shift.configuration.list = 
        mclapply(X=1:nItrs, FUN=function(itr){
                     Ystar   = YY;
                     set.seed( 101 + itr);
                     idx.vec = sample(1:nrow(YY), replace = TRUE);
                     for( idx in 1:ncol(YY) ){
                         YYstar        = YY[idx.vec, idx];
                         Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx]; 
                     }
                     eM  <-  tryCatch({
                         estimate_shift_configuration(tr, Ystar, l1ou.options = model$l1ou.options);
                     }, error = function(e) {
                         print("l1OU error, return NA");
                         return(NA); }  );

                     if(all(is.na(eM))) {return(NA);}
                     return(eM$shift.configuration);
                }, mc.cores = nCores);

    valid.count <- 0;
    for( i in 1:length(shift.configuration.list)){
        if( all(is.na( shift.configuration.list[[i]] )) ){
            next;
        }
        valid.count <- valid.count + 1;
        detection.vec[ shift.configuration.list[[i]] ] = 
            detection.vec[ shift.configuration.list[[i]] ] + 1;
    }

    return(detection.vec/valid.count);
}
