## generate the W matrix, feature vectors, as described in the doc.
generate_prediction_vec  <-  function(tr, 
                                      shift.configuration, 
                                      conv.regimes, 
                                      alpha, 
                                      ageMatrix=NULL, 
                                      designMatrix=F){

    ## In fact ageMatrix is the approximate design matrix.
    if(is.null(ageMatrix)){
        if ( is.na(alpha) ){
            X   <-  generate_design_matrix(tr, "apprX")
        }else{
            X   <-  generate_design_matrix(tr, "orgX",  alpha = alpha)
        }

        if(designMatrix){
            Cinvh <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot", normalize.tree.height=TRUE)$sqrtInvSigma )
            X     <- Cinvh%*%X
        }
        preds <- cbind(1, X[,shift.configuration])
        colnames(preds) <- c(0, shift.configuration)
        Z           <- generate_design_matrix(tr, "simpX")
        template.Z  <- cbind(1, Z[,shift.configuration])
        colnames(template.Z) <- c(0, shift.configuration) 
    }else{
        stopifnot(ncol(ageMatrix)==length(shift.configuration))
        stopifnot(alpha>0)

        preds <- cbind(1, 1-exp(-alpha*ageMatrix))
        colnames(preds) <- c(0, shift.configuration)

        template.Z  <- cbind(1, ageMatrix)
        colnames(template.Z) <- c(0, shift.configuration) 
    }

    ## now the coefficients in the linear regression represent the optimum values.
    ## rather than the shift values.

    for( i in 1:ncol(preds) ){
        for( j in 1:ncol(preds) ){
            if ( i == j )  
                next
            set1 <- which( template.Z[,i] > 0)
            set2 <- which( template.Z[,j] > 0)
            ## if edge j is an ancestor of i
            if ( all(set1 %in% set2) ) 
                preds[ ,j] <- preds[ ,j] - preds[ ,i]
        }
    } 

    W <- numeric()
    if ( length( conv.regimes ) > 0 ){
        for( i in 1:length(conv.regimes) ){
            set1 <- paste(conv.regimes[[i]])
            stopifnot( length(set1) > 0 )

            ## The equality constrain of two optimum values translates 
            ## to combining the corresponding predictors.
            if ( length(set1) > 1){
                W <- cbind( W, rowSums( preds[, paste(set1)] ) ) 
            } else{
                W <- cbind( W,  preds[, paste(set1)] ) 
            }
            colnames(W)[length(W[1,])] <- i
        }
    }

    return(W)
}

phylolm_interface_CR  <-  function(tr, Y, conv.regimes = list(), alpha=NA, fixed.alpha=FALSE, opt){

    shift.configuration <- opt$shift.configuration
    preds <- generate_prediction_vec(tr, shift.configuration, conv.regimes, alpha, ageMatrix=opt$ageMatrix)
    prev.val <- options()$warn
    options(warn = -1)

    if(fixed.alpha){
	    preds <- ifelse(preds>0,1,0)
	    fit <-  phylolm(Y~preds-1, phy  = tr, model = "OUfixedRoot",
			    starting.value = alpha,
			    lower.bound = alpha, 
			    upper.bound = alpha)

    }else{
	    fit <-  phylolm_CR(Y~preds-1, 
			       phy  = tr, 
			       model = "OUfixedRoot",
			       sc=shift.configuration,
			       cr=conv.regimes,
			       starting.value=alpha,
			       lower.bound=alpha/100
			       )
    }
    options(warn = prev.val)
    return(fit)
}

## compute the AICc score
cmp_AICc_CR  <-  function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    stopifnot( length(alpha) == ncol(Y) )

    nShifts    <- length( shift.configuration )
    nShiftVals <- length( conv.regimes ) -1## conv.regimes has intercept as an optimum value
    nTips      <- length( tree$tip.label )

    p <- nShifts + (nShiftVals + 3)*ncol(Y)
    N <- nTips*ncol(Y)
    df.1 <- 2*p + (2*p*(p+1))/(N-p-1) 
    if( p > N-2)  ##  for this criterion we should have p < N.
        return(Inf)
    df.2 <- 0
    score <- df.1
    for( i in 1:ncol(Y)){
        fit   <- phylolm_interface_CR(tree, matrix(Y[,i]), conv.regimes, alpha=alpha[[i]], opt=opt)
        if ( all( is.na( fit) ) ){ return(Inf) } 
        score <- score  -2*fit$logLik + df.2
    }
    return(score)
}

## compute the BIC score
cmp_BIC_CR <- function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    stopifnot( length(alpha) == ncol(Y) )

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nShiftVals <- length( conv.regimes ) - 1 
    nVariables <- ncol(Y)

    df.1  <- log(nTips)*(nShiftVals)
    score <- df.1
    #alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        df.2 <- log(nTips)*(nShifts + 3)
        fit  <- phylolm_interface_CR(tree, matrix(Y[,i]), conv.regimes, alpha=alpha[[i]], opt=opt)
        if ( all(is.na(fit)) ){ return(Inf) } 
        score <- score  -2*fit$logLik + df.2
    }
    return( score )
}



## compute the pBIC score
cmp_pBIC_CR  <-  function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1   <- 0
    df.1   <- 2*(nShifts)*log(nEdges-1)
    score  <- df.1
    #alpha  <- sigma2 <- logLik <- rep(0, ncol(Y))

    for(i in 1:ncol(Y)){
        fit   <- phylolm_interface_CR(tree, matrix(Y[,i]), conv.regimes, alpha=alpha[[i]], opt=opt)
        fit2  <- phylolm_interface_CR(tree, matrix(Y[,i]), conv.regimes, alpha=alpha[[i]], fixed.alpha=TRUE, opt=opt)
        if( all( is.na(fit) ) ){
           return(Inf)
        } 
        varY  <- var(Y[,i])
        ld    <- as.numeric(determinant(fit2$vcov * (fit$n - fit$d)/(varY*fit$n), log=T)$modulus)
        df.2  <- 2*log(nTips) - ld
        score <- score  -2*fit$logLik + df.2
    }
    return( score )
}


cmp_model_score_CR <- function(tree, Y, regimes=NULL, alpha=NA, opt){

    shift.configuration <- opt$shift.configuration

    if(is.null(regimes)){
        if(is.null(names(shift.configuration))){
            stop("the convergent regimes must be indicated through the names of the shift.configuration vector.")
        }
        cr.names <- names(shift.configuration)
        regimes <- list()
        idx <- 1
        for( cr in cr.names){
           regimes[[idx]] <- sort( shift.configuration[which(cr.names==cr)] )
           idx <- idx + 1
        }
    }

    if( opt$criterion == "AICc"){
        score <- cmp_AICc_CR(tree, Y, conv.regimes = regimes, alpha=alpha, opt=opt)
    } else if( opt$criterion == "pBIC"){
        score <- cmp_pBIC_CR(tree, Y, conv.regimes = regimes, alpha=alpha, opt=opt)
    } else if( opt$criterion == "BIC"){
        score <- cmp_BIC_CR(tree, Y, conv.regimes = regimes, alpha=alpha,  opt=opt)
    } else
        stop("undefined criterion for convergent evolution!")

    return(score)
}


generate_relation  <- function(tr, shift.configuration){

    s.p     = shift.configuration
    n.s.p   = length(s.p)
    nEdges  = Nedge(tr)

    M  <- numeric()
    tmp.s.p <- s.p
    for ( s1 in s.p ){
        tmp.s.p <- setdiff(tmp.s.p, s1)
        #for ( s2 in setdiff(s.p, s1) )
        for ( s2 in tmp.s.p )
        {
            nr <- rep(0, nEdges) 
            ## B M ADDed
            nr[[s1]] <-  1
            nr[[s2]] <- -1
            nr       <- nr[s.p]

            ## removing rows with only one +-1, we assume that previous step took care of redundancy.
            if ( length( which( abs(nr)>0) ) < 2)
                next

            M <- rbind(M, nr)
            rownames(M)[[length(rownames(M))]] <- paste0(s1," ",s2)
        }
    }

    colnames(M) <- s.p
    stopifnot( nrow(M) > 0 )
    return( M )
}

find_convergent_regimes  <-  function(tr, Y, alpha, criterion, regimes){
    stopifnot(ncol(Y)==1)
    library("magic")
    stopifnot(all( row.names(Y) == tr$tip.label))

    #alpha <- eModel$alpha
    Cinvh   <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot", normalize.tree.height=TRUE)$sqrtInvSigma )
    #Cinvh   <- t( cmp.OU.covariance(tr, alpha=alpha)$D ) 
    Y  <- Cinvh%*%Y

    shift.configuration <- unlist(regimes)
    X   <-  generate_prediction_vec(tr, shift.configuration, alpha=alpha, conv.regimes=regimes, designMatrix=TRUE)
    #X   <-  X[,-1]
    X   <- cbind(X,1)

    M   <- generate_relation(tr, 1:length(regimes))
    M   <- cbind(M,0)

    ###I multipled YY by a number to scale it up. If I don't genlasso doesn't return the whole solution path :S
    out     <- genlasso(100*Y, X, D=M, svd=T,  eps=0, approx=F, verbose=F)

    ### adding lambda=0 and removing the intercept
    out$beta   <- cbind(out$beta, coef(out, lambda=0)$beta )
    out$lambda <- c(out$lambda, 0)

    rownames(out$beta) <- colnames(X)

    ## removing the intercpet value from coefficients and adding it to a new vector.
    spots         <- length( out$beta[,1] )
    out$intercept <- out$beta[    spots, ]
    out$beta      <- out$beta[   -spots, ]
    M             <- M       [ , -spots  ]

    out$M <- M
    return(out)
}

## This method finds the convergent evolution model that maximizes the 
## criterion such as AICc. To find the optimum, it combines shifts into a 
## convergent regime or splits a regime if that decreases the criterion until 
## no progress.  NOTE: CR is a short for convergent regime.
estimate_convergent_regimes_surface  <-  function(model, opt){

    criterion <- opt$criterion
    Y         <- as.matrix(model$Y)
    tr        <- model$tree

    sc.prev <- sc  <- model$shift.configuration
    prev.min.score <- min.score <- Inf
    min.regimes    <- as.list(c(0,sc))
    ## elist represents the edgelist format of the regimes graph.
    ## At the beginning each regime forms a vertex with a self-loop. 
    elist.ref  <-  numeric()
    for(u in c(0,sc)){ elist.ref <- rbind( elist.ref, as.character(c(u,u)) ) }
    #current.num.cc <- length(sc)
    current.num.cc <- length(elist.ref)

    for( iter in 1:(2*length(sc)) ){

        has.progress <- FALSE
        ##NOTE: merge, add the edge (u,v), regimes if the IC decreases the most.

        run.list <- list()
	list.idx <- 1

        for( u in c(0, sc) ){ 
            for( v in c(0, sc) ){

                ##NOTE: test if we can add (u, v)
                if( u == v ){ next }
                ##FIXME: The convergence of the immediate shifts after the root
                ##to the background is pointless so don't check them.

                ## add the edge (u,v) to the graph. 
                elist <- as.matrix( rbind(elist.ref, as.character(c(u,v)) ) )
                g     <- graph_from_edgelist(elist, directed = FALSE)
                cc    <- decompose.graph(g) 
                ## check if it connects two connected components
                ## of the graph. If not, then it is redundant.
                if( length(cc) >= current.num.cc ){ next }
                ## extract the connected components 
                regimes <-  sapply(cc, function(x) as.numeric(names(V(x)))  )

                ## name each cr as the smallest shift index in it 
                sc.tmp  <-  sort(sc)
                for( thelist in lapply(regimes, sort) ){
                    names(sc.tmp)[sc.tmp %in% thelist]  <-  thelist[[1]] 
                }
                if(identical(names(sc.tmp), names(sc.prev))){ next }
                sc.prev <- sc.tmp

		if(!opt$parallel.computing){
			score   <-  cmp_model_score_CR(tr, Y, regimes, model$alpha, opt=opt)
			if( min.score > score ){
				min.score    <- score
				min.regimes  <- regimes
				elist.min    <- elist
				has.progress <- TRUE
			}
		}else{
			run.list[[list.idx]] <- list(elist=elist, regimes=regimes)
			list.idx <- list.idx+1
		}


            }
        }

	if( length(run.list)>0 && opt$parallel.computing ){

		RE.list <- mclapply( run.list, FUN=function(X){
					    return ( cmp_model_score_CR(tr, Y, X$regimes, model$alpha, opt=opt) )
			       }, mc.cores=opt$nCores)

		for(idx in 1:length(RE.list) ){
			IN <- run.list[[idx]]
			score   <- RE.list[[idx]] 

			if( min.score > score ){
				min.score    <- score
				min.regimes  <- IN$regimes
				elist.min    <- IN$elist
				has.progress <- TRUE
			}
		}
	}


        current.num.cc <- length(min.regimes) 
        elist.ref      <- elist.min

        ## break a CR if it increases the score
        ##FIXME: we don't need to check the very last edge we added
        if( has.progress){
            for( e.idx in 1:length(elist.ref[,1]) ){

                u <- elist.ref[e.idx, 1]
                v <- elist.ref[e.idx, 2]
                if(u==v){next}

                elist <- elist.ref[ -e.idx, ]
                g     <- graph_from_edgelist(elist, directed = FALSE)
                cc    <- decompose.graph(g) 

                if( length(cc) <= current.num.cc ){ next }

                regimes <- sapply(cc, function(x) as.numeric(names(V(x)))  )
                score   <- cmp_model_score_CR(tr, Y, regimes, model$alpha, opt=opt)

                if( min.score > score ){
                    min.score   <- score
                    elist.min   <- elist
                    min.regimes <- regimes
                }
            }
        }

        elist.ref      <- elist.min
        current.num.cc <- length(min.regimes) 

        #if no progress then terminate
        if( !has.progress ){ break }
        prev.min.score <- min.score
    }
    
    sc <- model$shift.configuration
    counter <- 1
    for (reg in min.regimes) {
        ## 0 represents the background(intercept) 
        background <- (0 %in% reg)
        for (item in sort(reg)) {
            names(sc)[which(sc == item)] <- ifelse(background, 0, counter)
        }
        counter <- counter + ifelse(background, 0, 1)
    }

    model$shift.configuration <- sc
    model$score               <- min.score 
    return(model)
}



#' Detects convergent regimes under an OU model
#'
#' Takes a model previously estimated by \code{\link{estimate_shift_configuration}},
#' including one or more traits and a configuration of evolutionary shifts, and detect which of these regime shifts
#' are convergent.
#'
#'@param model fitted object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@param method search method for finding convergent regimes. ``rr'' is based on genlasso,
#'  a regularized linear regression estimation. Currenly, this method can only accept a single trait.
#'  The default ``backward'' method is a heuristic similar to \code{surface_backward}
#'  in the \code{surface} package,
#'  using backward steps to repeatedly merge similar regimes into convergent regimes.
#'@param nCores number of processes to be created for parallel computing. If nCores=1 then it will run sequentially. Otherwise, it creates nCores processes by using mclapply function. For parallel computing it, requires parallel package.
#'
#'@examples
#' 
#'library(l1ou)
#'data("lizard.traits", "lizard.tree")
#'Y <- lizard.traits[, 1:1]
#' ## first fit a model to find individual shifts (no convergence assumed):
#'fit_ind <- estimate_shift_configuration(lizard.tree, Y, criterion="AICc")
#'fit_ind
#' ## then detect which of these shifts are convergent:
#'fit_conv <- estimate_convergent_regimes(fit_ind, criterion="AICc")
#'fit_conv
#'plot(fit_conv)
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
estimate_convergent_regimes  <-  function(model, 
                                        criterion=c("AICc", "pBIC", "BIC"),
                                        method=c("backward", "rr"),
					nCores=1
                                     ){
    opt <- list()
    opt$method <- match.arg(method)
    opt$criterion <- match.arg(criterion)
    opt$shift.configuration <- model$shift.configuration

    opt$nCores <- nCores
    opt$parallel.computing <- FALSE
    if( opt$nCores > 1){
        if(!require(parallel)){
            warning("parallel package is not available. The process will run sequentially.", immediate=TRUE)
            opt$nCores <- 1
        }else{
	    opt$parallel.computing <- TRUE
        }
    }
    opt$ageMatrix <- generate_design_matrix(model$tree, "apprX")[, model$shift.configuration]

    if(opt$method == "backward"){
        return(estimate_convergent_regimes_surface(model, opt=opt))
    }


    ## doing this for genlasso. FIXME: change the solver as it seems unstable.
    Y   <-  32*model$Y/lnorm(model$Y,l=2)
    Y   <-  as.matrix(Y)
    tr  <-  model$tree

    stopifnot( ncol(Y) == 1 ) # this method only works for univariate trait

    c.regimes <- prev.regimes <- all.regimes <- list()
    c.regimes[1:length(model$shift.configuration)] <- model$shift.configuration
    prev.min.score <- min.score <- Inf
    ar.counter <- 1
    
    ## similar to "backward" method. But here we may combine several shifts into convergent regimes
    ## at the same time therefore it is faster.
    for(iter in 1:length(model$shift.configuration) ){
        out  <-  find_convergent_regimes(tr, Y, model$alpha, opt$criterion, regimes = c.regimes )
        for(num.digits in c(12,13,15,16)){
            for( idx in 1:length(out$beta[1,]) ){
    
                new.est  <- round( out$M%*%out$beta[,idx], digits = num.digits)
                conv.reg <- rownames(out$M)[ which( new.est == 0) ]
    
                elist    <- numeric()
                for(e in conv.reg){## converting to numerical matrix.
                    elist <- rbind( elist, unlist(strsplit(e," ") ) )
                }
                for(e in paste(1:length(c.regimes)) ){
                    elist <- rbind( elist, c(e,e))
                }
    
                ## extracting connected components as the convergent regimes 
                g       <- graph.edgelist(elist, directed = FALSE)
                cc      <- decompose.graph(g, min.vertices = 0) 
                s.p.tmp <- model$shift.configuration
                regimes <- list()
                counter <- 1
    
                for(i in 1:length(cc) ){
    
                    regimes[[counter]] <- numeric()
                    the.cc <- as.numeric( setdiff( names(V(cc[[i]])), "r" ) )
                    for( vv in the.cc ){
                        regimes[[counter]] <- c(regimes[[counter]], c.regimes[[vv]])
                    }
                    s.p.tmp <- setdiff( s.p.tmp, regimes[[counter]] )
                    counter <- counter + 1
                }
    
                ## s.p.tmp contains shifts that are not in the graph
                for( s in s.p.tmp){
                    regimes[[counter]] <- s
                    counter <- counter + 1
                }
    
                if( identical(prev.regimes, regimes) )
                    next
                prev.regimes <- regimes
    
                res <- lapply( X=all.regimes, FUN=function(x){ return(identical(regimes,x)) }  )
                if( any( unlist(res)) )
                    next
    
                all.regimes[[ ar.counter ]] <- regimes
                ar.counter <- ar.counter + 1
    
                score <- cmp_model_score_CR(tr, Y, regimes, model$alpha, opt=opt)
    
                if( min.score > score ){
                    min.score      <- score
                    min.cr.regimes <- regimes
                    min.digits     <- num.digits
                }
            }
        }

        c.regimes <- min.cr.regimes
        if( min.score == prev.min.score )
            break;
        prev.min.score <- min.score

      }

      sc <- model$shift.configuration
      counter <- 1
      for( reg in c.regimes ){
          for( item in sort(reg) ){
              names(sc)[which(sc==item)] <- counter 
          }
          counter <- counter + 1
      }

      model$shift.configuration <- sc
      model$score  <-  cmp_model_score_CR(tr, model$Y, c.regimes, model$alpha, opt=opt)

      return(model)
}






 
