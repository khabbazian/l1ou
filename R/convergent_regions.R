generate_prediction_vec  <-  function(tr, 
                                      shift.configuration, 
                                      conv.regimes, 
                                      alpha, 
                                      ageMatrix=NULL, 
                                      designMatrix=F){

    if(is.null(ageMatrix)){
        if ( is.na(alpha) ){
            X   <-  generate_design_matrix(tr, "apprX")
        }else{
            X   <-  generate_design_matrix(tr, "orgX",  alpha = alpha)
        }

        if(designMatrix){
            Cinvh <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
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


    for( i in 1:ncol(preds) ){
        for( j in 1:ncol(preds) ){
            if ( i == j )  
                next
            set1 <- which( template.Z[,i] > 0)
            set2 <- which( template.Z[,j] > 0)
            if ( all(set1 %in% set2) ) 
                preds[ ,j] <- preds[ ,j] - preds[ ,i]
        }
    } 

    preds.2 <- numeric()
    if ( length( conv.regimes ) > 0 ){
        for( i in 1:length(conv.regimes) ){
            set1 <- paste(conv.regimes[[i]])
            stopifnot( length(set1) > 0 )

            if ( length(set1) > 1){
                preds.2 <- cbind( preds.2, rowSums( preds[, paste(set1)] ) ) 
            } else{
                preds.2 <- cbind( preds.2,  preds[, paste(set1)] ) 
            }
            colnames(preds.2)[length(preds.2[1,])] <- i
        }
    }

    preds <- preds.2
    return(preds)
}

phylolm_interface_CR  <-  function(tr, Y, shift.configuration, conv.regimes = list(), alpha=NA){

    preds <- generate_prediction_vec(tr, shift.configuration, conv.regimes, alpha)
    prev.val <- options()$warn
    options(warn = -1)
    #fit <-  phylolm(Y~preds-1, phy  = tr, model = "OUfixedRoot",
    #                  starting.value = alpha,
    #                  lower.bound = alpha, 
    #                  upper.bound = alpha)

    fit <-  phylolm_CR(Y~preds-1, 
                       phy  = tr, 
                       model = "OUfixedRoot",
                       sc=shift.configuration,
                       cr=conv.regimes,
                       starting.value=alpha,
                       lower.bound=alpha/100
                       )
    options(warn = prev.val)
    return(fit)
}


cmp_AICc_CR  <-  function(tr, Y, shift.configuration, conv.regimes, alpha){

    stopifnot( length(alpha) == ncol(Y) )

    nShifts    <- length( shift.configuration )
    ## -1, cause "conv.regimes" contains the intercept
    nShiftVals <- length( conv.regimes ) - 1 
    nTips      <- length( tr$tip.label )

    p <- nShifts + (nShiftVals + 3)*ncol(Y)
    N <- nTips*ncol(Y)
    df.1 <- 2*p + (2*p*(p+1))/(N-p-1) 
    if( p > N-2)  ##  for this criterion we should have p < N.
        return(Inf)
    df.2 <- 0
    score <- df.1
    for( i in 1:ncol(Y)){
        fit   <- phylolm_interface_CR(tr, matrix(Y[,i]), shift.configuration, conv.regimes, alpha=alpha[[i]])
        if ( all( is.na( fit) ) ){
            return(Inf)
        } 
        score <- score  -2*fit$logLik + df.2
    }
    return(score)
}


cmp_pBIC_CR  <-  function(tr, Y, shift.configuration, conv.regimes, alpha){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tr)
    nTips   = length(tr$tip.label)

    df.1   <- 0
    df.1    <- (nShifts)*log(nEdges-1)
    for( i in 1:nShifts){
        df.1 <- df.1 + log(nEdges-1-i)
        df.1 <- df.1 - log(i)
    }
    df.1 <- df.1 - sum( log( factorial(unlist(lapply(conv.regimes, length))) ) )
    df.1 <- 2*df.1

    score   <- df.1
    for(i in 1:ncol(Y)){
        fit   <- phylolm_interface_CR(tr, matrix(Y[,i]), shift.configuration, conv.regimes, alpha=alpha[[i]])
        if( all( is.na(fit) ) ){
           return(Inf)
        } 
        varY <- var(Y[,i])
        ld    <- as.numeric(determinant(fit$vcov * (fit$n - fit$d)/(varY*fit$n), log=T)$modulus)
        df.2  <- 3*log(nTips) - ld
        score <- score  -2*fit$logLik + df.2
    }
    return( score )
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
    Cinvh   <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
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


cmp_model_score_CR  <-  function(tr, Y, sc, regimes, criterion, alpha){
     if( criterion == "AICc"){
         score <- cmp_AICc_CR(tr, Y, sc, conv.regimes = regimes, alpha=alpha)
     } else { #if( criterion == "pBIC")
         score <- cmp_pBIC_CR(tr, Y, sc, conv.regimes = regimes, alpha=alpha)
     }
     return(score)
}

estimate_convergent_regimes_surface  <-  function(model, 
                                        criterion = c("AICc", "pBIC")
                                        ){
    criterion <- match.arg(criterion)
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
        ##NOTE: this part can be parallelized.
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

                ## this part has the computational overhead.
                score   <-  cmp_model_score_CR(tr, Y, model$shift.configuration, regimes, criterion, model$alpha)

                if( min.score > score ){
                    min.score    <- score
                    min.regimes  <- regimes
                    elist.min    <- elist
                    has.progress <- TRUE
                }

            }
        }

        current.num.cc <- length(min.regimes) 
        elist.ref      <- elist.min

        ## break a cr if it reduces the score
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
                score   <- cmp_model_score_CR(tr, Y, model$shift.configuration, regimes, criterion, model$alpha)

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
#' Given estimated evolutionary shift positions this function finds convergent regimes.
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@param method method for finding convergent regimes. ``rr'' is based on genlasso, a regularized linear regression estimation. ``backward'' is a heuristic method similar to \code{surface_backward} function.
#'@examples
#' 
#'library(l1ou)
#'data("lizard.traits", "lizard.tree")
#'Y <- lizard.traits[, 1:1]
#'
#'eModel <- estimate_shift_configuration(lizard.tree, Y, criterion="AICc")
#'eModel <- estimate_convergent_regimes(eModel, criterion="AICc")
#'eModel
#'plot(eModel)
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
estimate_convergent_regimes  <-  function(model, 
                                        criterion=c("AICc", "pBIC"),
                                        method=c("backward", "rr")
                                     ){
    method  <-   match.arg(method)
    if(method == "backward"){
        return(estimate_convergent_regimes_surface(model, criterion))
    }

    criterion   <-   match.arg(criterion)
    Y   <-  32*model$Y/norm(model$Y)
    Y   <-  as.matrix(Y)
    tr  <-  model$tree

    stopifnot( ncol(Y) == 1 ) # this method only works for univariate trait

    c.regimes <- prev.regimes <- all.regimes <- list()
    c.regimes[1:length(model$shift.configuration)] <- model$shift.configuration
    prev.min.score <- min.score <- Inf
    ar.counter <- 1
    
    for(iter in 1:length(model$shift.configuration) ){
        out  <-  find_convergent_regimes(tr, Y, model$alpha, criterion, regimes = c.regimes )
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
    
                score <- cmp_model_score_CR(tr, Y, model$shift.configuration, regimes, criterion, model$alpha)
    
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
      model$score  <-  cmp_model_score_CR(tr, model$Y, model$shift.configuration, c.regimes, criterion, model$alpha)

      return(model)
}






 
