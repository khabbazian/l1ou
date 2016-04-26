

generate_prediction_vec <- function(tr, shift.configuration, conv.regimes, alpha, designMatrix=F){

    nTips           = length(tr$tip.label)
    if ( is.na(alpha) ){
        #FIXME in package
        X           = l1ou:::generate_design_matrix(tr, "apprX")
    }else{
        #FIXME in package
        X           = l1ou:::generate_design_matrix(tr, "orgX",  alpha = alpha)
    }

    if(designMatrix){
        Cinvh   = t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
        #Cinvh   = t( cmp.OU.covariance(tr, alpha=alpha)$D ) 
        X  = Cinvh%*%X
    }

    Z               = l1ou:::generate_design_matrix(tr, "simpX")
    preds           = X[,shift.configuration] 
    template.Z      = Z[,shift.configuration] 

    colnames(preds)      = shift.configuration 
    colnames(template.Z) = shift.configuration 

    for( i in 1:ncol(preds) ){
        for( j in 1:ncol(preds) ){
            if ( i == j )  
                next
            set1 = which( template.Z[,i] > 0)
            set2 = which( template.Z[,j] > 0)
            if ( all(set1 %in% set2) ) 
                preds[ ,j] = preds[ ,j] - preds[ ,i]
        }
    } 

    preds.2 = numeric()
    if ( length( conv.regimes ) > 0 ){
        for( i in 1:length(conv.regimes) ){
            set1 = paste(conv.regimes[[i]])
            stopifnot( length(set1) > 0 )

            if ( length(set1) > 1){
                preds.2 = cbind( preds.2, rowSums( preds[,set1] ) ) 
            } else{
                preds.2 = cbind( preds.2,  preds[,set1] ) 
            }
            colnames(preds.2)[length(preds.2[1,])] = i
        }
    }

    #preds = cbind(1, preds)
    preds = cbind(1, preds.2)
    return(preds)

}

my.phylolm.interface.new <- function(tr, Y, shift.configuration, conv.regimes = list(), alpha=NA){

    preds = generate_prediction_vec(tr, shift.configuration, conv.regimes, alpha)
    options(warn = -1)
    fit    <-   phylolm(Y~preds-1, phy  = tr, model = "OUfixedRoot",
                        starting.value = alpha,
                        lower.bound = alpha, 
                        upper.bound = alpha)
    options(warn = 0)
    return(fit)
}



cmp.AICc.new <- function(tr, Y, shift.configuration, conv.regimes, alpha){

    stopifnot( length(alpha) == ncol(Y) )

    nShifts    = length( shift.configuration )
    nShiftVals = length( conv.regimes ) 
    nTips      = length( tr$tip.label )

    p = nShifts + (nShiftVals + 3)*ncol(Y)
    N = nTips*ncol(Y)
    df.1 = 2*p + (2*p*(p+1))/(N-p-1) 
    if( p > N-2)  ##  for this criterion we should have p < N.
        return(Inf)
    df.2 = 0
    score = df.1
    for( i in 1:ncol(Y)){
        fit   = my.phylolm.interface.new(tr, matrix(Y[,i]), shift.configuration, conv.regimes, alpha=alpha[[i]])
        if ( all( is.na( fit) ) ){
            return(Inf)
        } 
        score = score  -2*fit$logLik + df.2
    }
    return(score)
}


cmp.pBIC.new <- function(tr, Y, shift.configuration, conv.regimes, alpha){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tr)
    nTips   = length(tr$tip.label)

    df.1   = 0
    df.1    = (nShifts)*log(nEdges-1)
    for( i in 1:nShifts){
        df.1 = df.1 + log(nEdges-1-i)
        df.1 = df.1 - log(i)
    }
    df.1 = df.1 - sum( log( factorial(unlist(lapply(conv.regimes, length))) ) )
    df.1 = 2*df.1

    score   = df.1
    for(i in 1:ncol(Y)){
        fit   = my.phylolm.interface.new(tr, matrix(Y[,i]), shift.configuration, conv.regimes, alpha=alpha[[i]])
        if( all( is.na(fit) ) ){
           return(Inf)
        } 
        ld    = as.numeric(determinant(fit$vcov * (fit$n - fit$d)/fit$n, log=T)$modulus)
        df.2  = 3*log(nTips) - ld
        score = score  -2*fit$logLik + df.2
    }
    return( score )
}


generate_relation <-function(tr, shift.configuration){

    s.p     = shift.configuration
    n.s.p   = length(s.p)
    nEdges  = Nedge(tr)

    M  = numeric()
    tmp.s.p = s.p
    for ( s1 in s.p ){
        tmp.s.p = setdiff(tmp.s.p, s1)
        #for ( s2 in setdiff(s.p, s1) )
        for ( s2 in tmp.s.p )
        {
            nr = rep(0, nEdges) 
            ## B M ADDed
            nr[[s1]] =  1
            nr[[s2]] = -1
            nr       = nr[s.p]

            ## removing rows with only one +-1, we assume that previous step took care of redundancy.
            if ( length( which( abs(nr)>0) ) < 2)
                next

            M = rbind(M, nr)
            rownames(M)[[length(rownames(M))]] = paste0(s1," ",s2)
        }
    }

    colnames(M) = s.p
    stopifnot( nrow(M) > 0 )
    return( M )
}

find_convergent_regions <- function(tr, Y, alpha, criterion, regimes){
    stopifnot(ncol(Y)==1)
    library("magic")
    stopifnot(all( row.names(Y) == tr$tip.label))

    #alpha = eModel$alpha
    Cinvh   = t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
    #Cinvh   = t( cmp.OU.covariance(tr, alpha=alpha)$D ) 
    Y  = Cinvh%*%Y

    shift.configuration = unlist(regimes)
    X   =  generate_prediction_vec(tr, shift.configuration, alpha=alpha, conv.regimes=regimes, designMatrix=TRUE)
    X   =  X[,-1]
    X   = cbind(X,1)

    M   = generate_relation(tr, 1:length(regimes))
    M   = cbind(M,0)

    ###I multipled YY by a number to scale it up. If I don't genlasso doesn't return the whole solution path :S
    out     = genlasso(100*Y, X, D=M, svd=T,  eps=0, approx=F, verbose=F)

    ### adding lambda=0 and removing the intercept
    out$beta   = cbind(out$beta, coef(out, lambda=0)$beta )
    out$lambda = c(out$lambda, 0)

    rownames(out$beta) = colnames(X)

    ## removing the intercpet value from coefficients and adding it to a new vector.
    spots         = length( out$beta[,1] )
    out$intercept = out$beta[    spots, ]
    out$beta      = out$beta[   -spots, ]
    M             = M       [ , -spots  ]

    out$M = M
    return(out)
}


#' Detects convergent regions under an OU model
#'
#' Given estimated evolutionary shifts this function finds convergent regions
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@examples
#' 
#'library(l1ou)
#'data("lizard.traits", "lizard.tree")
#'Y <- lizard.traits[, 1:1]
#'
#'eModel <- estimate_shift_configuration(lizard.tree, Y, criterion = "AICc")
#'eModel <- estimate_convergent_regions(eModel, criterion="AICc")
#'plot(eModel)
#'
#'@export
estimate_convergent_regions <- function(model, 
                                        criterion = c("AICc", "pBIC")
                                        ){

    criterion  <-  match.arg(criterion)
    Y  <- model$Y
    Y  <- as.matrix(Y)
    tr <- model$tree

    stopifnot( ncol(Y) == 1 ) # this method only works for univariate trait

    c.regimes = prev.regimes = all.regimes = list()
    c.regimes[1:length(model$shift.configuration)] = model$shift.configuration
    prev.min.score = min.score    = Inf
    ar.counter   = 1
    
    for(iter in 1:length(model$shift.configuration) ){
        out <- find_convergent_regions(tr, Y, model$alpha, criterion, regimes = c.regimes )
        for(num.digits in c(12,13,15,16)){
            for( idx in 1:length(out$beta[1,]) ){
    
                new.est  = round( out$M%*%out$beta[,idx], digits = num.digits)
                conv.reg = rownames(out$M)[ which( new.est == 0) ]
    
                elist    = numeric()
                for(e in conv.reg){## converting to numerical matrix.
                    elist = rbind( elist, unlist(strsplit(e," ") ) )
                }
                for(e in paste(1:length(c.regimes)) ){
                    elist = rbind( elist, c(e,e))
                }
    
                g     = graph.edgelist(elist, directed = FALSE)
                ## extracting connected components as the convergent regions 
                cc    = decompose.graph(g, min.vertices = 0)
    
                s.p.tmp = model$shift.configuration
                regimes = list()
                counter = 1
    
                for(i in 1:length(cc) ){
    
                    regimes[[counter]] = numeric()
                    the.cc = as.numeric( setdiff( names(V(cc[[i]])), "r" ) )
                    for( vv in the.cc ){
                        regimes[[counter]] = c(regimes[[counter]], c.regimes[[vv]])
                    }
                    s.p.tmp = setdiff( s.p.tmp, regimes[[counter]] )
                    counter = counter + 1
                }
    
                ## s.p.tmp contains shifts that are not in the graph
                for( s in s.p.tmp){
                    regimes[[counter]] = s
                    counter = counter + 1
                }
    
                if( identical(prev.regimes, regimes) )
                    next
                prev.regimes = regimes
    
                res = lapply( X=all.regimes, FUN=function(x){ return(identical(regimes,x)) }  )
                if( any( unlist(res)) )
                    next
    
                all.regimes[[ ar.counter ]] = regimes
                ar.counter <- ar.counter + 1
    
    
                if( criterion == "AICc"){
                    score = cmp.AICc.new(tr, Y, model$shift.configuration, conv.regimes = regimes, alpha=model$alpha)
                }
                if( criterion == "pBIC"){
                    score = cmp.pBIC.new(tr, Y, model$shift.configuration, conv.regimes = regimes, alpha=model$alpha)
                }
    
                if( min.score > score ){
                    min.score      = score
                    min.cr.regimes = regimes
                    min.digits     = num.digits
                }
    
            }
        }

        c.regimes = min.cr.regimes
        if( min.score == prev.min.score )
            break;
        prev.min.score = min.score

      }

      sc = model$shift.configuration
      for( len in 1:length(c.regimes) ){
          for( item in c.regimes[[len]] ){
              names(sc)[which(sc==item)] = len
          }
      }

      model$shift.configuration = sc
      model$score = min.score 
      return(model)
}



 
