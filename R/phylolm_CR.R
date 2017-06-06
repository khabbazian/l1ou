phylolm_CR <- function(formula, data=list(), phy, 
                         model=c("OUrandomRoot","OUfixedRoot"),
                         lower.bound=NULL, upper.bound=NULL, starting.value=NULL,
                         sc, cr, ...) 
{
    ## initialize	
    if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
    model = match.arg(model)	
    mytree=phy

    if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
    if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	

    tol = 1e-10	
    phy = reorder(phy,"pruningwise")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    externalEdge <- des <=n

    mf = model.frame(formula=formula,data=data)
    if (nrow(mf)!=length(phy$tip.label))
        stop("number of rows in the data does not match the number of tips in the tree.")
    if (is.null(rownames(mf))) {
        warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
        data.names = phy$tip.label 
    }
    else data.names = rownames(mf)
    order = match(data.names, phy$tip.label)
    if (sum(is.na(order))>0) {
        warning("data names do not match with the tip labels.\n")
        rownames(mf) = data.names
    } else {
        tmp = mf
        rownames(mf) = phy$tip.label
        mf[order,] = tmp[1:nrow(tmp),]
    }
    X = model.matrix(attr(mf, "terms"), data=mf)
    y = model.response(mf)
    d = ncol(X)
    OU = c("OUrandomRoot","OUfixedRoot")
    flag = 0 # flag and D are used for OU model if tree is not ultrametric:
    D = NULL #            for the generalized 3-point structure

    ## preparing for OU model
    D = numeric(n)
    if (!is.ultrametric(phy)) stop("the tree has to be ultrametric")
    dis = pruningwise.distFromRoot(phy)
    Tmax = max(dis[1:n])
    D = Tmax - dis[1:n]
    D = D - mean(D)
    phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
    ## phy is now ultrametric, with height Tmax:


    ## calculate Tmax = average distance from root to tips,
    ## to choose appropriate starting values later # fixit
    dis = pruningwise.distFromRoot(phy)[1:n]
    Tmax = mean(dis)

    ## Default bounds
    bounds.default = matrix(c(1e-7/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0), ncol=2, byrow=TRUE)
    rownames(bounds.default) = c("alpha","lambda","kappa","delta","rate")
    colnames(bounds.default) = c("min","max")

    ## Default starting values
    starting.values.default = c(0.5/Tmax,0.5,0.5,0.5,-1/Tmax) 
    names(starting.values.default) = c("alpha","lambda","kappa","delta","rate")

    ## User defined bounds and starting values
    if (is.null(lower.bound)) {
        lower.bound = bounds.default[1,1]
    }
    if (is.null(upper.bound)) {
        upper.bound = bounds.default[1,2]
    }
    if (is.null(starting.value)) {
        starting.value = starting.values.default[1]
    }	

    ## preparing for general use of "parameter" for branch length transformation
    prm = list(myname = starting.value)
    names(prm) = "alpha"

    ## log-likelihood, computation using the three-point structure
    ole= 4 + 2*d + d*d # output length
    loglik <- function(parameters,y,X) {

        tree = transf.branch.lengths(phy,model,parameters=parameters,
                                     check.pruningwise=F,check.ultrametric=F,D=D,check.names=F)$tree

        tmp <- .C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
               as.integer(1),as.integer(d),as.integer(ROOT),as.double(tree$root.edge),as.double(tree$edge.length),
               as.integer(des), as.integer(anc), as.double(as.vector(y)), as.double(as.vector(X)),
               result=double(ole))$result # tmp has, in this order:

        ## logdetV, 1'V^{-1}1, y'V^{-1}1, y'V^{-1}y, X'V^{-1}1, X'V^{-1}X, X'V^{-1}y
        comp = list(vec11=tmp[2], y1=tmp[3], yy=tmp[4], X1=tmp[5:(4+d)],
                    XX=matrix(tmp[(5+d):(ole-d)], d,d),Xy=tmp[(ole-d+1):ole],logd=tmp[1])

	if(rcond(comp$XX) < .Machine$double.eps)
        	return(list(n2llh=-Inf, betahat = NA, sigma2hat=NA, vcov=NA))

        invXX <- solve(comp$XX, tol=.Machine$double.eps) 
        betahat <- invXX%*%comp$Xy
        sigma2hat <- as.numeric((comp$yy - 2*t(betahat)%*%comp$Xy + t(betahat)%*%comp$XX%*%betahat)/n)
        if (sigma2hat<0) { stop("sigma2hat<0") }
        n2llh = as.numeric( n*log(2*pi) + n + n*log(sigma2hat) + comp$logd) # -2 log-likelihood

        ## because diag matrix used for generalized 3-point structure is exp(alpha diag(D))
        vcov = sigma2hat*invXX*n/(n-d)
        return(list(n2llh=n2llh, betahat = as.vector(betahat), sigma2hat=sigma2hat,vcov=vcov))
    }

    ## Fitting
    lower = lower.bound
    upper = upper.bound
    start = starting.value

    {
        ## Optimization of phylogenetic correlation parameter is needed
        if ((lower>start)||(upper<start))
            stop("The starting value is not within the bounds of the parameter.")

        ageMatrix <- generate_design_matrix(mytree, type="apprX")[,sc]
        minus2llh <- function(logvalue) {
            prm[[1]]=exp(logvalue)
            preds <- generate_prediction_vec(mytree, sc,
                                                    cr,  alpha=prm[[1]],
                                                     ageMatrix=ageMatrix)
            loglik(prm, y, X=preds)$n2llh
            #loglik(prm, y, X=X)$n2llh
        }

        if (lower==upper) {
            MLEvalue = lower
        } else {
            logstart = log(start)
            loglower = log(lower)
            logupper = log(upper)

            opt <- optim(logstart, fn = minus2llh,
                         method = "L-BFGS-B", lower=loglower, upper = logupper, ...)
            MLEvalue = as.numeric(exp(opt$par))
        }
        if ((isTRUE(all.equal(MLEvalue,lower, tol=tol)))||(isTRUE(all.equal(MLEvalue,upper,tol=tol)))) {
            matchbound = 1
            if ((model %in% c("lambda","kappa"))&&(MLEvalue == 1)) matchbound=0
            if ((model == "EB")&&(MLEvalue == 0)) matchbound=0
            if (matchbound)
                warning(paste("the estimation of", names(prm), 'matches the upper/lower bound for this parameter.
                              You may change the bounds using options "upper.bound" and "lower.bound".\n'))
        }
        prm[[1]] = MLEvalue
        BMest = loglik(prm, y, X)
        if (model %in% OU)
            BMest$sigma2hat = 2*MLEvalue * BMest$sigma2hat # was "gamma" originally: sigma2 = 2 alpha gamma
        results <- list(coefficients=BMest$betahat, sigma2=BMest$sigma2hat, optpar=MLEvalue,
                        logLik=-BMest$n2llh/2, p=2+d, aic=2*(2+d)+BMest$n2llh, vcov = BMest$vcov)
    }

    names(results$coefficients) = colnames(X)
    colnames(results$vcov) = colnames(X)
    rownames(results$vcov) = colnames(X)
    results$fitted.values = X %*% results$coefficients
    results$residuals = y - results$fitted.values
    results$mean.tip.height = Tmax
    results$y = y
    results$X = X
    results$n = n
    results$d = d
    results$formula = formula
    results$call = match.call()
    results$model = model
    class(results) = "phylolm"
    return(results)
}

################################################
################################################

print.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("Call:\n")
    print(x$call)
    cat("\n")
    aiclogLik = c(x$aic,x$logLik)
    names(aiclogLik) = c("AIC","logLik")
    print(aiclogLik, digits = digits)
    cat("\nParameter estimate(s) using ML:\n")
    if (!is.null(x$optpar)) {
        if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
        if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
        if (x$model=="EB") cat("rate:",x$optpar)
        cat("\n")
    }
    cat("sigma2:",x$sigma2,"\n")
    cat("\nCoefficients:\n")
    print(x$coefficients)
}
################################################
summary.phylolm <- function(object, ...) {
    se <- sqrt(diag(object$vcov))
    tval <- coef(object) / se
    TAB <- cbind(Estimate = coef(object), StdErr = se, t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$n - object$d))
    res <- list(call=object$call, coefficients=TAB,
                residuals = object$residuals, sigma2 = object$sigma2,
                optpar=object$optpar, logLik=object$logLik,
                df=object$p, aic=object$aic, model=object$model,
                mean.tip.height=object$mean.tip.height)
    class(res) = "summary.phylolm"
    res
}
################################################
print.summary.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    aiclogLik = c(x$aic,x$logLik)
    names(aiclogLik) = c("AIC","logLik")
    print(aiclogLik, digits = digits)
    r <- zapsmall(quantile(x$residuals), digits + 1)
    names(r) <- c("Min", "1Q", "Median", "3Q", "Max")
    cat("\nRaw residuals:\n")
    print(r, digits = digits)

    cat("\nMean tip height:",x$mean.tip.height)
    cat("\nParameter estimate(s) using ML:\n")
    if (!is.null(x$optpar)) {
        if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
        if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
        if (x$model=="EB") cat("rate:",x$optpar)
        cat("\n")
    }
    cat("sigma2:",x$sigma2,"\n")
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
    if (!is.null(x$optpar)) {
        cat("\nNote: p-values are conditional on ")
        if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha=",x$optpar,".",sep="")
        if (x$model %in% c("lambda","kappa","delta")) cat(x$model,"=",x$optpar,".",sep="")
        if (x$model=="EB") cat("rate=",x$optpar,".",sep="")
    }
    cat("\n")
}
################################################
residuals.phylolm <-function(object,type=c("response"), ...){
    type <- match.arg(type)
    object$residuals	 
}
################################################
vcov.phylolm <- function(object, ...){
    object$vcov
}
################################################
logLik.phylolm <- function(object, ...){
    res = list(logLik = object$logLik, df = object$p)
    class(res) = "logLik.phylolm"
    res
}
print.logLik.phylolm <- function (x, ...) {
    cat("'log Lik.' ",x$logLik," (df=",x$df,")\n", sep = "")
}
AIC.logLik.phylolm <- function(object, k=2, ...) {
    return(k*object$df - 2*object$logLik)
}
AIC.phylolm <- function(object, k=2, ...) {
    return(AIC(logLik(object),k))
}
extractAIC.phylolm <- function(fit, scale, k=2, ...) {
    c(fit$p, - 2*fit$logLik + k * fit$p)
}
nobs.phylolm <- function(object, ...){
    return(object$n)
}
################################################
predict.phylolm <- function(object, newdata=NULL, ...){
    if (object$model=="trend")
        stop("Predicting for trend model has not been implemented.")
    if(is.null(newdata)) y <- fitted(object)
    else{			
        X = model.matrix(delete.response(terms(formula(object))),data = newdata)
        y <- X %*% coef(object)
    }
    y
}
################################################
plot.phylolm <-function(x, ...){
    plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value", ...)
}
################################################
