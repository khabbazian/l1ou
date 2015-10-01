# 
#'Detecting Evolutionary Shifts
#'
#'@param tr input x axis 
#'@param Y input y axis 
#'@param max.nShifts  maximum number of shifts by default it is half of the number of tips
#'@param criterion  which can be c("pBIC", "pBICess", "mBIC", "BIC", "AIC", "AICc")
#'@param root.model which can be c("OUrandomRoot", "OUfixedRoot")
#'@param silence  a flag for writing to output
#'@param alpha.upper by default is .alpha.upper.bound(tr)
#'@param alpha.lower by default is 0
#'@param standardize if TRUE then in multivariate case the inpute traits will be normalized
#'@param num.top.placements  by default it is max.nShifts/2
#'@param edge.length.threshold the edge length that should be considerd zero by default 10*.Machine$double.eps
#'@param grp.delta  input parameter to figure out the lambda sequence for grplasso by defualt it is 1/16
#'@param grp.seq.ub inpute parameter to figure out the lambda sequence for grplasso by defualt it is 5
#'
#'@return eModel estimated model
#'
#'@examples
#' 
#' library("l1ou"); 
#' library("ape");
#'
#' trFileName = "GA_Anolis_MCC.tre"
#' 
#' tr <- read.tree(paste0("data/lizards/", trFileName) );
#' tr <- reorder(tr, "postorder");
#' tr <- normalize.tree(tr);
#' 
#' responseMatrix <- read.csv("data/lizards/table_s2.csv");
#' responseMatrix <- responseMatrix[order(responseMatrix[,1]),  ]; 
#' responseMatrix <- responseMatrix[order(order(tr$tip.label)), ];
#' 
#' Y           <- as.matrix( responseMatrix[, 2] );
#' rownames(Y) <- responseMatrix[, 1];
#' 
#' eModel <- est.shift.placement(tr, Y);
#'
#' print(eModel$shift.placement);
#'
#'@export

est.shift.placement <- function(tr, Y, 
           max.nShifts           = floor(length(tr$tip.label)/2), 
           criterion             = c("pBIC", "pBICess", "mBIC", "BIC", "AIC", "AICc"), 
           root.model            = c("OUrandomRoot", "OUfixedRoot"),
           silence               = TRUE,
           alpha.upper           = .alpha.upper.bound(tr), 
           alpha.lower           = 0,
           standardize           = TRUE,
           num.top.placements    = max.nShifts/2,
           edge.length.threshold = 10*.Machine$double.eps,
           grp.delta             = 1/16,
           grp.seq.ub            = 5
     ){

## loading required package
    library("ape");
    library("phylolm");
    library("igraph");

    source("linear_alg_OU_covariance.R");
    source("tools.R");

    ## unifying the types
    Y  <-  as.matrix(Y);

    ## setting up options
    l1ou.options  <-  list();

    l1ou.options$silence               <- silence;
    l1ou.options$use.saved.scores      <- TRUE;

    l1ou.options$max.nShifts           <- max.nShifts;
    l1ou.options$criterion             <- match.arg(criterion);
    l1ou.options$root.model            <- match.arg(root.model);

    ##TODO: return warning if est alpha is close to its upperbound
    l1ou.options$alpha.upper.bound     <- alpha.upper;
    l1ou.options$alpha.lower.bound     <- alpha.lower;

    l1ou.options$edge.length.threshold <- edge.length.threshold;
    l1ou.options$num.top.placements    <- num.top.placements;

    l1ou.options$standardize           <- standardize;
    l1ou.options$grp.seq.ub            <- grp.seq.ub;
    l1ou.options$grp.delta             <- grp.delta;
    l1ou.options$Z                     <- generate.design.matrix(tr, "simpX");

    ##NOTE: checking the assumptions 
    stopifnot(is.ultrametric(tr));
    #stopifnot(is.binary.tree(tr));
    stopifnot(identical(tr$edge , reorder(tr, "postorder")$edge));
    stopifnot(nrow(Y) == length(tr$tip.label));
    stopifnot(all( row.names(Y) == tr$tip.label));

    if(l1ou.options$use.saved.scores){
        library("Rcpp");
        sourceCpp("placement_score_db.cpp");
    }

    if( ncol(Y) == 1 ){  #univariate l1ou
        eModel1 = .est.shift.placement.known.alpha(tr, Y, est.alpha=TRUE,      opt=l1ou.options);
        eModel  = .est.shift.placement.known.alpha(tr, Y, alpha=eModel1$alpha, opt=l1ou.options);
        ## in case the first step finds a better solution
        if( eModel$score>eModel1$score )
            eModel = eModel1;

    } 
    if( ncol(Y) > 1 ){ #multivariate l1ou
        eModel1 = .est.shift.placement.known.alpha.multivariate(tr, Y, est.alpha=TRUE,      opt=l1ou.options);
        eModel  = .est.shift.placement.known.alpha.multivariate(tr, Y, alpha=eModel1$alpha, opt=l1ou.options);
        ## in case the first step finds a better solution
        if( eModel$score>eModel1$score )
            eModel = eModel1;
    } 

    ## I don't like this way of coding but I have to clear the static db object in the c++ code manually.:S
    ## TODO: can we implement it differently?

    if( l1ou.options$use.saved.scores){
        erase_placement_score_db();
    }

    return(eModel);
}

.est.shift.placement.known.alpha <- function(tr, Y, alpha=0, est.alpha=FALSE, opt){  

    stopifnot( alpha >=0 );

    library("lars");
    if ( est.alpha ){ ## BM model
        X   = generate.design.matrix(tr, "apprX");
        Cinvh   = t( cmp.OU.covariance(tr, alpha=0)$D ); 
    } else{           ## OU model
        X   = generate.design.matrix(tr, "orgX", alpha=alpha );
        Cinvh   = t( cmp.OU.covariance(tr, alpha=alpha)$D ); 
    }

    to.be.removed   = c(length(tr$edge.length), which(tr$edge.length < opt$edge.length.threshold));

    YY  = Cinvh%*%Y;
    XX  = Cinvh%*%X;

    nP  = ncol(XX);
    XX  = XX[,-to.be.removed];

    sol.path = lars(XX, YY, type="lasso", normalize=FALSE, intercept=TRUE, max.steps=opt$max.nShifts);

    Tmp = matrix(0, nrow(sol.path$beta), nP);
    Tmp[,-to.be.removed] = sol.path$beta;
    sol.path$beta = Tmp;

    result  = select.best.solution(sol.path, Y, opt);
    eModel  = assign.model(tr, Y, result$shift.placement, opt);

    .print.out.param(eModel, opt$silence);
    return(eModel);
}


.est.shift.placement.known.alpha.multivariate <- function(tr, Y, alpha=0, est.alpha=FALSE, opt){
    library("grplasso");
    library("magic");

    stopifnot( alpha >=0 );

    if ( est.alpha == FALSE ){
        stopifnot(ncol(Y) == length(alpha));
    }

    ## standardazing
    ##NOTE: it turns out that this step is necessary at least for lizard data set.
    if(opt$standardize==TRUE){
        Y  = standardize.matrix(Y);
    }

    nVariables    = ncol(Y);
    X             = generate.design.matrix(tr, "apprX");

    ##NOTE: I used to add a column of 1 to take care of intercept. But it turns out that the new design matrix 
    ##NOTE: cases some problem for grplasso. It takes a lot of memory and it become slower. At some point I should
    ##NOTE: change the solver.

    ##X             = cbind(X,1);
    ##ncolX         = ncol(X);
    to.be.removed = c(ncol(X), which(tr$edge.length < opt$edge.length.threshold));
    ##to.be.removed = c(ncolX-1, which(tr$edge.length < opt$edge.length.threshold));

    offset        = rep(ncol(X)*(0:(ncol(Y)-1)), each=length(to.be.removed));
    to.be.removed = rep(to.be.removed, ncol(Y)) + offset;

    YY  = Y;
    grpX = matrix(0,0,0); #empty matrix;
    for( i in 1:ncol(YY)){
        X  = matrix(0,0,0);
        if ( est.alpha == TRUE ){
            X   = generate.design.matrix(tr, "apprX");
            RE  = cmp.OU.covariance(tr,     alpha = 0 );
        } else {
            X   = generate.design.matrix(tr, "orgX", alpha=alpha[[i]] );
            RE  = cmp.OU.covariance(tr,     alpha = alpha[[i]] );
        }
        Cinvh   = t(RE$D); #\Sigma^{-1/2}
        YY[,i]  = Cinvh%*%YY[,i];
        grpX    = adiag(grpX, Cinvh%*%X);  
    }

    np     = ncol(grpX);
    grpY   = c(YY);
    grpX   = grpX[,-to.be.removed];
    grpIdx = rep(1:ncol(X), ncol(Y))[-to.be.removed];


    ##including the intercept.
    #grpX = cbind(grpX,1);
    #grpIdx = c(grpIdx, NA);

    sol    = run.grplasso(grpX, grpY, nVariables, grpIdx, opt);

    #Tmp                  = matrix(0, np+1, ncol(sol$coefficients));
    Tmp                  = matrix(0, np, ncol(sol$coefficients));
    Tmp[-to.be.removed,] = matrix(sol$coefficients);
    sol$coefficients     = Tmp;

    ##removing the intercept results
    #sol$coefficients     = sol$coefficients[-ncol(grpX), ];

    result  = select.best.solution(sol, Y, opt=opt);
    eModel  = assign.model(tr, Y, result$shift.placement, opt=opt);

    .print.out.param(eModel, opt$silence);
    return(eModel);
}

cmp.uncertainity <- function(tr, model, nItrs = 100, ...){

    library("parallel");
    source("linear_alg_OU_covariance.R");
    RE    = cmp.OU.covariance(tr, alpha=model$alpha);
    C.IH  = t(RE$D);
    C.H   = RE$B;

    ##NOTE: make sure the same Y used for finding shifts is beeing used here too.

    Y     = model$Y;
    YY    = C.IH%*%(Y - model$mu );

    detection.vec = rep(0, nrow(tr$edge));

    #for(itr in 1:nItrs){
    #    YYstar = sample(YY, replace = TRUE);
    #    Ystar  = (C.H%*%YYstar) + model$mu; 
    #    eM     = est.shift.placement(tr, Ystar,  
    #            criterion = model$criterion, 
    #            root.model     = model$root.model, ...);
    #    detection.vec[eM$shift.placement] = detection.vec[eM$shift.placement] + 1;
    #}
    #return(detection.vec/nItrs);

    shift.placement.list = 
        mclapply(X=1:nItrs, FUN=function(itr){

        YYstar = sample(YY, replace = TRUE);
        Ystar  = (C.H%*%YYstar) + model$mu ; 

        eM  <-  tryCatch({
            est.shift.placement(tr, Ystar,  
                                criterion = model$criterion, 
                                root.model     = model$root.model, ...);
        }, error = function(e) {
            print("l1OU error, return NA");
            return(NA); }  );

        if(all(is.na(eM))) {return(NA);}
        return(eM$shift.placement);
    }, mc.cores = 20);

    valid.count <- 0;
    for( i in 1:length(shift.placement.list)){
        if( all(is.na( shift.placement.list[[i]] )) ){
            next;
        }
        valid.count <- valid.count + 1;
        detection.vec[ shift.placement.list[[i]] ] = 
            detection.vec[ shift.placement.list[[i]] ] + 1;
    }

    return(detection.vec/valid.count);
}

cmp.uncertainity.multivariate <- function(tr, model, nItrs = 100, ...){

    source("linear_alg_OU_covariance.R");
    library("parallel");

    Y = as.matrix(model$Y);
    stopifnot( length(model$alpha)     == ncol(Y) );

    YY       = Y;
    C.Hlist  = list();
    for( idx in 1:ncol(Y) ){
        RE    = cmp.OU.covariance(tr, alpha = model$alpha[[idx]] ); 
        C.IH  = t(RE$D); 
        C.Hlist[[idx]] = RE$B;
        YY[, idx]      = C.IH%*%(Y[, idx] - model$mu[ ,idx]);
    }

    detection.vec = rep(0, nrow(tr$edge));

#    for(itr in 1:nItrs){
#
#        Ystar   = YY;
#        idx.vec = sample(1:nrow(YY), replace = TRUE);
#        for( idx in 1:ncol(YY) ){
#            YYstar        = YY[idx.vec, idx];
#            Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx]; 
#        }
#        eM  <-  tryCatch({
#            est.shift.placement(tr, Ystar,  
#                                criterion = model$criterion, ...);
#        }, error = function(e) {
#            print("l1OU error, return NA");
#            return(NA); }  );
#
#        if(all(is.na(eM))) {next;}
#
#        detection.vec[eM$shift.placement] = detection.vec[eM$shift.placement] + 1;
#    }

    shift.placement.list = 
        mclapply(X=1:nItrs, FUN=function(itr){

        Ystar   = YY;
        set.seed( 101 + itr);
        idx.vec = sample(1:nrow(YY), replace = TRUE);
        for( idx in 1:ncol(YY) ){
            YYstar        = YY[idx.vec, idx];
            Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx]; 
        }
        eM  <-  tryCatch({
            est.shift.placement(tr, Ystar,  
                                criterion = model$criterion, ...);
        }, error = function(e) {
            print("l1OU error, return NA");
            return(NA); }  );

        if(all(is.na(eM))) {return(NA);}
        #detection.vec[eM$shift.placement] = detection.vec[eM$shift.placement] + 1;
        return(eM$shift.placement);
    }, mc.cores = 20);

    valid.count <- 0;
    for( i in 1:length(shift.placement.list)){
        if( all(is.na( shift.placement.list[[i]] )) ){
            next;
        }
        valid.count <- valid.count + 1;
        detection.vec[ shift.placement.list[[i]] ] = 
            detection.vec[ shift.placement.list[[i]] ] + 1;
    }

    return(detection.vec/valid.count);
}


generate.design.matrix <- function(tr, type="apprX", alpha){
    library("igraph");

    stopifnot( is.ultrametric(tr) );
    stopifnot( sum( 1:length(tr$tip.label) %in% tr$edge[,1]) == 0);

    nTips  = length(tr$tip.label);
    rNode  = nTips + 1; 
    nEdges = Nedge(tr);

    g  = graph.edgelist(tr$edge, directed = TRUE);
    X  = matrix(0, nTips, nEdges);

    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode = "out", output="epath")$epath;
    ## there must be always a path.
    stopifnot( all( lapply( root2tip, length ) > 0) ); 
    ## since it is ultrametric.
    Tval  = sum(tr$edge.length[root2tip[[1]] ]); 

    if(type == "orgX"){
        for(i in 1:nTips){
            lvec    = c(0, tr$edge.length[root2tip[[i]]]);
            timeVec = Tval - cumsum(lvec);
            timeVec = timeVec[1:length(timeVec)-1];
            X[i, root2tip[[i]] ] = 1 - exp(-alpha*timeVec);
        }
    }else if(type == "apprX"){
        for(i in 1:nTips){
            lvec    = c(0, tr$edge.length[root2tip[[i]]]);
            timeVec = Tval - cumsum(lvec);
            timeVec = timeVec[1:length(timeVec)-1];
            X[i, root2tip[[i]] ] = timeVec;
        }
    }else if(type == "simpX"){
        for(i in 1:nTips){
            X[i, root2tip[[i]] ] = 1;
        }
    }else
        stop("Undefined design matrix type");
    return(X);
}

select.best.solution <- function(sol.path, Y, opt){

    nSols   = get.num.solutions(sol.path);
    stopifnot( nSols > 0 );

    score.vec  = idx.vec = numeric();
    prevshift.placement = NA;
    placement.list = list();
    for(idx in 1:nSols) {

        shift.placement = get.shift.placement(sol.path, idx, Y);
        shift.placement = correct.unidentifiability(tr, shift.placement, opt);

        if ( length(shift.placement) >= opt$max.nShifts    )  { break;}
        if ( setequal(shift.placement, prevshift.placement ) ){ next; }

        score = cmp.model.score(tr, Y, shift.placement, opt);

        placement.list[[idx]] = shift.placement;
        score.vec             = c(score.vec, score);
        idx.vec               = c(idx.vec, idx);

        prevshift.placement   = shift.placement;
    }

    idx.vec   = idx.vec[sort(score.vec, index.return=TRUE)$ix];
    min.score = Inf;   
    min.idx   = NA;

    for( i in 1:min(opt$num.top.placements, length(idx.vec)) ){ 
        if ( is.na(idx.vec[[i]]) ){ break; }
        res = do.backward.selection(tr, Y, placement.list[[ idx.vec[[i]]  ]], opt);
        if ( min.score > res$score){
            min.score       = res$score;
            shift.placement = res$shift.placement;
        }
    }
    return ( list(score=min.score, shift.placement=shift.placement) );
}

do.backward.selection <- function(tr, Y, shift.placement, opt){

    shift.placement = sort(shift.placement, decreasing = TRUE);
    org.score       = cmp.model.score(tr, Y, shift.placement, opt);

    if( length(shift.placement) < 3 ) { 
        return(list(score=org.score, shift.placement=shift.placement)); 
    }  
    for( sp in shift.placement){
        new.placement = setdiff(shift.placement, sp);
        new.score     = cmp.model.score(tr, Y, new.placement, opt);      
        if ( new.score <= org.score){
            shift.placement = new.placement;
            org.score       = new.score;
        }
    }
    return(list(score=org.score, shift.placement=shift.placement));
}


cmp.model.score <-function(tr, Y, shift.placement, opt){

    shift.placement = correct.unidentifiability(tr, shift.placement, opt);

    if(opt$use.saved.scores){
        ##if it's been already computed
        score = get.placement.score.from.list(shift.placement);
        if(!is.na(score)){
            return(score);
        }
    }

    Y       = as.matrix(Y);
    nEdges  = length(tr$edge.length);
    nTips   = length(tr$tip.label);
    nShifts = length(shift.placement);

    if( opt$criterion == "BIC"){
        df.1  = log(nEdges-1)*(nShifts);
        df.2  = log(nTips)*(nShifts + 3);
    } else if( opt$criterion == "AIC"){
        df.1  = 2*nShifts;
        df.2  = 2*3;
    } else if( opt$criterion == "AICc"){
        ## AICc implemented in SURFACE
        p = nShifts + (nShifts + 2)*ncol(Y);
        N = nTips*ncol(Y);
        df.1 = 2*p + (2*p*(p+1))/(N-p-1); 
        ## FIXME I am not sure about the following ...
        if( p > N-2)
            return(Inf);
        df.2 = 0;
    } else if( opt$criterion == "mBIC"){
        res =  cmp.mBIC.df(tr, shift.placement, opt);  
        df.1 = res$df.1;
        df.2 = res$df.2;
    } else if( opt$criterion == "pBICess"){
        score = cmp.pBICess(tr, Y, shift.placement, opt) ;
        if( opt$use.saved.scores){
            add.placement.score.to.list(shift.placement, score);
        }
        return( score );
    } else if( opt$criterion == "pBIC"){
        score = cmp.pBIC(tr, Y, shift.placement, opt) ;
        if( opt$use.saved.scores){
            add.placement.score.to.list(shift.placement, score);
        }
        return( score );
    } 

    score = df.1;
    for( i in 1:ncol(Y)){
        fit   = my.phylolm.interface(tr, Y[,i], shift.placement, opt);
        if ( all( is.na( fit) ) ){
            return(Inf);
            #return(NA);
        } 
        score = score  -2*fit$logLik + df.2;
    }

    if( opt$use.saved.scores){
        add.placement.score.to.list(shift.placement, score);
    }
    return(score);
}

my.phylolm.interface <- function(tr, Y, shift.placement, opt){

    preds = cbind(1, opt$Z[ ,shift.placement]);

    options(warn = -1);
    #fit <- try( phylolm(Y~preds-1, phy=tr, model=opt$root.model) );
    #fit <- try( phylolm( Y~preds-1, phy=tr, model=opt$root.model, 
    #        starting.value = max(1, opt$alpha.lower.bound),
    #        upper.bound    = opt$alpha.upper.bound, 
    #        lower.bound    = opt$alpha.lower.bound ) );  

    fit    <-  try( phylolm(Y~preds-1, phy=tr, model=opt$root.model,
                            lower.bound    = opt$alpha.lower.bound, 
                            upper.bound    = opt$alpha.upper.bound ) );

    options(warn = 0);

    if(class(fit) == "try-error"){ 
      warning( paste0( "phylolm internal error. returning NA; \n
                     num shifts: ", length(shift.placement)) );
      return(NA);
    }

    return(fit);
}

cmp.mBIC.df <- function(tr, shift.placement, opt){
## we assume that tree is of post order.

    shift.placement = sort(shift.placement);
    nTips           = length(tr$tip.label);
    nShifts         = length(shift.placement);

    df.1 =  0; 
    ## pen for the alpha sigma2 and intercept
    df.2 =  3*log(nTips);

    if(nShifts > 0 ){
        ## pen for shift placement
        df.1 = (2*nShifts - 1) *log(nTips);
        ## pen for alpha sigma2 and intercept
        df.2 = 3*log(nTips);

        all.covered.tips = numeric();
        for(eIdx in shift.placement){
            covered.tips = which( opt$Z[,eIdx] > 0 );
            nUniqueTips  = length( setdiff(covered.tips, all.covered.tips) );
            all.covered.tips = union(covered.tips, all.covered.tips);

            ## this must not happen if the input is an 
            ## identifiable placement and the tree is in post order.
            stopifnot( nUniqueTips > 0);
            df.2 = df.2 + log(nUniqueTips); 
        }
        nUniqueTips = length( setdiff(1:nTips, all.covered.tips) );
        df.2 = df.2 + log(nUniqueTips); 
    } 

    return( list(df.1=df.1, df.2=df.2) );
}

cmp.pBICess <- function(tr, Y, shift.placement, opt){

    nShifts = length(shift.placement);
    nEdges  = length(tr$edge[,1]);
    nTips   = length(tr$tip.label);

    df.1  = 2*(nShifts)*log(nEdges-1);
    score = df.1;
    for(i in 1:ncol(Y)){
        fit  = my.phylolm.interface(tr, Y[,i], shift.placement, opt);
        if( all( is.na(fit) ) ){
           #return(NA);
           return(Inf);
        }
        ess  = effective.sample.size(tr, edges=shift.placement, model="OUfixedRoot", 
                 parameters=list(alpha=fit$optpar), FALSE, FALSE);

        df.2  = 3*log(nTips+1) + sum(log(ess+1));
        score = score  -2*fit$logLik + df.2 ;
    }
    return( score );
}

cmp.pBIC <- function(tr, Y, shift.placement, opt){

    nShifts = length(shift.placement);
    nEdges  = length(tr$edge[,1]);
    nTips   = length(tr$tip.label);

    df.1    = 2*(nShifts)*log(nEdges-1);
    score   = df.1;

    for(i in 1:ncol(Y)){
        fit   = my.phylolm.interface(tr, Y[,i], shift.placement, opt);
        if( all( is.na(fit) ) ){
           #return(NA);
           return(Inf);
        } 
        ld    = as.numeric(determinant(fit$vcov * (fit$n - fit$d)/fit$n, log=T)$modulus);
        df.2  = 3*log(nTips) - ld;
        score = score  -2*fit$logLik + df.2 ;
    }
    return( score );
}


assign.model <- function(tr, Y, shift.placement, opt){

    Y       = as.matrix(Y);
    nEdges  = length(tr$edge.length);
    nTips   = length(tr$tip.label);

    mu = alpha = sigma2 = numeric();
    shift.values = optimums = numeric();
    intercept    = optimums.tmp = numeric();

    for(i in 1:ncol(Y)){

        nShifts = length(shift.placement);
        fit     = my.phylolm.interface(tr, as.matrix(Y[,i]), shift.placement, opt);
        if ( all( is.na(fit) ) ){
            stop("model score is NA in assign.model function! this should not happen");
        }

        alpha   = c(alpha,  fit$optpar);
        sigma2  = c(sigma2, fit$sigma2);
        ## E[Y]
        mu      = cbind(mu, fit$fitted.values);

        intercept      = c(intercept, fit$coefficients[[1]]);
        shift.values   = cbind(shift.values, fit$coefficients[2:(nShifts+1)]);

        optimums.tmp = rep(fit$coefficients[[1]], nEdges);
        if( length(shift.placement) > 0 )
            optimums.tmp = convert.shifts2regions(tr, shift.placement, 
                                       fit$coefficients[2:(nShifts+1)]) + fit$coefficients[[1]]; 

        optimums = cbind(optimums, optimums.tmp);
    }
    score = cmp.model.score (tr, Y, shift.placement, opt);

    ##NOTE: adding the trait which used to detect shift positions
    return( list(Y=Y, shift.placement=shift.placement, shift.values=shift.values,
                optimums=optimums, nShifts=length(shift.placement), alpha=alpha, 
                sigma2=sigma2, intercept=intercept, mu = mu, score=score,
                l1ou.opt=opt) );
}

run.grplasso <- function(grpX, grpY, nVariables, grpIdx, opt){


    delta      = opt$grp.delta;
    seq.ub     = opt$grp.seq.ub;

    max.nTries = 5;

    lmbdMax   = 1.2*lambdamax(grpX, grpY, model = LinReg(), index = grpIdx, standardize = FALSE)+1;
    base.seq  = seq(0, seq.ub, delta);

    for( itrTmp in 1:max.nTries){  

        lmbd  = lmbdMax*(0.5^base.seq);
        sink("/dev/null");
        sol = grplasso(grpX, y = grpY, standardize = FALSE, center = FALSE, 
                lambda  = lmbd, model = LinReg(), index = grpIdx,
                control = grpl.control(tol= 0.01));
        sink();

        df.vec = apply(sol$coefficients, 2, function(x) length(which(abs(x)>0))/nVariables );
        df.missing = setdiff(0:max.nShifts, df.vec);

        if( length(df.missing) > 0 ){

            base.seq.idx = seq.tmp = numeric();
            if ( df.missing[[1]] == 0 ){ 
                lmbdMax    = lmbdMax + 2; 
                df.missing = df.missing[-1];
                if ( length( df.missing ) == 0)
                    next;
            }

            for ( mdf in df.missing ){
                if ( length(which(df.vec > mdf)) == 0 )
                    break;
                i = min( which(df.vec > mdf) );
                base.seq.idx = c(base.seq.idx, i);
            }

            idx = 1;
            for ( t in 2:length(base.seq) ){
                if ( idx <= length(base.seq.idx) && t == base.seq.idx[[idx]] ){
                    seq.tmp = c(seq.tmp, seq(base.seq[[t-1]], base.seq[[t]], delta/4) );
                    idx     = idx + 1;
                }
                seq.tmp = c(seq.tmp, base.seq[[t]]);
            }

            #cutting the extra part or adding more
            if( length( which(df.vec > opt$max.nShifts ) ) > 0 ){
                bs.up   = base.seq[ min( which(df.vec > opt$max.nShifts) ) ];
                seq.tmp = seq.tmp [ which(seq.tmp < bs.up) ];
                seq.ub  = min(seq.ub, bs.up+1);
            } else{
                seq.tmp = c(seq.tmp, seq(seq.ub, seq.ub+1, delta) );
                seq.ub  = seq.ub + 1;
            }

            base.seq = seq.tmp;
            if(length(base.seq.idx)>0) 
                delta = delta/4;
        } else { break; }
    }

## actual running of the grplasso
    lmbd  = lmbdMax*(0.5^base.seq);
    sink("/dev/null");
    sol = grplasso(grpX, y = grpY, standardize = FALSE, center = FALSE, 
                   lambda  = lmbd, model = LinReg(), index = grpIdx);
    sink();

    df.missing = setdiff(0:opt$max.nShifts, df.vec);
    for(dfm in df.missing){
        warning( paste0( "There are no solutions with ", dfm , " number of shifts
                        in the solution path of grplasso. you may  want to chanage grp.delta and 
                        grp.seq" ) );
    }
    return(sol);
}
