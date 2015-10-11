# 
#' detects evolutionary shifts in one or multiple traits.
#'
#'@param tree an ultrametric phylogenetic tree of class phylo with branch lengths in postorder.
#'@param Y a trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param max.nShifts an upper bound for the number of shifts. The default value is half the number of tips.
#'@param criterion an information criterion for the model selection.
#'@param root.model an ancestral state model at the root.
#'@param quietly logical. If FALSE, it writes to the output.
#'@param alpha.upper an upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum branch length connected to tips. 
#'@param alpha.lower a lower bound for the phylogenetic adaptation rate.
#'@param standardize logical. If TRUE, the columns of the trait matrix will be standardized.
#'@param num.top.configurations  an internal argument corresponding to the number of the shift configurations that is chosen for further improvement.
#'@param edge.length.threshold the minimum edge length that is considered non-zero.
#'@param grp.delta  an internal parameter. The input lambda sequence for grplasso will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param grp.seq.ub an internal parameter. The input lambda sequence for grplasso will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param l1ou.options if the option object is provided, all the default values will be ignored. 
#'@return 
#' \item{Y}{the input trait vector/matrix.}
#' \item{shift.configuration}{estimated position of shifts, i.e. indices of edges where the estimated shifts occur.}
#' \item{shift.values}{estimates of the shift values.}
#' \item{nShifts}{estimate of number of shifts.}
#' \item{optimums}{the optimum values of the trait along the edges. If the data is multivariate, it is a matrix where each row corresponds to an edge.}
#' \item{alpha}{the maximum likelihood estimate(s) of the adaptation rate \eqn{\alpha}.}
#' \item{sigma2}{the maximum likelihood estimate(s) of the variance rate \eqn{\sigma^2}.}
#' \item{mu}{the fitted values.}
#' \item{residuals}{residuals.}
#' \item{score}{the information criterion value of the estimated shift configuration.}
#' \item{l1ou.options}{the list of options that are used.}
#'
#'@examples
#' 
#' data("lizard.traits", "lizard.tree")
#' Y = lizard.traits[,1]
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' ew <- rep(1, 198) # the tree has 198 edges
#' ew[eModel$shift.configuration] <- 3
#' plot_l1ou(lizard.tree, eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#'@references
#' M. Khabbazian, R. Kriebel, K. Rohe, and Cecile Ane. Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models 
#'
#'@export
estimate_shift_configuration <- function(tree, Y, 
           max.nShifts           = floor(length(tree$tip.label)/2), 
           criterion             = c("pBIC", "pBICess", "mBIC", "BIC", "AIC", "AICc"), 
           root.model            = c("OUrandomRoot", "OUfixedRoot"),
           quietly               = TRUE,
           alpha.upper           = alpha_upper_bound(tree), 
           alpha.lower           = 0,
           standardize           = TRUE,
           num.top.configurations    = max.nShifts/2,
           edge.length.threshold = .Machine$double.eps,
           grp.delta             = 1/16,
           grp.seq.ub            = 5,
           l1ou.options          = NA
     ){

    ## unifying the types
    Y  <-  as.matrix(Y);

    ## setting up options
    if(all(is.na(l1ou.options))){
        l1ou.options                       <-  list();
        l1ou.options$use.saved.scores      <- TRUE;
        l1ou.options$max.nShifts           <- max.nShifts;
        l1ou.options$criterion             <- match.arg(criterion);
        l1ou.options$root.model            <- match.arg(root.model);
        l1ou.options$quietly               <- quietly;
        ##TODO: return warning if estimated alpha is close to its upperbound. What do I mean by close?
        l1ou.options$alpha_upper_bound     <- alpha.upper;
        l1ou.options$alpha.lower.bound     <- alpha.lower;
        l1ou.options$edge.length.threshold <- edge.length.threshold;
        l1ou.options$num.top.configurations    <- num.top.configurations;
        l1ou.options$standardize           <- standardize;
        l1ou.options$grp.seq.ub            <- grp.seq.ub;
        l1ou.options$grp.delta             <- grp.delta;
        l1ou.options$Z                     <- generate_design_matrix(tree, "simpX");
    }



    ##NOTE: checking the assumptions 
    stopifnot(is.ultrametric(tree));
    stopifnot(identical(tree$edge , reorder(tree, "postorder")$edge));
    stopifnot(nrow(Y) == length(tree$tip.label));
    stopifnot(all( row.names(Y) == tree$tip.label));

    if( ncol(Y) == 1 ){  #univariate l1ou
        eModel1 = estimate_shift_configuration_known_alpha(tree, Y, est.alpha=TRUE,      opt=l1ou.options);
        eModel  = estimate_shift_configuration_known_alpha(tree, Y, alpha=eModel1$alpha, opt=l1ou.options);
        ## in case the first step finds a better solution
        if( eModel$score>eModel1$score )
            eModel = eModel1;

    } 
    if( ncol(Y) > 1 ){ #multivariate l1ou
        eModel1 = estimate_shift_configuration_known_alpha_multivariate(tree, Y, est.alpha=TRUE,      opt=l1ou.options);
        eModel  = estimate_shift_configuration_known_alpha_multivariate(tree, Y, alpha=eModel1$alpha, opt=l1ou.options);
        ## in case the first step finds a better solution
        if( eModel$score>eModel1$score )
            eModel = eModel1;
    } 

    ## I really don't like this way of coding but I have to clear the static db object in the cpp code here. 
    ## TODO: can it be implemented differently?
    if( l1ou.options$use.saved.scores){
        erase_configuration_score_db();
    }
    return(eModel);
}

estimate_shift_configuration_known_alpha <- function(tree, Y, alpha=0, est.alpha=FALSE, opt){  

    stopifnot( alpha >=0 );

    #library("lars");
    if ( est.alpha ){ ## BM model
        X   = generate_design_matrix(tree, "apprX");
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=0)$sqrtInvSigma ); 
    } else{           ## OU model
        X   = generate_design_matrix(tree, "orgX", alpha=alpha );
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=alpha)$sqrtInvSigma ); 
    }

    to.be.removed   = c(length(tree$edge.length), which(tree$edge.length < opt$edge.length.threshold));

    YY  = Cinvh%*%Y;
    XX  = Cinvh%*%X;

    nP  = ncol(XX);
    XX  = XX[,-to.be.removed];

    sol.path = lars(XX, YY, type="lasso", normalize=FALSE, intercept=TRUE, max.steps=opt$max.nShifts);

    Tmp = matrix(0, nrow(sol.path$beta), nP);
    Tmp[,-to.be.removed] = sol.path$beta;
    sol.path$beta = Tmp;

    result  = select_best_solution(tree, Y, sol.path, opt);
    eModel  = assign_model(tree, Y, result$shift.configuration, opt);

    print_out(eModel, opt$quietly);
    return(eModel);
}


estimate_shift_configuration_known_alpha_multivariate <- function(tree, Y, alpha=0, est.alpha=FALSE, opt){
    library("grplasso");
    library("magic");

    stopifnot( alpha >=0 );

    if ( est.alpha == FALSE ){
        stopifnot(ncol(Y) == length(alpha));
    }

    ## standardizing
    if(opt$standardize==TRUE){
        Y  = standardize_matrix(Y);
    }

    nVariables    = ncol(Y);
    X             = generate_design_matrix(tree, "apprX");

    ##NOTE: I used to add a column of 1 to take care of intercept. But it turns out that the new design matrix 
    ##NOTE: causes some problem for grplasso(?). It takes a lot of memory and it become slower (have no idea why). 
    ##NOTE: At some point I should update the group lasso solver. 

    ##X             = cbind(X,1);
    ##ncolX         = ncol(X);
    to.be.removed = c(ncol(X), which(tree$edge.length < opt$edge.length.threshold));
    ##to.be.removed = c(ncolX-1, which(tree$edge.length < opt$edge.length.threshold));

    offset        = rep(ncol(X)*(0:(ncol(Y)-1)), each=length(to.be.removed));
    to.be.removed = rep(to.be.removed, ncol(Y)) + offset;

    YY  = Y;
    grpX = matrix(0,0,0); #empty matrix;
    for( i in 1:ncol(YY)){
        X  = matrix(0,0,0);
        if ( est.alpha == TRUE ){
            X   = generate_design_matrix(tree, "apprX");
            RE  = sqrt_OU_covariance(tree,     alpha = 0 );
        } else {
            X   = generate_design_matrix(tree, "orgX", alpha=alpha[[i]] );
            RE  = sqrt_OU_covariance(tree,     alpha = alpha[[i]] );
        }
        Cinvh   = t(RE$sqrtInvSigma); #\Sigma^{-1/2}
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

    sol    = run_grplasso(grpX, grpY, nVariables, grpIdx, opt);

    #Tmp                  = matrix(0, np+1, ncol(sol$coefficients));
    Tmp                  = matrix(0, np, ncol(sol$coefficients));
    Tmp[-to.be.removed,] = matrix(sol$coefficients);
    sol$coefficients     = Tmp;

    ##removing the intercept results
    #sol$coefficients     = sol$coefficients[-ncol(grpX), ];

    result  = select_best_solution(tree, Y, sol, opt=opt);
    eModel  = assign_model(tree, Y, result$shift.configuration, opt=opt);

    print_out(eModel, opt$quietly);
    return(eModel);
}



generate_design_matrix <- function(tree, type="apprX", alpha){
    #library("igraph");

    stopifnot( is.ultrametric(tree) );
    stopifnot( sum( 1:length(tree$tip.label) %in% tree$edge[,1]) == 0);

    nTips  = length(tree$tip.label);
    rNode  = nTips + 1; 
    nEdges = Nedge(tree);

    g  = graph.edgelist(tree$edge, directed = TRUE);
    X  = matrix(0, nTips, nEdges);

    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode = "out", output="epath")$epath;
    ## there must be always a path.
    stopifnot( all( lapply( root2tip, length ) > 0) ); 
    ## since it is ultrametric.
    Tval  = sum(tree$edge.length[root2tip[[1]] ]); 

    if(type == "orgX"){
        for(i in 1:nTips){
            lvec    = c(0, tree$edge.length[root2tip[[i]]]);
            timeVec = Tval - cumsum(lvec);
            timeVec = timeVec[1:length(timeVec)-1];
            X[i, root2tip[[i]] ] = 1 - exp(-alpha*timeVec);
        }
    }else if(type == "apprX"){
        for(i in 1:nTips){
            lvec    = c(0, tree$edge.length[root2tip[[i]]]);
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

select_best_solution <- function(tree, Y, sol.path, opt){

    nSols   = get_num_solutions(sol.path);
    stopifnot( nSols > 0 );

    score.vec  = idx.vec = numeric();
    prevshift.configuration = NA;
    configuration.list = list();
    for(idx in 1:nSols) {

        shift.configuration = get_shift_configuration(sol.path, idx, Y);
        shift.configuration = correct_unidentifiability(tree, shift.configuration, opt);

        if ( length(shift.configuration) >= opt$max.nShifts    )  { break;}
        if ( setequal(shift.configuration, prevshift.configuration ) ){ next; }

        score = cmp_model_score(tree, Y, shift.configuration, opt);

        configuration.list[[idx]] = shift.configuration;
        score.vec             = c(score.vec, score);
        idx.vec               = c(idx.vec, idx);

        prevshift.configuration   = shift.configuration;
    }

    idx.vec   = idx.vec[sort(score.vec, index.return=TRUE)$ix];
    min.score = Inf;   
    min.idx   = NA;

    for( i in 1:min(opt$num.top.configurations, length(idx.vec)) ){ 
        if ( is.na(idx.vec[[i]]) ){ break; }
        res = do_backward_selection(tree, Y, configuration.list[[ idx.vec[[i]]  ]], opt);
        if ( min.score > res$score){
            min.score       = res$score;
            shift.configuration = res$shift.configuration;
        }
    }
    return ( list(score=min.score, shift.configuration=shift.configuration) );
}

do_backward_selection <- function(tree, Y, shift.configuration, opt){

    shift.configuration = sort(shift.configuration, decreasing = TRUE);
    org.score       = cmp_model_score(tree, Y, shift.configuration, opt);

    if( length(shift.configuration) < 3 ) { 
        return(list(score=org.score, shift.configuration=shift.configuration)); 
    }  
    for( sp in shift.configuration){
        new.configuration = setdiff(shift.configuration, sp);
        new.score     = cmp_model_score(tree, Y, new.configuration, opt);      
        if ( new.score <= org.score){
            shift.configuration = new.configuration;
            org.score       = new.score;
        }
    }
    return(list(score=org.score, shift.configuration=shift.configuration));
}


#
#' computes the information criterion score for the given configuration
#'
#'@param tree an ultrametric phylogenetic tree of class phylo with branch lengths in postorder.
#'@param Y a trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param shift.configuration shift positions, i.e. indices of the edges where the shifts occur.
#'@param criterion an information criterion.
#'@param root.model an ancestral state model at the root.
#'
#'@return the information criterion value of the shift configuration.
#'
#'@examples
#' 
#' library("l1ou") 
#' data("lizard.traits", "lizard.tree")
#' Y <- lizard.traits[,1] 
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' ic.score <- configuration_ic(lizard.tree, eModel$Y, eModel$shift.configuration, criterion="pBIC")
#' print(ic.score)
#'
#'@seealso \code{\link{estimate_shift_configuration}} 
#'
#'@references
#' M. Khabbazian, R. Kriebel, K. Rohe, and Cecile Ane. Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models 
#'@export
configuration_ic <- function(tree, Y, shift.configuration, 
                     criterion = c("pBIC", "pBICess", "mBIC", "BIC", "AIC", "AICc"), 
                     root.model = c("OUrandomRoot", "OUfixedRoot")
                     ){
    opt = list();

    opt$criterion         <- match.arg(criterion);
    opt$root.model        <- match.arg(root.model);
    opt$alpha_upper_bound <- alpha_upper_bound(tree);
    opt$alpha.lower.bound <- 0;
    opt$Z                 <- generate_design_matrix(tree, "simpX");
    opt$use.saved.scores  <- FALSE;

    score = cmp_model_score(tree, Y, shift.configuration, opt);
    return(score);
}

cmp_model_score <-function(tree, Y, shift.configuration, opt){

    shift.configuration = correct_unidentifiability(tree, shift.configuration, opt);

    if(opt$use.saved.scores){
        ##if it's been already computed
        score = get_configuration_score_to_list(shift.configuration);
        if(!is.na(score)){
            return(score);
        }
    }
    ic = opt$criterion;

    Y       = as.matrix(Y);
    nEdges  = length(tree$edge.length);
    nTips   = length(tree$tip.label);
    nShifts = length(shift.configuration);

    if( ic == "BIC"){
        df.1  = log(nEdges-1)*(nShifts);
        df.2  = log(nTips)*(nShifts + 3);
    } else if( ic == "AIC"){
        df.1  = 2*nShifts;
        df.2  = 2*3;
    } else if( ic == "AICc"){
        ## AICc implemented in SURFACE
        p = nShifts + (nShifts + 2)*ncol(Y);
        N = nTips*ncol(Y);
        df.1 = 2*p + (2*p*(p+1))/(N-p-1); 
        ## FIXME I am not sure about the following ...
        if( p > N-2)
            return(Inf);
        df.2 = 0;
    } else if( ic == "mBIC"){
        res =  cmp_mBIC_df(tree, shift.configuration, opt);  
        df.1 = res$df.1;
        df.2 = res$df.2;
    } else if( ic == "pBICess"){
        score = cmp_pBICess(tree, Y, shift.configuration, opt) ;
        if( opt$use.saved.scores){
            add_configuration_score_to_list(shift.configuration, score);
        }
        return( score );
    } else if( ic == "pBIC"){
        score = cmp_pBIC(tree, Y, shift.configuration, opt) ;
        if( opt$use.saved.scores){
            add_configuration_score_to_list(shift.configuration, score);
        }
        return( score );
    } 

    score = df.1;
    for( i in 1:ncol(Y)){
        fit   = my_phylolm_interface(tree, Y[,i], shift.configuration, opt);
        if ( all( is.na( fit) ) ){
            return(Inf);
            #return(NA);
        } 
        score = score  -2*fit$logLik + df.2;
    }

    if( opt$use.saved.scores){
        add_configuration_score_to_list(shift.configuration, score);
    }
    return(score);
}

my_phylolm_interface <- function(tree, Y, shift.configuration, opt){

    preds = cbind(1, opt$Z[ ,shift.configuration]);

    options(warn = -1);
    #fit <- try( phylolm(Y~preds-1, phy=tree, model=opt$root.model) );
    #fit <- try( phylolm( Y~preds-1, phy=tree, model=opt$root.model, 
    #        starting.value = max(1, opt$alpha.lower.bound),
    #        upper.bound    = opt$alpha_upper_bound, 
    #        lower.bound    = opt$alpha.lower.bound ) );  

    fit    <-  try( phylolm(Y~preds-1, phy=tree, model=opt$root.model,
                            lower.bound    = opt$alpha.lower.bound, 
                            upper.bound    = opt$alpha_upper_bound ), silent = TRUE);
    options(warn = 0);

    if(class(fit) == "try-error"){ 
      warning( paste0( "phylolm internal error. returning NA; num shifts: ", length(shift.configuration)) );
      return(NA);
    }

    return(fit);
}

cmp_mBIC_df <- function(tree, shift.configuration, opt){
## we assume that tree is of post order.

    shift.configuration = sort(shift.configuration);
    nTips           = length(tree$tip.label);
    nShifts         = length(shift.configuration);

    df.1 =  0; 
    ## pen for the alpha sigma2 and intercept
    df.2 =  3*log(nTips);

    if(nShifts > 0 ){
        ## pen for shift configuration
        df.1 = (2*nShifts - 1) *log(nTips);
        ## pen for alpha sigma2 and intercept
        df.2 = 3*log(nTips);

        all.covered.tips = numeric();
        for(eIdx in shift.configuration){
            covered.tips = which( opt$Z[,eIdx] > 0 );
            nUniqueTips  = length( setdiff(covered.tips, all.covered.tips) );
            all.covered.tips = union(covered.tips, all.covered.tips);

            ## this must not happen if the input is an 
            ## identifiable configuration and the tree is in post order.
            stopifnot( nUniqueTips > 0);
            df.2 = df.2 + log(nUniqueTips); 
        }
        nUniqueTips = length( setdiff(1:nTips, all.covered.tips) );
        df.2 = df.2 + log(nUniqueTips); 
    } 

    return( list(df.1=df.1, df.2=df.2) );
}

cmp_pBICess <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration);
    nEdges  = length(tree$edge[,1]);
    nTips   = length(tree$tip.label);

    df.1  = 2*(nShifts)*log(nEdges-1);
    score = df.1;
    for(i in 1:ncol(Y)){
        fit  = my_phylolm_interface(tree, Y[,i], shift.configuration, opt);
        if( all( is.na(fit) ) ){
           #return(NA);
           return(Inf);
        }
        ess  = effective.sample.size(tree, edges=shift.configuration, model="OUfixedRoot", 
                 parameters=list(alpha=fit$optpar), FALSE, FALSE);

        df.2  = 3*log(nTips+1) + sum(log(ess+1));
        score = score  -2*fit$logLik + df.2 ;
    }
    return( score );
}

cmp_pBIC <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration);
    nEdges  = length(tree$edge[,1]);
    nTips   = length(tree$tip.label);

    df.1    = 2*(nShifts)*log(nEdges-1);
    score   = df.1;

    for(i in 1:ncol(Y)){
        fit   = my_phylolm_interface(tree, Y[,i], shift.configuration, opt);
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


assign_model <- function(tree, Y, shift.configuration, opt){

    Y       = as.matrix(Y);
    nEdges  = length(tree$edge.length);
    nTips   = length(tree$tip.label);

    resi = mu = alpha = sigma2 = numeric();
    shift.values = optimums = numeric();
    intercept    = optimums.tmp = numeric();

    for(i in 1:ncol(Y)){

        nShifts = length(shift.configuration);
        fit     = my_phylolm_interface(tree, as.matrix(Y[,i]), shift.configuration, opt);
        if ( all( is.na(fit) ) ){
            stop("model score is NA in assign_model function! This should not happen");
        }

        alpha   = c(alpha,  fit$optpar);
        sigma2  = c(sigma2, fit$sigma2);
        ## E[Y]
        mu      = cbind(mu, fit$fitted.values);
        resi    = cbind(resi, fit$residuals);

        intercept      = c(intercept, fit$coefficients[[1]]);
        shift.values   = cbind(shift.values, fit$coefficients[2:(nShifts+1)]);

        optimums.tmp = rep(fit$coefficients[[1]], nEdges);
        if( length(shift.configuration) > 0 )
            optimums.tmp = convert_shifts2regions(tree, shift.configuration, 
                                       fit$coefficients[2:(nShifts+1)]) + fit$coefficients[[1]]; 

        optimums = cbind(optimums, optimums.tmp);
    }
    score = cmp_model_score (tree, Y, shift.configuration, opt);

    ##NOTE: adding the trait which used to detect shift positions
    return( list(Y=Y, 
                 shift.configuration=shift.configuration, 
                 shift.values=shift.values,
                 nShifts=length(shift.configuration), 
                 optimums=optimums, 
                 alpha=alpha, 
                 sigma2=sigma2, 
                 intercept=intercept, 
                 mu = mu, 
                 residuals = resi,
                 score=score,
                 l1ou.options=opt) );
}


function (grpx, grpy, nvariables, grpidx, opt) {
    delta  = opt$grp.delta
    seq.ub = opt$grp.seq.ub
    max.ntries = 10;
    lmbdmax = 1.2 * lambdamax(grpx, grpy, model = linreg(), index = grpidx, standardize = false) + 1

    base.seq = seq(0, seq.ub, delta)
    for (itrtmp in 1:max.ntries) {
        lmbd = lmbdmax * (0.5^base.seq)
        sink("/dev/null")
        sol = grplasso(grpx, y = grpy, standardize = false, center = false, 
            lambda = lmbd, model = linreg(), index = grpidx, 
            control = grpl.control(tol = 0.01))
        sink()
        df.vec = apply(sol$coefficients, 2, function(x) length(which(rowsums(matrix(abs(x), ncol=nvariables)) > 0)))
        df.missing = setdiff(0:(opt$max.nshifts+1), df.vec)
        names(df.vec) <- null;
        print(df.vec);

        cutextra = true;
        if (length(df.missing) > 0) {
            base.seq.tmp =numeric();
            idx = prev.idx = 1;
            for (mdf in df.missing) {

                if( length( which(df.vec > mdf) ) > 0 ){ ##that means at the begining or in the middle of the sequence
                    idx = min(which(df.vec > mdf))
                    if( idx == 1){ ## that means we should add to the begining
                        ##todo: i don't like it
                        lmbdmax    = lmbdmax + 2
                        next;
                    }
                    if( prev.idx == idx){ ## not to repeat the same sequence.
                        next;
                    }
                    lower = base.seq[[idx - 1]] + delta/4;
                    upper = base.seq[[idx]]
                    base.seq.tmp=c( base.seq.tmp, base.seq[prev.idx:idx], seq(lower, upper, delta/4));
                } else{ ## that means we should add to the end
                    lower = base.seq.tmp[[length(base.seq.tmp)]]+delta;
                    upper = max(lower, max(base.seq)) + 1;
                    base.seq.tmp = c(base.seq.tmp, seq(lower, upper, delta));
                    cutextra  = false;
                }
                prev.idx = idx;
            }

            if( cutextra ){
                indices = which(df.vec > (opt$max.nshifts+4));
                if( length(indices) > 0 ){
                    upper.idx = min(indices);
                }else{
                    upper.idx = length(base.seq);
                }

                if( idx < upper.idx){
                    base.seq.tmp = c(base.seq.tmp, base.seq[idx:upper.idx] ); 
                }
            }
            delta = delta/4
            base.seq = unique(sort(base.seq.tmp));
        }else { break }
    }
    lmbd = lmbdmax * (0.5^base.seq)
    sink("/dev/null")
    sol = grplasso(grpx, y = grpy, standardize = false, center = false, 
        lambda = lmbd, model = linreg(), index = grpidx)
    sink()
    df.missing = setdiff(0:opt$max.nshifts, df.vec)
    for (dfm in df.missing) {
        warning(paste0("there are no solutions with ", dfm, " number of shifts  
                in the solution path of grplasso. you may want to change grp.delta and grp.seq"))
    }
    return(sol)
}



#run_grplasso <- function(grpX, grpY, nVariables, grpIdx, opt){
#
#
#    delta      = opt$grp.delta;
#    seq.ub     = opt$grp.seq.ub;
#
#    max.nTries = 5;
#
#    lmbdMax   = 1.2*lambdamax(grpX, grpY, model = LinReg(), index = grpIdx, standardize = FALSE)+1;
#    base.seq  = seq(0, seq.ub, delta);
#
#    for( itrTmp in 1:max.nTries){  
#
#        lmbd  = lmbdMax*(0.5^base.seq);
#        sink("/dev/null");
#        sol = grplasso(grpX, y = grpY, standardize = FALSE, center = FALSE, 
#                lambda  = lmbd, model = LinReg(), index = grpIdx,
#                control = grpl.control(tol= 0.01));
#        sink();
#
#        df.vec = apply(sol$coefficients, 2, function(x) length(which(abs(x)>0))/nVariables );
#        df.missing = setdiff(0:opt$max.nShifts, df.vec);
#
#        if( length(df.missing) > 0 ){
#
#            base.seq.idx = seq.tmp = numeric();
#            if ( df.missing[[1]] == 0 ){ 
#                lmbdMax    = lmbdMax + 2; 
#                df.missing = df.missing[-1];
#                if ( length( df.missing ) == 0)
#                    next;
#            }
#
#            for ( mdf in df.missing ){
#                if ( length(which(df.vec > mdf)) == 0 )
#                    break;
#                i = min( which(df.vec > mdf) );
#                base.seq.idx = c(base.seq.idx, i);
#            }
#
#            idx = 1;
#            for ( t in 2:length(base.seq) ){
#                if ( idx <= length(base.seq.idx) && t == base.seq.idx[[idx]] ){
#                    seq.tmp = c(seq.tmp, seq(base.seq[[t-1]], base.seq[[t]], delta/4) );
#                    idx     = idx + 1;
#                }
#                seq.tmp = c(seq.tmp, base.seq[[t]]);
#            }
#
#            #cutting the extra part or adding more
#            if( length( which(df.vec > opt$max.nShifts ) ) > 0 ){
#                bs.up   = base.seq[ min( which(df.vec > opt$max.nShifts) ) ];
#                seq.tmp = seq.tmp [ which(seq.tmp < bs.up) ];
#                seq.ub  = min(seq.ub, bs.up+1);
#            } else{
#                seq.tmp = c(seq.tmp, seq(seq.ub, seq.ub+1, delta) );
#                seq.ub  = seq.ub + 1;
#            }
#
#            base.seq = seq.tmp;
#            if(length(base.seq.idx)>0) 
#                delta = delta/4;
#        } else { break; }
#    }
#
### actual running of the grplasso
#    lmbd  = lmbdMax*(0.5^base.seq);
#    sink("/dev/null");
#    sol = grplasso(grpX, y = grpY, standardize = FALSE, center = FALSE, 
#                   lambda  = lmbd, model = LinReg(), index = grpIdx);
#    sink();
#
#    df.missing = setdiff(0:opt$max.nShifts, df.vec);
#    for(dfm in df.missing){
#        warning( paste0( "There are no solutions with ", dfm , " number of shifts
#                        in the solution path of grplasso. You may want to change grp.delta and 
#                        grp.seq" ) );
#    }
#    return(sol);
#}
