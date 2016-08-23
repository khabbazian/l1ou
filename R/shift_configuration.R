# 
#' Detects evolutionary shifts under an OU model
#'
#'This function takes in one or multiple traits, and automatically detects the phylogenetic placement and 
#'the magnitude of shifts in the evolution of these traits. The model assumes an Ornstein-Uhlenbeck process
#'whose parameters are estimated (adaptation `strength' \eqn{\alpha}{alpha} and drift variance \eqn{\sigma^2}{sigma^2}).
#'Instantaneous shifts in the optimal trait value affect the traits over time.
#'
#'@param tree ultrametric tree of class phylo with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param max.nShifts upper bound for the number of shifts. The default value is half the number of tips.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@param root.model ancestral state model at the root.
#'@param candid.edges a vector of indices of candidate edges where the shifts may occur. If provided, shifts will only be allowed on these edges; otherwise all edges will be considered.
#'@param quietly logical. If FALSE, a basic summary of the progress and results is printed.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum branch length connected to tips. 
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param lars.alg model selection algorithm for LARS in univariate case. 
#'@param nCores number of processes to be created for parallel computing. If nCores=1 then it will run sequentially. Otherwise, it creates nCores processes by using mclapply function. For parallel computing it, requires parallel package.
#'@param rescale logical. If TRUE, the columns of the trait matrix are first rescaled so that all have the same l2-norm. If TRUE, the scores will be based on the rescale one.
#'@param edge.length.threshold minimum edge length that is considered non-zero. Branches with shorter length are considered as soft polytomies, disallowing shifts on such branches.
#'@param grp.delta internal (used when the data contain multiple traits). The input lambda sequence for the group lasso, in `grplasso', will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param grp.seq.ub (used for multiple traits). The input lambda sequence for grplasso will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param l1ou.options if provided, all the default values will be ignored. 
#'@return 
#' \item{Y}{input trait vector/matrix.}
#' \item{tree}{input tree.}
#' \item{shift.configuration}{estimated shift positions, i.e. vector of indices of edges where the estimated shifts occur.}
#' \item{shift.values}{estimates of the shift values.}
#' \item{nShifts}{estimated number of shifts.}
#' \item{optima}{optimum values of the trait at tips. If the data are multivariate, this is a matrix where each row corresponds to a tip.}
#' \item{edge.optima}{optimum values of the trait on the edges. If the data are multivariate, this is a matrix where each row corresponds to an edge.}
#' \item{alpha}{maximum likelihood estimate(s) of the adaptation rate \eqn{\alpha}{alpha}, one per trait.}
#' \item{sigma2}{maximum likelihood estimate(s) of the variance rate \eqn{\sigma^2}{sigma^2}, one per trait.}
#' \item{mu}{fitted values, i.e. estimated trait means.}
#' \item{residuals}{residuals. These residuals are phylogenetically correlated.}
#' \item{score}{information criterion value of the estimated shift configuration.}
#' \item{profile}{list of shift configurations sorted by their ic scores.}
#' \item{l1ou.options}{list of options that were used.}
#'
#'@details
#'For information criteria: see \code{\link{configuration_ic}}. 
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' # here lizard.traits already has row names:
#' rownames(lizard.traits)
#' # also, it is a matrix (not data frame) so columns retain row names:
#' names(lizard.traits[,1])
#' # If your trait data "dat" does not have row names but instead has
#' # species names in a column called "species", then you can
#' # create row names containing the species names like this:
#' # rownames(dat) <- dat$species
#' lizard <- adjust_data(lizard.tree, lizard.traits[,1])
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y)
#' eModel
#'  
#' ## use parallel computing to accelerate the computation
#' eModel.par <- estimate_shift_configuration(lizard$tree, lizard$Y, nCores=8)
#'
#' stopifnot( identical( sort(eModel.par$shift.configuration), sort(eModel$shift.configuration) ) ) ## TRUE
#'
#' nEdges <- Nedge(lizard.tree) # total number of edges
#' ew <- rep(1,nEdges)  # to set default edge width of 1
#' ew[eModel$shift.configuration] <- 3   # to widen edges with a shift 
#' plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#' # example to constrain the set of candidate branches with a shift
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y, criterion="AICc")
#' ce <- eModel$shift.configuration # set of candidate edges
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y, candid.edges = ce)
#' plot(eModel, edge.ann.cex=0.7, cex=0.5, label.offset=0.02)
#'
#'@references
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
estimate_shift_configuration <- function(tree, Y, 
           max.nShifts            = floor(length(tree$tip.label)/2), 
           criterion              = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
           root.model             = c("OUfixedRoot", "OUrandomRoot"),
           candid.edges           = NA,
           quietly                = TRUE,
           alpha.starting.value   = NA, 
           alpha.upper            = alpha_upper_bound(tree), 
           alpha.lower            = NA,
           lars.alg               = c("lasso", "stepwise"),
           nCores                 = 1,
           rescale                = TRUE,
           edge.length.threshold  = .Machine$double.eps,
           grp.delta              = 1/16,
           grp.seq.ub             = 5,
           l1ou.options           = NA
     ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if (is.null(tree$edge.length)) stop("the tree has no branch lengths.")
    if (is.null(tree$tip.label))   stop("the tree has no tip labels.")	
    if(!is.ultrametric(tree))      stop("the input tree is not ultrametric.")

    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the tree is not in postorder, use adjust_data function to reorder the tree!")

    Y <- as.matrix(Y)

    multivariate.missing <- FALSE
    if( any(is.na(Y)) ){
        if( ncol(Y) == 1){
            stop("some of the entries of the trait vector (Y) are missing.
                 you may use drop.tip to drop corresponding tips from the tree.\n")
        }else{
            ## the trait matrix can have some missing values as long as
            ## all the tips have at least one value.
            warning("some of the entries of the trait matrix (Y) are
                    missing.\n", immediate.=TRUE)
            if( length( which(rowSums(!is.na(Y))==0) ) != 0 ){
                stop("the trait matrix has a row with all missing values. you
                     may use drop.tip to drop corresponding tips from the tree.\n")
            }
            multivariate.missing <- TRUE
        }
    }

    if( class(Y) != "matrix"){
        Y <- as.matrix(Y)
    }

    if( nrow(Y) != length(tree$tip.label)){
       stop("the number of entries/rows of the trait vector/matrix (Y) 
            doesn't match the number of tips.\n") 
    }

    if( is.null(rownames(Y)) ){
        warning("no names provided for the trait(s) entries/rows. So it is assumed that 
                entries/rows match with the tip labels in the same order.\n", immediate.=TRUE)
        rownames(Y)  <- tree$tip.label
    } else{
        if( any(is.na(rownames(Y))) ){
            stop("some of the names in either trait vector/matrix or tree tip 
                 labels are unavailable.\n")
        }
    }

    if(!identical(rownames(Y), tree$tip.label)){
        diffres = setdiff(rownames(Y), tree$tip.label)
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the tip labels of the input tree.\n")
        }
        diffres = setdiff(tree$tip.label, rownames(Y))
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the input trait. You may want to use drop.tip(tree, setdiff(tree$tip.label,rownames(Y))) 
                 to drop extra tips in the tree.\n")
        }

        stop("the order of entries/rows of the trait vector/matrix (Y) does not matche 
             the order of the tip labels. Use adjust_data function to fix that.\n")
    }

    stopifnot(all(rownames(Y) == tree$tip.label))
    stopifnot(identical(rownames(Y), tree$tip.label))

    if( max.nShifts > length(tree$tip.label)){
        warning("max.nShifts should be a positive number less than number of tips. I set it to number of tips.\n")
        max.nShifts  <-  length(tree$tip.label)
    }
    if( max.nShifts < 0){
        warning("max.nShifts should be a positive number less than number of tips. I set it to 0.\n")
        max.nShifts  <- 1 
    }

    if(!is.na(alpha.upper))
        if( alpha.upper <= 0 ){
            stop("alpha.upper must be strictly positive.\n")
        }

    if(!is.na(alpha.lower))
        if( alpha.lower < 0 ){
            warning("alpha.lower must be greater than zero. I set it to zero.\n")
            alpha.lower <- 0 
        }
    if(!is.na(alpha.lower))
        if( alpha.upper < alpha.lower ){
            warning("alpha.upper must be equal or greater than alpha.lower. I set them equal.\n")
            alpha.upper  <- alpha.lower
        }

    stopifnot( nCores > 0 )
    parallel.computing <- FALSE
    if(nCores>1){
        if(!require(parallel)){
            warning("parallel package is not available. The process will run sequentially.", immediate=TRUE)
            nCores <- 1
        }else{
            parallel.computing <- TRUE
        }
    }

    if (all(is.na(l1ou.options))){
        l1ou.options                   <- list()
        l1ou.options$nCores             <- nCores
        l1ou.options$parallel.computing <- parallel.computing
        ##NOTE: The saving_score functions are unprotected. To avoid race,
        ## I simply disable them in parallel mode until later that I figure out how to fix it.
        if(parallel.computing)
            l1ou.options$use.saved.scores  <- FALSE
        else
            l1ou.options$use.saved.scores  <- TRUE

        l1ou.options$max.nShifts       <- max.nShifts
        l1ou.options$criterion         <- match.arg(criterion)
        l1ou.options$lars.alg          <- match.arg(lars.alg)
        l1ou.options$root.model        <- match.arg(root.model)
        l1ou.options$quietly           <- quietly

        l1ou.options$alpha.starting.value   <- alpha.starting.value
        l1ou.options$alpha.upper.bound      <- alpha.upper
        l1ou.options$alpha.lower.bound      <- alpha.lower
        l1ou.options$edge.length.threshold  <- edge.length.threshold
        l1ou.options$rescale                <- rescale

        l1ou.options$grp.seq.ub   <- grp.seq.ub
        l1ou.options$grp.delta    <- grp.delta
        l1ou.options$candid.edges <- candid.edges
        l1ou.options$Z            <- generate_design_matrix(tree, "simpX")
        ## each tree in tree.list represents a trait where the tips corresponding
        ## to NA values in the trail have been dropped.
        l1ou.options$multivariate.missing <- multivariate.missing
        if(multivariate.missing)
            l1ou.options$tree.list <- gen_tree_array(tree, Y)
        else
            l1ou.options$tree.list <- NULL
    }

    if(length(l1ou.options$candid.edges)==0 || l1ou.options$max.nShifts == 0){
        l1ou.options$use.saved.scores <- FALSE
        return( fit_OU(tree, Y, shift.configuration=c(), l1ou.options=l1ou.options) )
    }

    if (l1ou.options$use.saved.scores) { erase_configuration_score_db() }

    cat("Starting first LASSO (alpha=0) to find a list of candidate configurations.\n")
    if (ncol(Y) == 1) {
        eModel1 = estimate_shift_configuration_known_alpha(tree, 
            Y, est.alpha = TRUE, opt = l1ou.options)
        cat("Starting second LASSO (alpha=",round(eModel1$alpha,2),") for another list of candidates.\n")
        eModel = estimate_shift_configuration_known_alpha(tree, 
            Y, alpha = eModel1$alpha, opt = l1ou.options)
        if (eModel$score > eModel1$score) 
            eModel = eModel1
    }
    if (ncol(Y) > 1) {
        eModel1 = estimate_shift_configuration_known_alpha_multivariate(tree, 
            Y, est.alpha = TRUE, opt = l1ou.options)
        cat("Starting second LASSO (alpha=",round(eModel1$alpha,2),") for another list of candidates.\n")
        eModel = estimate_shift_configuration_known_alpha_multivariate(tree, 
            Y, alpha = eModel1$alpha, opt = l1ou.options)
        if (eModel$score > eModel1$score) 
            eModel = eModel1
    }

    eModel$profile = list_investigated_configs() 

    if (l1ou.options$use.saved.scores) { erase_configuration_score_db() }

    return(eModel)
}


estimate_shift_configuration_known_alpha <- function(tree, Y, alpha=0, est.alpha=FALSE, opt){  

    stopifnot( alpha >=0 )

    if ( est.alpha ){ ## BM model
        X   = generate_design_matrix(tree, "apprX")
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=0, root.model = "OUfixedRoot", 
                                        check.order=F, check.ultrametric=F)$sqrtInvSigma ) 
    } else{           ## OU model
        X   = generate_design_matrix(tree, "orgX", alpha=alpha )
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=alpha, root.model = opt$root.model, 
                                        check.order=F, check.ultrametric=F)$sqrtInvSigma ) 
    }

    if(!all(is.na(opt$candid.edges))){
        to.be.removed  = setdiff(1:length(tree$edge.length), opt$candid.edges)
    }else{
        to.be.removed  = c(length(tree$edge.length), which(tree$edge.length < opt$edge.length.threshold))
    }

    YY  = Cinvh%*%Y
    XX  = Cinvh%*%X

    nP  = ncol(XX)
    XX  = as.matrix(XX[,-to.be.removed])

    capture.output(
            sol.path  <- lars(XX, YY, type=opt$lars.alg, normalize=FALSE,
                              intercept=TRUE, max.steps=opt$max.nShifts)
        )

    Tmp = matrix(0, nrow(sol.path$beta), nP)
    Tmp[,-to.be.removed] = sol.path$beta
    sol.path$beta = Tmp

    result  = select_best_solution(tree, Y, sol.path, opt)
    eModel  = fit_OU_model(tree, Y, result$shift.configuration, opt)

    if(!opt$quietly){
        print(eModel)
        cat("-------\n")
    }
    return(eModel)
}


estimate_shift_configuration_known_alpha_multivariate <- function(tree, Y, alpha=0, est.alpha=FALSE, opt){

    stopifnot( alpha >=0 )

    if ( est.alpha == FALSE ){
        stopifnot(ncol(Y) == length(alpha))
    }

    Ymv <- Y
    if(opt$rescale==TRUE){
        Ymv  = rescale_matrix(Ymv)
    }

    nVariables    = ncol(Ymv)
    nEdges        = Nedge(tree)
    ##X             = generate_design_matrix(tree, "apprX")
    ##X             = cbind(X,1)
    ##ncolX         = ncol(X)
    ##to.be.removed = c(ncolX-1, which(tree$edge.length < opt$edge.length.threshold))

    if(!all(is.na(opt$candid.edges))){
        to.be.removed  = setdiff(1:length(tree$edge.length), opt$candid.edges)
    }else{
        to.be.removed = c(nEdges, which(tree$edge.length < opt$edge.length.threshold))
    }

    offset        = rep(nEdges*(0:(ncol(Ymv)-1)), each=length(to.be.removed))
    to.be.removed = rep(to.be.removed, ncol(Ymv)) + offset

    YY  = Ymv
    grpX = matrix(0,0,0) #empty matrix
    for( i in 1:ncol(Ymv)){
        X     = matrix(0,0,0)

        if ( est.alpha == TRUE ){
            X   = generate_design_matrix(tree, "apprX")
            RE  = sqrt_OU_covariance(tree, root.model = "OUfixedRoot",
                                     alpha = 0, check.order=F, check.ultrametric=F )
        } else {
            X   = generate_design_matrix(tree, "orgX", alpha=alpha[[i]])
            RE  = sqrt_OU_covariance(tree,  root.model = opt$root.model,   
                                     alpha = alpha[[i]], 
                                     check.order=F, check.ultrametric=F )
        }
        Cinvh   = t(RE$sqrtInvSigma) #\Sigma^{-1/2}

        y.ava        = !is.na(Ymv[,i])
        YY[y.ava, i] = Cinvh[y.ava, y.ava] %*% Ymv[y.ava, i]

        ##X     = cbin(Cinvh%*%X,1)
        ##grpX  = adiag(grpX, X)  
        grpX    = adiag(grpX, Cinvh%*%X)  
    }

    np     = ncol(grpX)
    grpY   = c(YY)
    grpX   = as.matrix(grpX[,-to.be.removed])
    grpIdx = rep(1:ncol(X), ncol(Ymv))[-to.be.removed]
    ##grpIdx = rep(c(1:(ncol(X)-1),NA), ncol(Ymv))[-to.be.removed]


    ###NOTE: handling NAs in grpY
    grpY.ava = !is.na(grpY)
    grpY     = grpY[grpY.ava]
    grpX     = grpX[grpY.ava, ]

    grpX.col.nZero.idx = which(colSums(abs(grpX))!=0)
    grpX.nCol          = ncol(grpX)
    grpX               = grpX[,grpX.col.nZero.idx]
    grpIdx             = grpIdx[grpX.col.nZero.idx]

    sol = run_grplasso(grpX, grpY, nVariables, grpIdx, opt)

    Tmp                      = matrix(0, grpX.nCol, ncol(sol$coefficients))
    Tmp[grpX.col.nZero.idx,] = matrix(sol$coefficients)
    sol$coefficients         = Tmp
    ###end NOTE:
    
    Tmp                  = matrix(0, np, ncol(sol$coefficients))
    Tmp[-to.be.removed,] = matrix(sol$coefficients)
    sol$coefficients     = Tmp

    ##removing the intercept results
    #sol$coefficients     = sol$coefficients[-ncol(grpX), ]

    ### use the original Y
    result  = select_best_solution(tree, Y, sol, opt=opt)
    eModel  = fit_OU_model(tree, Y, result$shift.configuration, opt=opt)

    if(!opt$quietly){
        print(eModel)
        print("-------")
    }
    return(eModel)
}



generate_design_matrix <- function(tree, type="apprX", alpha){
    stopifnot( is.ultrametric(tree) )
    stopifnot( sum( 1:length(tree$tip.label) %in% tree$edge[,1]) == 0)

    nTips  = length(tree$tip.label)
    rNode  = nTips + 1 
    nEdges = Nedge(tree)

    g  = graph.edgelist(tree$edge, directed = TRUE)
    X  = matrix(0, nTips, nEdges)

    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode="out", output="epath")$epath
    ## there must be always a path.
    stopifnot( all( lapply( root2tip, length ) > 0) ) 
    ## since it is ultrametric.
    Tval  = sum(tree$edge.length[root2tip[[1]]]) 

    if(type == "orgX"){
        for(i in 1:nTips){
            lvec    = c(0, tree$edge.length[root2tip[[i]]])
            timeVec = Tval - cumsum(lvec)
            timeVec = timeVec[1:length(timeVec)-1]
            X[i, root2tip[[i]] ] = 1 - exp(-alpha*timeVec)
        }
    }else if(type == "apprX"){
        for(i in 1:nTips){
            lvec    = c(0, tree$edge.length[root2tip[[i]]])
            timeVec = Tval - cumsum(lvec)
            timeVec = timeVec[1:length(timeVec)-1]
            X[i, root2tip[[i]] ] = timeVec
        }
    }else if(type == "simpX"){
        for(i in 1:nTips){
            X[i, root2tip[[i]] ] = 1
        }
    }else
        stop("Undefined design matrix type")
    return(X)
}

select_best_solution <- function(tree, Y, sol.path, opt){

    nSols = get_num_solutions(sol.path)
    stopifnot( nSols > 0 )

    all.shifts = numeric()
    prev.shift.configuration = NA
    min.score = Inf   

    candid.idx <- 1
    shift.configuration.list <- list()
    for (idx in 1:nSols) {

        shift.configuration = get_configuration_in_sol_path(sol.path, idx, Y)
        shift.configuration = correct_unidentifiability(tree, shift.configuration, opt)

        if ( length(shift.configuration) > opt$max.nShifts ){break}
        if ( setequal(shift.configuration, prev.shift.configuration) ){next}
        prev.shift.configuration = shift.configuration

        ## sorting shifts based on their age in the solution path
        all.shifts  = c(all.shifts, shift.configuration)
        freq.shifts = rep(0, length(shift.configuration))
        count = 1
        for( s in shift.configuration){
            freq.shifts[[count]] = length( which(all.shifts == s) )
            count = count + 1
        }
        names(shift.configuration)  <- freq.shifts
        shift.configuration <- shift.configuration[order(names(shift.configuration), decreasing=TRUE)]
        shift.configuration.list[[candid.idx]] <- shift.configuration
        candid.idx <- candid.idx + 1
    }

    search_ith_config <- function(sc){
        res <- do_backward_correction(tree, Y, sc, opt)
        return(res)
    }

    if(opt$parallel.computing){
        all.res <- mclapply(rev(shift.configuration.list), 
                            FUN=search_ith_config, 
                            mc.cores=opt$nCores)
        for (i in length(all.res):1 ){
            res <- all.res[[i]] 
            if (min.score > res$score){
                min.score = res$score
                best.shift.configuration = res$shift.configuration
            }
        }
    }else{
        for (i in 1:length(shift.configuration.list) ){
            res <- search_ith_config( shift.configuration.list[[i]] )
            if (min.score > res$score){
                min.score = res$score
                best.shift.configuration = res$shift.configuration
            }
        }
    }

    return ( list(score=min.score, shift.configuration=best.shift.configuration) )
}

do_backward_correction <- function(tree, Y, shift.configuration, opt){

    org.score = cmp_model_score(tree, Y, shift.configuration, opt)

    if( length(shift.configuration) < 3 ) { 
        return(list(score=org.score, shift.configuration=shift.configuration)) 
    }  

    #nShifts = length(shift.configuration)
    #removal.candids = shift.configuration[1:(nShifts-1)]
    #for( sp in removal.candids )
    for(sp in shift.configuration)
    {
        new.configuration = setdiff(shift.configuration, sp)
        new.score         = cmp_model_score(tree, Y, new.configuration, opt)      
        if ( new.score < org.score ){
            shift.configuration = new.configuration
            org.score           = new.score
        }
    }

    return(list(score=org.score, shift.configuration=shift.configuration))
}


#
#' Computes the information criterion score for a given configuration
#'
#'@param tree ultrametric tree of class phylo, with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param shift.configuration shift positions, i.e. vector of indices of the edges where the shifts occur.
#'@param criterion an information criterion (see Details).
#'@param root.model an ancestral state model at the root.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum length of external branches, corresponding to a half life greater or equal to the minimum external branch length.
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param fit.OU.model logical. If TRUE, it returns an object of class l1ou with all the parameters estimated.
#'@param l1ou.options if provided, all the default values will be ignored. 
#'
#'@return Information criterion value of the given shift configuration.
#'
#'@details
#'AICc gives the usual small-sample size modification of AIC. 
#'BIC gives the usual Bayesian information criterion, here penalizing each shift as 2 parameters. 
#'mBIC is the modified BIC proposed by Ho and Ané (2014).
#'pBIC is the phylogenetic BIC for shifts proposed by Khabbazian et al.
#'pBICess is a version of pBIC where the determinant term is replaced by a sum of the log of effective sample sizes (ESS), similar to the ESS proposed by Ané (2008). 
#' 
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' lizard <- adjust_data(lizard.tree, lizard.traits[,1])
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y)
#' configuration_ic(lizard$tree, eModel$Y, eModel$shift.configuration, criterion="pBIC")
#'
#' ### building l1ou object out of the second best score 
#' eModel2 = configuration_ic(eModel$tree, eModel$Y, eModel$profile$configurations[[2]], 
#'                           fit.OU.model=TRUE, l1ou.options=eModel$l1ou.options)
#' plot(eModel2)
#'
#'@seealso \code{\link{estimate_shift_configuration}} \code{\link{adjust_data}}
#'
#'@references
#'Cécile Ané, 2008. "Analysis of comparative data with hierarchical autocorrelation". Annals of Applied Statistics 2(3):1078-1102.
#'
#'Ho, L. S. T. and Ané, C. 2014.  "Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. 5(11):1133-1146.
#'
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
configuration_ic <- function(tree, Y, shift.configuration, 
                     criterion    = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
                     root.model   = c("OUfixedRoot", "OUrandomRoot"),
                     alpha.starting.value = NA,
                     alpha.upper  = alpha_upper_bound(tree), 
                     alpha.lower  = NA,
                     fit.OU.model = FALSE, 
                     l1ou.options = NA
                   ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the input phylogenetic tree is not in postorder. Use adjust_data function.")

    Y  = as.matrix(Y)
    if(!identical(rownames(Y), tree$tip.label)) stop("rownames of Y and tree$tip.label are not identical.")


    opt = list()
    if(!all(is.na(l1ou.options))){
        opt = l1ou.options
    }else{
        opt$criterion            <- match.arg(criterion)
        opt$root.model           <- match.arg(root.model)
        opt$alpha.starting.value <- alpha.starting.value
        opt$alpha.upper.bound    <- alpha.upper
        opt$alpha.lower.bound    <- alpha.lower
        opt$Z                    <- generate_design_matrix(tree, "simpX")
        opt$multivariate.missing <- FALSE
        opt$use.saved.scores     <- FALSE
    }

    s.c = correct_unidentifiability(tree, shift.configuration, opt)
    if( length(s.c) != length(shift.configuration) )
        stop(paste0("the input shift configuration is not a parsimony configuration. 
                    For instance,\n", s.c, "\n is an alternative configuration with fewer shifts."))

    if( fit.OU.model ){
        ##TODO: convergent evolution.
        eModel = fit_OU_model(tree, Y, shift.configuration, opt)
        return(eModel)
    }
    score = cmp_model_score(tree, Y, shift.configuration, opt)
    return(score)
}


#
#' Fits an OU model based on a given configuration
#'
#'@param tree ultrametric tree of class phylo, with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param shift.configuration shift positions, i.e. vector of indices of the edges where the shifts occur.
#'@param criterion an information criterion (see Details).
#'@param root.model model for the ancestral state at the root.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum length of external branches, corresponding to a half life greater or equal to the minimum external branch length.
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param l1ou.options if provided, all the default values will be ignored. 
#'
#'@return an object of class l1ou similar to \code{\link{estimate_shift_configuration}}.
#'
#'@details
#'AICc gives the usual small-sample size modification of AIC. 
#'BIC gives the usual Bayesian information criterion, here penalizing each shift as 2 parameters. 
#'mBIC is the modified BIC proposed by Ho and Ané (2014).
#'pBIC is the phylogenetic BIC for shifts proposed by Khabbazian et al.
#'pBICess is a version of pBIC where the determinant term is replaced by a sum of the log of effective sample sizes (ESS), similar to the ESS proposed by Ané (2008). 
#' 
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' lizard <- adjust_data(lizard.tree, lizard.traits[,1])
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y)
#'
#' ### building l1ou object out of the second best score 
#' eModel2 = fit_OU(eModel$tree, eModel$Y, eModel$profile$configurations[[2]], 
#'                           l1ou.options=eModel$l1ou.options)
#' plot(eModel2)
#' 
#' ### hypothesis testing
#'
#' data("lizard.traits", "lizard.tree")
#' Y <- lizard.traits[,1:1]
#' tr  <- lizard.tree
#' 
#' tr <- multi2di(tr)
#' tr <- reorder(tr, "postorder")
#' 
#' ### visualizing the tree with the edge indeces 
#' plot(tr)
#' edgelabels()
#' 
#' ## place the shift position based on the hypothesis
#' shift.config <- c(116, 77)
#' 
#' hModel <- fit_OU(tr, Y, shift.config, criterion="AICc")
#' plot(hModel)
#' print(hModel)
#'
#'@seealso \code{\link{estimate_shift_configuration}} \code{\link{adjust_data}}
#'
#'@references
#'Cécile Ané, 2008. "Analysis of comparative data with hierarchical autocorrelation". Annals of Applied Statistics 2(3):1078-1102.
#'
#'Ho, L. S. T. and Ané, C. 2014.  "Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. 5(11):1133-1146.
#'
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
fit_OU <- function(tree, Y, shift.configuration, 
                     criterion    = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
                     root.model   = c("OUfixedRoot", "OUrandomRoot"),
                     cr.regimes   = NULL,
                     alpha.starting.value = NA,
                     alpha.upper  = alpha_upper_bound(tree), 
                     alpha.lower  = NA,
                     l1ou.options = NA
                   ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the input phylogenetic tree is not in postorder. Use adjust_data function.")

    Y  = as.matrix(Y)
    if(!identical(rownames(Y), tree$tip.label)) stop("rownames of Y and tree$tip.label are not identical.")

    opt = list()
    if(!all(is.na(l1ou.options))){
        opt = l1ou.options
    }else{
        opt$criterion            <- match.arg(criterion)
        opt$root.model           <- match.arg(root.model)
        opt$alpha.starting.value <- alpha.starting.value
        opt$alpha.upper.bound    <- alpha.upper
        opt$alpha.lower.bound    <- alpha.lower
        opt$Z                    <- generate_design_matrix(tree, "simpX")
        opt$multivariate.missing <- FALSE
        opt$use.saved.scores     <- FALSE
    }

    s.c = correct_unidentifiability(tree, shift.configuration, opt)
    if( length(s.c) != length(shift.configuration) )
        stop(paste0("the input shift configuration is not parsimonious. For instance, shifts on these edges:\n",
                    s.c, "provides an alternative equivalent configuration with fewer shifts."))

     eModel = fit_OU_model(tree, Y, shift.configuration, opt)
     if(!is.null(cr.regimes) ){

        if( !( 0 %in% unlist(cr.regimes) ) ){
            stop("background/intercept is not included in the regimes! Represent the background by \"0\".")
        }
        if( !(identical(sort(c(0,shift.configuration)), sort(unlist(cr.regimes)) ) ) ){
            stop("convergent regimes do not match with the shift positions.")
        }
        cr.score <- cmp_model_score_CR(tree, Y, shift.configuration, regimes=cr.regimes, 
                           criterion=opt$criterion, alpha=eModel$alpha)
        eModel$cr.score <- cr.score
        for(idx in 1:length(cr.regimes)){
            if( 0 %in% cr.regimes[[idx]] ){
                names(eModel$shift.configuration)[which(eModel$shift.configuration %in% cr.regimes[[idx]])] <- 0
            }else{
                names(eModel$shift.configuration)[which(eModel$shift.configuration %in% cr.regimes[[idx]])] <- idx
            }
        }
     }
     return(eModel)
}

fit_OU_model <- function(tree, Y, shift.configuration, opt){

    Y      = as.matrix(Y)
    nEdges = Nedge(tree)
    nTips  = length(tree$tip.label)

    resi = mu = optima = matrix(data=NA, nrow=nTips, ncol=ncol(Y))
    shift.values = numeric() # possibly different # of shifts for different traits if missing values
    # edge.optima = matrix(NA, nEdges, ncol(Y))
    intercept = alpha = sigma2 = rep(NA, ncol(Y))
    logLik = rep(NA, ncol(Y))

    for(i in 1:ncol(Y)){
        s.c <- c()
        if(!is.null(opt$tree.list)){
            tr <- opt$tree.list[[i]]
            y  <- as.matrix(Y[!is.na(Y[,i]), i])
            if(length(shift.configuration > 0))
                augmented.s.c <- tr$old.order[shift.configuration] # index of shift edges in pruned tree, see gen_tree_array
            else
                augmented.s.c <- c()
            for(s in shift.configuration){
                n.s <- tr$old.order[[s]]
                if(!is.na(n.s)){ s.c <- c(s.c, n.s) }
            }
        } else{
            tr  <- tree
            y   <- as.matrix(Y[,i])
            s.c <- shift.configuration
            augmented.s.c <- shift.configuration
        }

        nShifts = length(s.c)

        fit <- my_phylolm_interface(tr, y, s.c, opt)
        if ( all(is.na(fit)) ){
            stop("model score is NA in fit_OU_model function! This should not happen.")
        }

        alpha[i]  <- fit$optpar
        sigma2[i]   <- fit$sigma2
        logLik[i] <- fit$logLik

        ## E[Y] and residuals: Y-EY
        mu[!is.na(Y[,i]), i]   = fit$fitted.values
        resi[!is.na(Y[,i]), i] = fit$residuals
        intercept[i] = fit$coefficients[[1]] # = y0 * e^-T + theta0_root * (1-e^-T), assumes ultrametric tree

        ## Now we have the alpha hat and we can form the true design matrix
        if( nShifts > 0 ){
            scale.values <- apply( generate_design_matrix(tr, type="orgX", alpha=alpha[i])[,s.c], 2, max)^-1
            fit$coefficients[2:(nShifts+1)] <- scale.values * fit$coefficients[2:(nShifts+1)]
        }

        if( length(shift.configuration) > 0 ){
            s.v = rep(NA, length(shift.configuration))
            s.v[!is.na(augmented.s.c)] = fit$coefficients[2:(nShifts+1)]
            shift.values = cbind(shift.values, s.v)
        }

        optima.tmp = rep(fit$coefficients[[1]], nTips)  # optima at the tips for one trait
        if( length(shift.configuration) > 0 )
            for(ish in 1:length(shift.configuration) ){ # i=index of Y column. ish=index of shift
                s <- shift.configuration[[ish]]         # edge number for shift 'ish' in full tree
                if(is.na(augmented.s.c[[ish]]))         # shift invisible in pruned tree
                    next;
                s.v <- fit$coefficients[which(s.c==augmented.s.c[[ish]])+1] # shift value
                optima.tmp = optima.tmp + opt$Z[,s] * s.v # requires Z to have 0/1 values. bug otherwise.
            }

        optima[,i] <- optima.tmp
    }

    rownames(optima) <- tree$tip.label

    ##NOTE: it reads the score from the database 
    ## and do not recompute the score. So it doesn't have any overhead.
    score = cmp_model_score (tree, Y, shift.configuration, opt) 

    model = list(
                 Y                   = Y, 
                 tree                = tree,
                 shift.configuration = shift.configuration, 
                 shift.values        = shift.values,
                 nShifts             = length(shift.configuration), 
                 optima              = optima, 
                 alpha               = alpha, 
                 sigma2              = sigma2, 
                 intercept           = intercept, 
                 mu                  = mu, 
                 residuals           = resi,
                 score               = score,
                 logLik              = logLik,
                 l1ou.options        = opt) 

    class(model) <- "l1ou"
    return( model )
}

cmp_model_score <-function(tree, Y, shift.configuration, opt){

    ##TODO: optimize it.
    shift.configuration <- correct_unidentifiability(tree, shift.configuration, opt)

    if(opt$use.saved.scores){ ##if it's been already computed
        score <- get_configuration_score_from_list(shift.configuration)
        if(!is.na(score)){ return(score) }
    }

    Y  <- as.matrix(Y)
    ic <- opt$criterion

    res <- NA
    if( ic == "AIC"){
        stop("undefined")
    } else if( ic == "BIC"){
        res <- cmp_BIC(tree, Y, shift.configuration, opt)
    } else if( ic == "AICc"){
        res <- cmp_AICc(tree, Y, shift.configuration, opt)
    } else if( ic == "mBIC"){
        res <- cmp_mBIC(tree, Y, shift.configuration, opt) 
    } else if( ic == "pBICess"){
        res <- cmp_pBICess(tree, Y, shift.configuration, opt) 
    } else if(ic == "pBIC"){
        res <- cmp_pBIC(tree, Y, shift.configuration, opt) 
    } 

    if(all(is.na(res))){return(Inf)}

    score <- res$score
    if(opt$use.saved.scores){
        add_configuration_score_to_list(shift.configuration, score,
             paste0(c(res$sigma2/(2*res$alpha),res$logLik),collapse=" "))
    }
    return(score)
}

get_data <- function(tree, Y, shift.configuration, opt, idx){

    if(!is.null(opt$tree.list)){
        tr    <- opt$tree.list[[idx]]
        y.ava <- !is.na(Y[,idx])
        y     <- as.matrix(Y[y.ava, idx])
        s.c   <- c()
        for(s in shift.configuration){
            n.s <- tr$old.order[[s]]
            if(!is.na(n.s)){
                s.c <- c(s.c, n.s)
            }
        }
        stopifnot(length(tr$tip.label)==nrow(y))
    }else{
        tr  <- tree
        y   <- as.matrix(Y[, idx])
        s.c <- shift.configuration
    }
    result     <- list()
    result$tr  <- tr
    result$y   <- y
    result$s.c <- s.c
    return(result)
}

cmp_BIC <- function(tree, Y, shift.configuration, opt){

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    df.1  <- log(nTips)*(nShifts)
    score <- df.1
    alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        df.2 <- log(nrow(y))*(length(s.c)+ 3)
        fit  <- my_phylolm_interface(tr, y, s.c, opt)
        if ( all(is.na(fit)) ){ return(NA) } 

        score <- score  -2*fit$logLik + df.2

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}


cmp_AICc <- function(tree, Y, shift.configuration, opt){

    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    p   <- nShifts + (nShifts + 2)*nVariables
    N   <- nTips*nVariables
    d.f <- 2*p + (2*p*(p+1))/(N-p-1) 
    if( p > N-2 )
        return(NA)

    score <- d.f
    alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit <- my_phylolm_interface(tr, y, s.c, opt)
        if ( all(is.na(fit)) ){ return(NA) } 
        score <- score  -2*fit$logLik 

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }

    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

cmp_mBIC <- function(tree, Y, shift.configuration, opt){

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    res =  cmp_mBIC_df(tree, shift.configuration, opt)  
    df.1 = res$df.1
    df.2 = res$df.2

    score <- df.1

    alpha <- sigma2 <- logLik <- rep(0, nVariables)
    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        if( nVariables > 1){
            res  = cmp_mBIC_df(tr, s.c, opt)  
            df.2 = res$df.2
        }

        fit <- my_phylolm_interface(tr, y, s.c, opt)
        if ( all(is.na(fit)) ){ return(NA) } 

        score <- score  -2*fit$logLik + df.2

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

cmp_mBIC_df <- function(tree, shift.configuration, opt){

    if(length(shift.configuration)>0){
        shift.configuration <- sort(shift.configuration)
    }
    nTips               <- length(tree$tip.label)
    nShifts             <- length(shift.configuration)

    df.1 <- 0 
    ## pen for the alpha sigma2 and intercept
    df.2 <- 3*log(nTips)

    if(nShifts > 0 ){
        ## pen for shift configuration
        df.1 = (2*nShifts - 1) *log(nTips)
        ## pen for alpha sigma2 and intercept
        df.2 = 3*log(nTips)

        all.covered.tips = numeric()
        for(eIdx in shift.configuration){
            covered.tips = which( opt$Z[,eIdx] > 0 )
            nUniqueTips  = length( setdiff(covered.tips, all.covered.tips) )
            all.covered.tips = union(covered.tips, all.covered.tips)
            ## this must not happen if the input is an 
            ## identifiable configuration (parsimonious) and the tree is in post order.
            stopifnot( nUniqueTips > 0)
            df.2 = df.2 + log(nUniqueTips) 
        }
        nUniqueTips = length( setdiff(1:nTips, all.covered.tips) )
        df.2 = df.2 + log(nUniqueTips) 
    } 

    return( list(df.1=df.1, df.2=df.2) )
}

cmp_pBICess <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1  = 2*(nShifts)*log(nEdges-1)
    score = df.1

    alpha = sigma2  = logLik = numeric()

    for(i in 1:ncol(Y)){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit  = my_phylolm_interface(tr, y, s.c, opt)
        if( all(is.na(fit)) ){
           return(NA)
        }
        ess  = effective.sample.size(tr, edges=s.c, model="OUfixedRoot", 
                 parameters=list(alpha=fit$optpar), FALSE, FALSE)

        df.2  = 3*log(nrow(y)+1) + sum(log(ess+1))
        score = score  -2*fit$logLik + df.2 

        alpha  = c(alpha, fit$optpar)
        sigma2 = c(sigma2, fit$sigma2)
        logLik = c(logLik, fit$logLik)
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}


cmp_pBIC <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1    = 2*(nShifts)*log(nEdges-1)
    score   = df.1
    alpha   = sigma2  = logLik = rep(0, ncol(Y))

    for(i in 1:ncol(Y)){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit   = my_phylolm_interface(tr, y, s.c, opt)
        if( all(is.na(fit)) ){
           return(NA)
        } 
        varY  = c(var(y))
        ld    = as.numeric(determinant(fit$vcov * (fit$n - fit$d)/(varY*fit$n), log=T)$modulus)
        df.2  = 2*log(nrow(y)) - ld
        score = score  -2*fit$logLik + df.2 

        alpha [[i]] = fit$optpar
        sigma2[[i]] = fit$sigma2
        logLik[[i]] = fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

my_phylolm_interface <- function(tree, Y, shift.configuration, opt, recmp.preds=FALSE, alpha=0){

    if(recmp.preds){
        Z <- generate_design_matrix(tree, type="orgX", alpha=alpha)
    }else{
        if(opt$multivariate.missing){
            if(is.null(tree$Z))
                stop("internal error: in multivariate.missing mode but tree$Z is null.")
            Z <- tree$Z
        }else{
            Z <- opt$Z 
        }
    }


    preds = cbind(1, Z[ ,shift.configuration])
    prev.val <-options()$warn 
    options(warn = -1)
    if( is.na(opt$alpha.lower.bound) & is.na(opt$alpha.starting.value) ){
        fit    <-  try( phylolm(Y~preds-1, phy=tree, model=opt$root.model,
                 upper.bound  = opt$alpha.upper.bound), silent = opt$quietly)
    } else {
        l = ifelse(is.na(opt$alpha.lower.bound), 0, opt$alpha.lower.bound)
        u = opt$alpha.upper.bound
        s = ifelse(is.na(opt$alpha.starting.value), max(0.5, l), opt$alpha.starting.value)

        fit    <-  try( phylolm(Y~preds-1, phy = tree, 
                                model          = opt$root.model,
                                starting.value = s, 
                                lower.bound    = l, 
                                upper.bound    = u), silent = opt$quietly)
    }
    options(warn = prev.val )

    if(class(fit) == "try-error"){ 
        if(!opt$quietly){
            warning( paste0( "phylolm returned error with a shift configuration
                            of size ", length(shift.configuration), ". You may
                            want to change alpha.upper/alpha.lower!") )
        }
        return(NA)
    }
    return(fit)
}


run_grplasso  <- function (grpX, grpY, nVariables, grpIdx, opt){
    delta  = opt$grp.delta
    seq.ub = opt$grp.seq.ub
    max.nTries = 7
    suppressMessages(
      lmbdMax  <-  1.2 * lambdamax(grpX, grpY, model = LinReg(), index = grpIdx, rescale = FALSE) + 1
    )

    base.seq = seq(0, seq.ub, delta)
    for (itrTmp in 1:max.nTries) {
        lmbd = lmbdMax * (0.5^base.seq)

        capture.output(
          sol  <-  grplasso(grpX, y = grpY, rescale = FALSE, center = FALSE, 
                            lambda = lmbd, model = LinReg(), index = grpIdx, 
                            control = grpl.control(tol = 0.01))
          )

        #df.vec = apply(sol$coefficients, 2, function(x) length(which(abs(x) > 0))/nVariables)
        df.vec = apply( sol$coefficients, 2, 
           function(x) length(which(rowSums(matrix(ifelse(abs(x[!is.na(grpIdx)])>0,1,0),
           ncol=nVariables))>nVariables-1)) )

        df.missing = setdiff(0:(opt$max.nShifts+1), df.vec)

        cutExtra = TRUE
        if (length(df.missing) > 0) {
            base.seq.tmp =numeric()
            idx = prev.idx = 1
            for (mdf in df.missing) {

                if( length( which(df.vec > mdf) ) > 0 ){ ##that means at the begining or in the middle of the sequence
                    idx = min(which(df.vec > mdf))
                    if( idx == 1){ ## that means we should add to the begining
                        ##TODO: I don't like it
                        lmbdMax    = lmbdMax + 2
                        next
                    }
                    if( prev.idx == idx){ ## not to repeat the same sequence.
                        next
                    }
                    lower = base.seq[[idx - 1]] + delta/4
                    upper = base.seq[[idx]]
                    base.seq.tmp=c( base.seq.tmp, base.seq[prev.idx:idx], seq(lower, upper, delta/4))
                } else{ ## that means we should add to the end
                    if( length(base.seq.tmp) > 0){
                        lower = base.seq.tmp[[length(base.seq.tmp)]]+delta
                    }else{lower = 0}
                    upper = max(lower, max(base.seq)) + 1
                    base.seq.tmp = c(base.seq.tmp, seq(lower, upper, delta))
                    cutExtra  = FALSE
                }
                prev.idx = idx
            }

            if( cutExtra ){
                indices = which(df.vec > (opt$max.nShifts+4))
                if( length(indices) > 0 ){
                    upper.idx = min(indices)
                }else{
                    upper.idx = length(base.seq)
                }

                if( idx < upper.idx){
                    base.seq.tmp = c(base.seq.tmp, base.seq[idx:upper.idx] ) 
                }
            }
            delta = delta/4
            base.seq = unique(sort(base.seq.tmp))
        }else { 
            indices = which(df.vec > (opt$max.nShifts + 4))
            if (length(indices) > 0) {
                upper.idx = min(indices)
                base.seq  = base.seq[1:upper.idx]
            }
            break 
        }
    }
    lmbd = lmbdMax * (0.5^base.seq)
     capture.output(
             sol  <-  grplasso(grpX, y = grpY, rescale = FALSE, center = FALSE, 
                               lambda = lmbd, model = LinReg(), index = grpIdx, 
                               control = grpl.control(tol = 1e-6) )
             )
    df.missing = setdiff(0:opt$max.nShifts, df.vec)
    for (dfm in df.missing) {
        warning(paste0("There are no solutions with ", dfm, " number of shifts  
                in the solution path of grplasso. You may want to change grp.delta and grp.seq"))
    }
    return(sol);
}


