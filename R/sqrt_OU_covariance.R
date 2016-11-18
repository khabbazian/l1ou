#
#' (inverse) square root of the phylogenetic covariance
#'
#' Computes an inverse square root and square root of the phylogenetic covariance matrix,
#' under the Brownian motion (BM) or the Ornstein-Uhlenbeck (OU) model.
#' The algorithm traverses the tree only once, hence the algorithm is very fast
#' and can be applied to very big trees.
#'
#'@param tree tree of class phylo with branch lengths. If alpha>0, i.e. under the OU model, the tree has to be ultrametric.
#'@param alpha adaptation rate for the OU model. The default is 0, which corresponds to the BM mode with a fixed ancestral state at the root.
#'@param root.model ancestral state model at the root.
#'@param check.order logical. If TRUE, the order will be checked to be in postorder traversal.
#'@param check.ultrametric logical. If TRUE, the tree will be checked to ultrametric.
#'@param normalize.tree.height logical. If TRUE, it class normalize_tree function after transf.branch.lengths.
#'
#'@return 
#' \item{sqrtInvSigma}{inverse square root of the phylogenetic covariance matrix.}
#' \item{sqrtSigma}{square root of the phylogenetic covariance matrix.}
#'
#'@examples
#'
#' data(lizard.tree)
#' res <- sqrt_OU_covariance(lizard.tree) # alpha not provided: so BM model.
#' Sigma <- vcv(lizard.tree)
#' dimnames(Sigma) <- NULL
#' all.equal(res$sqrtSigma %*% t(res$sqrtSigma), Sigma) # TRUE
#' all.equal(res$sqrtInvSigma %*% t(res$sqrtInvSigma), solve(Sigma)) # TRUE
#' 
#' 
#' ##Here's the example from "Eric A. Stone. 2011." (See references)
#'
#' tr <-  read.tree(text="((((Homo:.21,Pongo:.21):.28,Macaca:.49):.13,Ateles:.62):.38,Galago:1);") 
#' RE <- sqrt_OU_covariance(tr) 
#' B <- round( RE$sqrtSigma, digits=3)
#' D <- round( RE$sqrtInvSigma, digits=3)
#' print(B)
#' print(D)
#' 
#' 
#' ##Here is the examples on how to get the contrasts using sqrt_OU_covariance
#' data(lizard.tree, lizard.traits)
#' lizard <- adjust_data(lizard.tree, lizard.traits)
#' eModel <- estimate_shift_configuration(lizard$tree, lizard$Y)
#' theta <- eModel$intercept + l1ou:::convert_shifts2regions(eModel$tree,
#'                              eModel$shift.configuration, eModel$shift.values)
#' REf <- sqrt_OU_covariance(eModel$tree, alpha=eModel$alpha,
#'                                          root.model = "OUfixedRoot",normalize.tree.height=T,
#'                                          check.order=F, check.ultrametric=F)
#' REr <- sqrt_OU_covariance(eModel$tree, alpha=eModel$alpha,
#'                                          root.model = "OUrandomRoot",normalize.tree.height=T,
#'                                          check.order=F, check.ultrametric=F)
#' # `covInverseSqrt` represents the transpose of square root of  the inverse matrix of covariance.
#' # `covSqrt` represents the square root of the covariance matrix.
#'
#'  covInverseSqrtf  <- t(REf$sqrtInvSigma)
#'  covSqrtf   <- REf$sqrtSigma
#'  covInverseSqrtr  <- t(REr$sqrtInvSigma)
#'  covSqrtr   <- REr$sqrtSigma
#'  tij=vcv(eModel$tree)  # the time spent on each edge
#'  treeheight=max(tij)   
#'  dij=2*(treeheight-tij) 
#'  vrandom=1/(2*eModel$alpha)*exp(-eModel$alpha*dij) 
#'  vfix= 1/(2*eModel$alpha)*exp(-eModel$alpha*dij)*(1-exp(-2*eModel$alpha*tij))
#'  Ind_random=covInverseSqrtr%*%vrandom%*%t(covInverseSqrtr)
#'  Ind_fix=covInverseSqrtf%*%vfix%*%t(covInverseSqrtf)
#'  all.equal(Ind_random,diag(100))
#'  all.equal(Ind_fix,diag(100))
#'  Y  <- rTraitCont(eModel$tree, "OU", theta=theta, 
#'                                      alpha=eModel$alpha, 
#'                                      sigma=eModel$sigma, root.value=eModel$intercept)
#'  contrast    <-  covInverseSqrt%*%(Y - eModel$mu)
#'
#'
#'@references
#' Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#' Eric A. Stone. 2011. "Why the phylogenetic regression appears robust to tree misspecification". Systematic Biology, 60(3):245-260.
#'
#'@export
sqrt_OU_covariance <- function(tree, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot"), 
                               check.order=TRUE, check.ultrametric=TRUE, normalize.tree.height=FALSE){
    if( ! is.binary.tree(tree) ){
        tree         <- multi2di(tree, random=FALSE)
        check.order  <- TRUE 
    }
    root.model <- match.arg(root.model) 
    ##NOTE: the function assumes reordering does not change the order of the 
    ##nodes and it just change the order of edges, so that column i in each 
    ##matrix still corresponds to internal node 
    ##NOTE:  in case the tree is not binary; the order will change and it is no longer postorder. 
    if( check.order ){
        tree <- reorder(tree, "post")
    }

    if ( alpha > 0){
        ##NOTE: this step requires that the tree be ultrametric tree. 
        ##NOTE: If the tree is not ultrametric, the function returns a wrong result with no warning
        if(check.ultrametric){
            if(!is.ultrametric(tree)){
                stop("alpha>0, the tree has to be ultrametric") 
            }
        }
        tre <- transf.branch.lengths(tree, model=root.model, parameters=list(alpha=alpha))$tree
	if(normalize.tree.height){
		tre <- normalize_tree(tre)
	}
    }else{
        tre <- tree
        if( root.model == "OUrandomRoot"){
            warning("alpha=0, BM model, the ancestral state model is changed to the OUfixedRoot")
        }
    }

    my.edge.list <- cbind(tre$edge-1, tre$edge.length) 
    tre$root.edge <- ifelse(is.null(tre$root.edge), 0, tre$root.edge)
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tre$tip.label), tre$root.edge)
    tij=vcv(tre)  # the time spent on each edge
    treeheight=max(tij) 
    if (root.model == "OUrandomRoot"){
      result$sqrtSigma=result$sqrtSigma/sqrt(2*alpha)
      result$sqrtInvSigma=result$sqrtInvSigma*sqrt(2*alpha)
    } else if (root.model == "OUfixedRoot"){
      result$sqrtSigma=result$sqrtSigma/(sqrt((2*alpha)/(1-exp(-2*alpha*treeheight))))
      result$sqrtInvSigma=result$sqrtInvSigma*(sqrt((2*alpha)/(1-exp(-2*alpha*treeheight))))
                             
    }
    return(result)
}

