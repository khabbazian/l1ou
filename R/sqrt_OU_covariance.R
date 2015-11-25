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
#'@param check.ultramteric logical. If TRUE, the tree will be checked to ultrametric.
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
#'@references
#' Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané. "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models". In review.
#'
#' Eric A. Stone. 2011. "Why the phylogenetic regression appears robust to tree misspecification". Systematic Biology, 60(3):245-260.
#'
#'@export
sqrt_OU_covariance <- function(tree, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot"), 
                               check.order=TRUE, check.ultramteric=TRUE){

    tree       <- multi2di(tree, random=FALSE)
    root.model <- match.arg(root.model) 
    ##NOTE: the function assumes reordering does not change the order of the 
    ##nodes and it just change the order of edges, so that column i in each 
    ##matrix still corresponds to internal node 
    ##NOTE:  in case the tree is not binary; the order will change and it is no longer postorder. 
    #if( check.order ){
        tree <- reorder(tree, "post")
    #}

    if ( alpha > 0){
        ##NOTE: this step requires that the tree be ultrametric tree. 
        ##NOTE: If the tree is not ultrametric, the function returns a wrong result with no warning
        if(check.ultramteric){
            if(!is.ultrametric(tree)){
                stop("alpha>0, the tree has to be ultrametric") 
            }
        }
        tre <- transf.branch.lengths(tree, model=root.model, parameters=list(alpha=alpha))$tree
        tre <- normalize_tree(tre, check.ultramteric=FALSE)
    }else{
        tre <- tree
        if( root.model == "OUrandomRoot"){
            warning("alpha=0, BM model, the ancestral state model is changed to the OUfixedRoot")
        }
    }

    my.edge.list <- cbind(tre$edge-1, tre$edge.length) 
    tre$root.edge <- ifelse(is.null(tre$root.edge), 0, tre$root.edge)
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tree$tip.label), tre$root.edge)
    return(result)
}

