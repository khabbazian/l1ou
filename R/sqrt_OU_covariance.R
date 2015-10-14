
#
#' (inverse) square root of the phylogenetic covariance
#'
#' Computes an inverse square root and square root of the phylogenetic covariance matrix,
#' under the Brownian motion (BM) or the Ornstein-Uhlenbeck (OU) model.
#' The algorithm traverses the tree only once, hence the algorithm is very fast
#' and can be applied to very big trees.
#'
#'@param tree ultrametric tree of class phylo with branch lengths.
#'@param alpha adaptation rate for the OU model. The default is 0, which corresponds to the BM mode with a fixed ancestral state at the root.
#'@param root.model ancestral state model at the root.
#'
#'@return 
#' \item{sqrtInvSigma}{inverse square root of the phylogenetic covariance matrix.}
#' \item{sqrtSigma}{square root of the phylogenetic covariance matrix.}
#'
#'@examples
#'
#' library(l1ou)
#' data(lizard.tree)
#' res <- sqrt_OU_covariance(lizard.tree) # alpha not provided: so BM model.
#' Sigma <- vcv(lizard.tree)
#' dimnames(Sigma) <- NULL
#' all.equal(res$sqrtSigma %*% t(res$sqrtSigma) , Sigma) # TRUE
#' all.equal(res$sqrtInvSigma %*% t(res$sqrtInvSigma) , solve(Sigma)) # TRUE
#'
#'@references
#' Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané. Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models 
#'
#' Eric A. Stone. 2011. "Why the phylogenetic regression appears robust to tree misspecification". Systematic Biology, 60(3):245-260.
#'
#'@export
sqrt_OU_covariance <- function(tree, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot") ){

    tree         <- multi2di(tree, random=FALSE)
    root.model <- match.arg(root.model) 
    tree         <- reorder(tree, "prun")

    if ( alpha > 0){
        tre <- transf.branch.lengths(tree, model=root.model, parameters=list(alpha=alpha))$tree
    } else{
        tre <- tree
        if( root.model == "OUrandomRoot"){
            warning("when alpha is zero the model should be OUfixedRoot, so I switched to the OUfixedRoot")
        }
    }


    my.edge.list <- cbind(tre$edge-1, tre$edge.length) 
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tree$tip.label))
    return(result)
}

