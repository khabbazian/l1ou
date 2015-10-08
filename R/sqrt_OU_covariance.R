
#
#' computes the negative square root and square root of the phylogeny covariance matrix. 
#'
#'@param tree an ultrametric phylogenetic tree of class phylo with branch lengths.
#'@param alpha the adaptation rate for the OU model.
#'@param root.model an ancestral state model at the root.
#'
#'@return 
#' \item{sqrtSigma}{square root of the phylogeny covariance matrix.}
#' \item{sqrtInvSigma}{inverse square root of the phylogeny covariance matrix.}
#'
#'@examples
#'
#' library("l1ou")
#' data("lizard.tree")
#' res <- sqrt_OU_covariance(lizard.tree)
#' Sigma <- vcv(lizard.tree)
#' dimnames(Sigma) <- NULL
#' all.equal(res$sqrtSigma %*% t(res$sqrtSigma) , Sigma) # TRUE
#' all.equal(res$sqrtInvSigma %*% t(res$sqrtInvSigma) , solve(Sigma)) # TRUE
#'
#'@export
#'
sqrt_OU_covariance <- function(tree, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot") ){

    tree         <- multi2di(tree, random=FALSE);
    root.model <- match.arg(root.model); 
    tree         <- reorder(tree, "prun");

    if ( alpha > 0){
        tre <- transf.branch.lengths(tree, model=root.model, parameters=list(alpha=alpha))$tree;
    } else{
        tre <- tree;
        if( root.model == "OUrandomRoot"){
            warning("when alpha is zero the model should be OUfixedRoot, so I switched to the OUfixedRoot");
        }
    }


    my.edge.list <- cbind(tre$edge-1, tre$edge.length); 
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tree$tip.label));
    return(result);
}

