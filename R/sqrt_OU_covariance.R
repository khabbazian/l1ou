
#
#' computes the nagative square root and square root of the phylogeny covaiance matrix. 
#'
#'@param tr the input ultrametric phylogeny.
#'@param alpha the adaptation rate for OU model.
#'@param root.model the model of phylogeny ancestoral state.
#'
#'@return returns negative square root and square root of the phylogeny covariance matrix, i.e. \eqn{ \Sigma^{-1/2} } and \eqn{ \Sigma^{1/2}}.
#'
#'@examples
#'
#' library("l1ou")
#' data("lizard.tree")
#' res <- sqrt_OU_covariance(lizard.tree);
#' Sigma <- vcv(lizard.tree)
#' dimnames(Sigma) <- NULL
#' all.equal(res$sqrtSigma %*% t(res$sqrtSigma) , Sigma) # TRUE
#' all.equal(res$sqrtInvSigma %*% t(res$sqrtInvSigma) , solve(Sigma)) # TRUE
#'
#'@export
#'
sqrt_OU_covariance <- function(tr, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot") ){

    tr         <- multi2di(tr, random=FALSE);
    root.model <- match.arg(root.model); 
    tr         <- reorder(tr, "prun");

    if ( alpha > 0){
        tre <- transf.branch.lengths(tr, model=root.model, parameters=list(alpha=alpha))$tree;
    } else{
        tre <- tr;
        if( root.model == "OUrandomRoot"){
            warning("when alpha is zero the model should be OUfixedRoot, so I switched to the OUfixedRoot");
        }
    }


    my.edge.list <- cbind(tre$edge-1, tre$edge.length); 
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tr$tip.label));
    return(result);
}

