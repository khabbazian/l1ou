
#
#' Nagative square root and square root of covaiance matrix. 
#'
#'@param tr The input phylogeny
#'@param alpha The adaptation rate.
#'@param root.model The model of phylogeny ancestoral state.
#'
#'@return Negative square root and square root of the phylogeny covariance matrix.
#'
#'@export
#'
sqrt.ou.covariance <- function(tr, alpha=0, root.model = c("OUrandomRoot", "OUfixedRoot")){

    tr         <- multi2di(tr, random=FALSE);
    root.model <- match.arg(root.model); 
    tr         <- reorder(tr, "prun");

    if ( alpha > 0){
        tre <- transf.branch.lengths(tr, model=root.model, parameters=list(alpha=alpha))$tree;
    } else{
        tre <- tr;
    }
    my.edge.list <- cbind(tre$edge-1, tre$edge.length); 
    result <- cmp_sqrt_OU_covariance(my.edge.list, length(tr$tip.label));
    return(result);
}

