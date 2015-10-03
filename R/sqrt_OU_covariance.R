library("Rcpp");
library("ape");
library("phylolm");

sqrt.ou.covariance <- function(tre0, alpha=0, root.model = c("OUrandomRoot", "OUfixedRoot")){

    tre0       <- multi2di(tre0, random=FALSE);
    root.model <- match.arg(root.model); 
    tre0       <- reorder(tre0, "prun");

    if ( alpha > 0){
        tre <- transf.branch.lengths(tre0, model=root.model, parameters=list(alpha=alpha))$tree;
    } else{
        tre <- tre0;
    }
    my.edge.list <- cbind(tre$edge-1, tre$edge.length); 
    result <- cmp_sqrt_OU_covariance(my.edge.list, length(tre0$tip.label));
    return(result);
}
