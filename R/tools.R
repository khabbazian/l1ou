


lnorm          <- function(v,l=1)   { return( (sum(abs(v)^l))^(1/l) ) }
#my.time.format <-function()         { return(format(Sys.time(),"%y_%m_%d_%H_%M")) }


add_configuration_score_to_list  <- function(shift.configuration, score){
    shift.configuration = sort(shift.configuration)
    add_configuration_score_to_db(paste0(shift.configuration, collapse=" "), score)
}

get_configuration_score_to_list <- function(shift.configuration){
    shift.configuration = sort(shift.configuration)
    res = get_score_of_configuration(paste0(shift.configuration, collapse=" "))
    if( res$valid == FALSE){
        return(NA)
    }
    return(res$value)
}



print_out <- function(eModel, silence){
    if ( silence == FALSE)
        print( paste0( "EST: alpha: ", eModel$alpha, " sigma2: ",  
                 eModel$sigma2, " gamma: ", eModel$sigma2/(2*eModel$alpha),
                 " score: ", eModel$score ) )
}





standardize_matrix <- function(Y){
    for(i in 1:ncol(Y)){
        Y[,i] = Y[,i] - mean(Y[,i])
    }
    Y   = Y%*%(0.1*nrow(Y)*diag(apply(Y,2,lnorm,l=2)^-1))
    return(Y)
}



correct_unidentifiability <- function(tree, shift.configuration, opt){

    if( length(shift.configuration) < 2)  { return(shift.configuration) }
    shift.configuration = sort(shift.configuration)
    nTips    = length(tree$tip.label)
    nN       = nrow(opt$Z)

    all.covered.tips = numeric()
    for(sp in shift.configuration){
        covered.tips = which( opt$Z[,sp] > 0 )
        nUniqueTips = length( setdiff(covered.tips, all.covered.tips) )
        if ( nUniqueTips == 0 )
            shift.configuration = setdiff(shift.configuration, sp)
        all.covered.tips = union(covered.tips, all.covered.tips)
    }

    while ( length(shift.configuration) > 1 ) {
        coverage = c()
        for(sp in shift.configuration)
            coverage = union( coverage, which(opt$Z[,sp] > 0) )

        if( length( setdiff(1:nN, coverage) ) == 0 ){
            shift.configuration = shift.configuration[-which.max(shift.configuration)]
        } else { break }
    }
    return(shift.configuration)
}



alpha_upper_bound <- function(tree){
    nTips       = length(tree$tip.label)
    #eLenSorted  = sort(tree$edge.length[which(tree$edge[,2] < nTips)]) 
    #topMinLen   = ceiling( length(eLenSorted)*(5/100) )
    #return( log(2)/median(eLenSorted[1:topMinLen]) )
    return( log(2)/min(tree$edge.length[which(tree$edge[,2] < nTips)]) )
}



get_num_solutions <- function(sol.path){
    if ( grepl("lars",sol.path$call)[[1]]  ){
        return ( length(sol.path$beta[,1]) )
    } 

    if ( any( grepl("grplasso",sol.path$call) ) ){
        return ( length(sol.path$coefficients[1,]) )
    }
    stop(paste0(match.call(), ":undefined solver!"))
}



get_shift_configuration <- function(sol.path, index, Y, tidx=1){
    if ( grepl("lars",sol.path$call)[[1]]  ){
        beta      = sol.path$beta[index,]
        shift.configuration  = which( abs(beta) > 0 )
    } else if( any( grepl("grplasso",sol.path$call) ) ){
        beta  = sol.path$coefficients[, index]
        lIdx  = length(beta)/ncol(Y)
        #shift.configuration  = which( abs(beta[ (1+(tidx-1)*lIdx):(tidx*lIdx) ]) > 0 )
        ##NOTE: in case, in a group of variables some are zero and some non-zero i consider all as non-zero
        shift.configuration  = which( rowSums(matrix(abs(beta),nrow=lIdx)) > 0 ) 
    } else {  
        stop(paste0(match.call(), ":undefined solver!"))
    }

    return(shift.configuration)
}



convert_shifts2regions <-function(tree, shift.configuration, shift.values){

    stopifnot( length(shift.configuration) == length(shift.values) )

    nTips   = length(tree$tip.label)
    nEdges  = length(tree$edge.length)
    g       = graph.edgelist(tree$edge, directed = TRUE)
    o.vec = rep(0, nEdges)

    options(warn = -1)
    if( length(shift.configuration) > 0)
    for(itr in 1:length(shift.configuration) ){
        eIdx     = shift.configuration[[itr]]
        vIdx     = tree$edge[eIdx, 2]

        o.vec.tmp = rep(0, nEdges)

        path2tips = get.shortest.paths(g, vIdx, to=1:nTips, mode ="out", output="epath")$epath
        o.vec.tmp[eIdx] = shift.values[[itr]]
        for(i in 1:nTips){
            o.vec.tmp[path2tips[[i]]] =  shift.values[[itr]]
        }
        o.vec = o.vec + o.vec.tmp 
    }
    options(warn = 0)
    return( o.vec )
}

#' Normalizes branch lengths to a unit tree height
#'
#' Normalizes all branch lengths by the same factor, so that the distance from the root to all tips is equal to one. 
#'@param tree ultrametric tree of class phylo with branch lengths, and edges in postorder.
#'
#'@return normalized phylogenetic tree, of class phylo.
#'
#'@export
normalize_tree <- function(tree){
    stopifnot(is.ultrametric(tree))

    nTips  = length(tree$tip.label)
    rNode  = nTips + 1 
    nEdges = Nedge(tree)

    g        = graph.edgelist(tree$edge, directed = TRUE)
    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode="out", output="epath")$epath

    Tval     = sum(tree$edge.length[root2tip[[1]] ])
    tree$edge.length = tree$edge.length / Tval
    return(tree)
}


#'
#' Visualizes a shift configuration: tree and trait(s).
#'
#' plots the tree annotated to show the edges with a shift, and the associated trait data side by side.
#'
#'@param tree phylogenetic tree of class phylo.
#'@param model object returned by \code{\link{estimate_shift_configuration}}.
#'@param palette vector of colors, of size the number of shifts plus one. The last element is the color for the background regime (regime at the root).
#'@param edge.ann logical. If TRUE, annotates edges with shift values or by labels in tree$edge.label, if non-empty. 
#'@param edge.ann.cex amount by which the annotation text should be magnified relative to the default.
#'@param plot.bar logical. If TRUE, the bars corresponding to the trait values will be plotted.
#'@param bar.axis logical. If TRUE, the axis of of trait(s) range will be plotted. 
#'@param ... further arguments to be passed on to plot.phylo 
#'
#'@return none.
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1]
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' nEdges <- length(lizard.tree$edge[,1]);
#' ew <- rep(1,nEdges) 
#' ew[eModel$shift.configuration] <- 3
#' plot_l1ou(lizard.tree, eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#'@export
#'
plot_l1ou <- function (tree, model, palette = NA, edge.ann = TRUE, 
                       edge.ann.cex = 1, plot.bar = TRUE, bar.axis = TRUE, ...) 
{
    stopifnot(identical(tree$edge, reorder(tree, "postorder")$edge))
    shift.configuration = sort(model$shift.configuration, decreasing = T)
    nShifts = model$nShifts
    nEdges = length(tree$edge.length)
    if (bar.axis) 
        par(oma = c(3, 0, 0, 3))
    Y = as.matrix(model$Y)
    stopifnot(identical(rownames(Y), tree$tip.label))

    if( plot.bar){
        layout(matrix(1:(1 + ncol(Y)), 1, (1 + ncol(Y))), width = c(2, 1, 1, 1, 1))
    }

    if (is.na(palette)) {
        palette = c(sample(rainbow(nShifts)), "gray")
    }

    stopifnot(length(palette) == model$nShifts + 1)
    edgecol = rep(palette[nShifts + 1], nEdges)
    counter = 1
    Z = model$l1ou.options$Z
    for (shift in model$shift.configuration) {
        edgecol[[shift]] = palette[[counter]]
        tips = which(Z[, shift] > 0)
        for (tip in tips) {
            edgecol[which(Z[tip, 1:shift] > 0)] = palette[[counter]]
        }
        counter = counter + 1
    }
    plot.phylo(tree, edge.color = edgecol, no.margin = TRUE, ...)

    if( edge.ann ){
        if( length(tree$edge.label) == 0){
            tree$edge.label = rep(NA, nEdges)
            counter <- 1
            for (shift in shift.configuration) {
                tree$edge.label[shift] = paste(round(model$shift.values[counter, ], digits = 2), collapse = ",")
                counter <- counter + 1
            }
        }
        edgelabels(tree$edge.label, cex=edge.ann.cex, adj = c(0.5, -0.25), frame = "none")
    }

    if(plot.bar){
        nTips = length(tree$tip.label)
        barcol = rep("gray", nTips)
        for (i in 1:nTips) {
            barcol[[i]] = edgecol[which(tree$edge[, 2] == i)]
        }
        if (bar.axis) 
            par(mar = c(0, 0, 0, 3))
        for (i in 1:ncol(Y)) {
            normy = (Y[, i] - mean(Y[, i]))/sd(Y[, i])
            barplot(as.vector(normy), border = FALSE, col = barcol, 
                    horiz = TRUE, names.arg = "", xaxt = "n")
            if (bar.axis) 
                axis(1, at = range(normy), labels = round(range(normy), 
                                                          digits = 2))
        }
    }
}
