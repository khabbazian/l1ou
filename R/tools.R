#' Adjusts the tree and traits to meet the requirements of \code{estimate_shift_configuration}
#'
#' Returns a new tree and new data matrix, where the tree edges are in postorder, and the data row names match the order of the tree tip labels.
#'
#'@param tree ultrametric tree of class phylo with branch lengths.
#'@param Y trait vector/matrix without missing entries.
#'@param normalize logical. If TRUE, normalizes branch lengths to a unit tree height.
#'@param quietly logical. If FALSE, changes in tree/trait are printed.
#'
#'@return 
#' \item{tree}{tree of class phylo, with the same topology as the input \code{tree} but adjusted edge order.}
#' \item{Y}{trait vector/matrix with adjusted row names and row order.}
#'@examples
#' data(lizard.tree, lizard.traits)
#' # here, lizard.traits is a matrix, so columns retain row names:
#' names(lizard.traits[,1])
#' lizard <- adjust_data(lizard.tree, lizard.traits[,1])
#' 
#' # for a data frame, make sure to retain row names if a single column is selected:
#' lizard.traits <- as.data.frame(lizard.traits)
#' lizard <- adjust_data(lizard.tree, subset(lizard.traits, select=1))
#'@export
adjust_data <- function(tree, Y, normalize = TRUE, quietly=FALSE){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if( !identical(tree$edge, reorder(tree, "postorder")$edge)){
        if(!quietly)
            cat("the new tree edges are ordered differently: in postorder.\n")
        tree  <- reorder(tree, "postorder")
    }

    if( normalize ){
        tree <- normalize_tree(tree)
        if(!quietly)
            cat("the new tree is normalized: each tip at distance 1 from the root.\n")
    }

    if( class(Y) != "matrix"){
        Y <- as.matrix(Y)
        if(!quietly)
            cat(paste("new Y: matrix of size", nrow(Y), "x", ncol(Y), "\n" ))
    }

    if( nrow(Y) != length(tree$tip.label)){
       stop("the number of entries/rows of the trait vector/matrix (Y) 
            doesn't match the number of tips.\n") 
    }

    if( is.null(rownames(Y)) ){
        warning("no names provided for the trait(s) entries/rows.\nAssuming that rows match the tip labels in the same order.",
                immediate.=TRUE)
        rownames(Y)  <- tree$tip.label
    } else{
        if( any(is.na(rownames(Y))) ){
            stop("some of the names in the trait vector/matrix or in the tree's 
                 tip.label are unavailable.\n")
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
            stop(" do(es) not exist in the input trait. you may want to use
                 drop.tip(tree, setdiff(tree$tip.label,rownames(Y))) 
                 to drop extra tips in the tree.\n")
        }

        if(!quietly)
            cat("reordered the rows of the trait vector/matrix (Y) to match the order of the tip labels.\n")
 
        #Y  <-  Y[order(rownames(Y)),  ] 
        #Y  <-  Y[order(order(tree$tip.label)), ]
        Y <- as.matrix(Y[tree$tip.label, ])

    }


    stopifnot(all(rownames(Y) == tree$tip.label))
    stopifnot(identical(rownames(Y), tree$tip.label))

    return(list(tree=tree, Y=Y))
}

lnorm <- function(v,l=1)   { return( (sum(abs(v)^l, na.rm=TRUE))^(1/l) ) }

## This function is useful for handling missing values in multivariate regression. 
## It generates a list of design matrices and trees considering according to the missing values.
gen_tree_array <- function(tree, Y){ 
    ## here I assume that the tree tip labels match the Y matrix rows 
    ## in the same order.
    tree.list <- list()
    for(trait.idx in 1:ncol(Y)){
        availables <- rownames(Y)[!is.na(Y[,trait.idx])]

        tr <- drop.tip(tree, setdiff(tree$tip.label, availables))
        tr <- reorder(tr, "postorder")

        X.1 <- generate_design_matrix(tree, type="simpX")
        rownames(X.1) <- tree$tip.label
        X.2 <- generate_design_matrix(tr, type="simpX")
        rownames(X.2) <- tr$tip.label

        ## when we drop some tips of a tree edges order changes so we need 
        ## to have a universal mapping for shift configurations.
        old.order <- rep(NA, Nedge(tree))
        for(i in 1:Nedge(tree)){
            tip.set  <- rownames(X.1)[which(X.1[,i]>0)]
            tip.set  <- intersect(tip.set, rownames(X.2)) 
            if(length(tip.set)==0)
                next
            if(length(tip.set)==1)
                edge.set <- which(X.2[tip.set,]==1)
            else
                edge.set <- which(colSums(X.2[tip.set,])==length(tip.set))

            if(length(edge.set) > 1)
                e.idx    <- edge.set[ which(colSums(X.2[,edge.set])==length(tip.set)) ]
            else
                e.idx    <- edge.set 

            stopifnot(length(e.idx)==1)
            old.order[[i]] <- e.idx
        }
        tr$old.order <- old.order
        tr$Z <- X.2
        tree.list[[trait.idx]]  <-  tr
    }
    return(tree.list)
}

add_configuration_score_to_list  <- function(shift.configuration, score, moreInfo){
    shift.configuration = sort(shift.configuration)
    add_configuration_score_to_db( paste0(shift.configuration, collapse=" "), 
                                  score, moreInfo )
}

get_configuration_score_from_list <- function(shift.configuration){
    if(length(shift.configuration) > 0){
        shift.configuration <- sort(shift.configuration)
    }
    res <- get_score_of_configuration(paste0(shift.configuration, collapse=" "))
    if( res$valid == FALSE){
        return(NA)
    }
    return(res$value)
}

list_investigated_configs <- function(){
    tmpList = get_stored_config_score()
    c.s = list()
    c.s$scores = tmpList$scores
    c.s$configurations = lapply(tmpList$configurations, 
                                FUN=function(x) as.numeric(unlist(strsplit(x, split=" ")) ) ) 
    #for( i in 1:length(c.s$scores)){
    #    c.s$configurations[[i]] = as.numeric(unlist(strsplit(tmpList$configurations[[i]], split=" ")) )
    #    #c.s$moreInfo      [[i]] = as.numeric(unlist(strsplit(tmpList$moreInfo      [[i]], split=" ")) )
    #}
    return(c.s)
}

print_out <- function(eModel, quietly){
    if (quietly == FALSE){
        print( paste0( "EST: alpha: ", eModel$alpha, " sigma2: ",  
                 eModel$sigma2, " gamma: ", eModel$sigma2/(2*eModel$alpha),
                 " score: ", eModel$score, " nShifts: ", eModel$nShifts ) )
        print("-------------------")
    }
}


rescale_matrix <- function(Y){
    #for(i in 1:ncol(Y)){
    #    Y[,i] = Y[,i] - mean(Y[,i])
    #}
    #Y   = Y%*%(0.1*nrow(Y)*diag(apply(Y,2,lnorm,l=2)^-1))
  
    Y  <- 0.1*nrow(Y)*scale(Y, center=TRUE, scale=apply(Y,2,lnorm,l=2))
    return(Y)
}



effective.sample.size <- function(phy, edges=NULL,
             model = c("BM","OUrandomRoot","OUfixedRoot","lambda","delta","EB"),
             parameters = NULL, check.pruningwise = TRUE,
             check.ultrametric = TRUE){
        # requires ultrametric tree. Kappa disallowed: causes non-ultrametric tree
        #                            both OU models result in same n_e
        # removes every edge in 'edges' to split phy into m+1 subtrees, then
        # sums log(n_e+1) over all subtrees. n_e = max(V) * one' V^{-1} one
        # where V = phylogenetic covariance matrix for the subtree.
        model = match.arg(model)
        if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
        if (check.pruningwise) phy = reorder(phy,"pruningwise")
        if (check.ultrametric)
            if (!is.ultrametric(phy))
                stop("ultrametric tree required to calculate effective sample sizes.")
        Di <- numeric(length(phy$tip.label)) # zeros
        phy <- transf.branch.lengths(phy,model,parameters=parameters,
                                     check.pruningwise=check.pruningwise,check.ultrametric=FALSE,
                                     D=Di,check.names=F)$tree
        rootedge <- dim(phy$edge)[1]+1
        if (is.null(edges)){ sortededges <- rootedge }
        else{
            o <- order(edges) 
            r <- rank(edges)
            sortededges <- c(edges[o],rootedge)
        }
        tmp <- .C("effectiveSampleSize", as.integer(dim(phy$edge)[1]), # edges
                  as.integer(length(phy$tip.label)), as.integer(phy$Nnode), # tips and nodes
                  as.integer(length(phy$tip.label)+1), # root index
                  as.double(phy$root.edge),as.double(phy$edge.length),
                  as.integer(phy$edge[, 2]), as.integer(phy$edge[, 1]), # descendents and ancestors
                  as.integer(sortededges), # edges to cut, including root edge
                  result=double(length(edges)+1))$result # tmp has, in this order:
        if (is.null(edges))
            res <- tmp
        else res <- tmp[c(length(tmp),r)]
        return(res)
}


correct_unidentifiability <- function(tree, shift.configuration, opt){

    if( length(shift.configuration) < 2)  { return(shift.configuration) }
    shift.configuration = sort(shift.configuration)
    nTips    = length(tree$tip.label)
    nN       = nrow(opt$Z)

    all.covered.tips = numeric()
    identifiable = TRUE
    for(sp in shift.configuration){
        covered.tips = which( opt$Z[,sp] > 0 )
        nUniqueTips = length( setdiff(covered.tips, all.covered.tips) )
        if ( nUniqueTips == 0 ){
            shift.configuration = setdiff(shift.configuration, sp)
            identifiable = FALSE
        }
        all.covered.tips = union(covered.tips, all.covered.tips)
    }

    if( identifiable ){ return(shift.configuration) }

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



get_configuration_in_sol_path <- function(sol.path, index, Y, tidx=1){
    if ( grepl("lars",sol.path$call)[[1]]  ){
        beta      = sol.path$beta[index,]
        shift.configuration  = which( abs(beta) > 0 )
    } else if( any( grepl("grplasso",sol.path$call) ) ){
        #beta  = sol.path$coefficients[, index]
        #lIdx  = length(beta)/ncol(Y)
        #shift.configuration  = which( abs(beta[ (1+(tidx-1)*lIdx):(tidx*lIdx) ]) > 0 )
        ##NOTE: in case, in a group of variables some are zero and some non-zero i consider all as non-zero
        #shift.configuration  = which( rowSums(matrix(abs(beta),nrow=lIdx)) > 0 ) 

        beta = sol.path$coefficients[, index]
        nVariables = ncol(Y)
        MM = matrix(ifelse(abs(beta)>0,1,0), ncol = nVariables)
        shift.configuration = which(rowSums(MM) >= nVariables/2)

    } else {  
        stop(paste0(match.call(), ":undefined solver!"))
    }

    return(shift.configuration)
}



#' Converts shift values to optimum values on the edges.
#'
#' Converts a model indicated with shift values to a model with optimum values on the edges.
#'
#'@param tree ultrametric tree of class phylo with branch lengths.
#'@param shift.configuration vector of edge indices with shifts.
#'@param shift.values vector of shift values.
#'
#'@return vector of size number of edges with optimum value of the trait on the corresponding edge.
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' 
#' sc <- c(55, 98, 118, 74, 14, 77,  32, 164)
#' sv <- c(2 ,  3,   4,  4,  1,  2, 0.5,   1)
#'
#' root.value <- -2
#'
#' optimum.values <- convert_shifts2regions(lizard.tree, sc, sv) + root.value
#'
#'
#'@export
convert_shifts2regions <-function(tree, shift.configuration, shift.values){

    stopifnot( length(shift.configuration) == length(shift.values) )

    nTips   = length(tree$tip.label)
    nEdges  = Nedge(tree)
    g       = graph.edgelist(tree$edge, directed = TRUE)
    o.vec = rep(0, nEdges)

    prev.val <-options()$warn 
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
    options(warn = prev.val)
    return( o.vec )
}

#' Normalizes branch lengths to a unit tree height
#'
#' Normalizes all branch lengths by the same factor, so that the distance from the root to all tips is equal to one. 
#'@param tree ultrametric tree of class phylo with branch lengths, and edges in postorder.
#'@param check.ultrametric logical. If TRUE, it checks if the input tree is ultrametric.
#'
#'@return normalized phylogenetic tree, of class phylo.
#'
#'@export
normalize_tree <- function(tree, check.ultrametric=TRUE){

    if(check.ultrametric){
        if(!is.ultrametric(tree)) 
            stop("the input tree is not ultrametric")
    }

    nTips  = length(tree$tip.label)
    rNode  = nTips + 1 
    nEdges = Nedge(tree)

    g        = graph.edgelist(tree$edge, directed = TRUE)
    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode="out", output="epath")$epath

    root.edge <- ifelse(is.null(tree$root.edge), 0, tree$root.edge)
    Tval     = root.edge + sum(tree$edge.length[root2tip[[1]] ])
    #Tval = mean ( sapply( 1:nTips, FUN=function(x) sum(tree$edge.length[root2tip[[x]]])   )  )

    tree$edge.length = tree$edge.length / Tval
    if(!is.null(tree$root.edge)){
	    tree$root.edge <- tree$root.edge / Tval 
    }
    return(tree)
}


#'
#' Visualizes a shift configuration: tree and trait(s)
#'
#' plots the tree annotated to show the edges with a shift, and the associated trait data side by side.
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param palette vector of colors, of size the number of shifts plus one. The last element is the color for the background regime (regime at the root).
#'@param edge.shift.ann logical. If TRUE, annotates edges by shift values. 
#'@param edge.shift.adj adjustment argument to give to edgelabel() for labeling edges by shift values.
#'@param asterisk logical. If TRUE, the shift positions will be annotated by "*". It is useful for gray scale plots.
#'@param edge.label vector of size number of edges.
#'@param edge.label.ann logical. If TRUE, annotates edges by labels in tree$edge.label, if non-empty, or edge.label. 
#'@param edge.label.adj adjustment argument to give to edgelabel() for labeling edges.
#'@param edge.label.pos relative position of the edge.label on the edge. 0 for the beginning of the edge and 1 for the end of the edge. 
#'@param edge.ann.cex amount by which the annotation text should be magnified relative to the default.
#'@param plot.bar logical. If TRUE, the bars corresponding to the trait values will be plotted.
#'@param bar.axis logical. If TRUE, the axis of of trait(s) range will be plotted. 
#'@param ... further arguments to be passed on to plot.phylo. 
#'
#'@return none.
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1]
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' nEdges <- Nedge(lizard.tree)
#' ew <- rep(1,nEdges) 
#' ew[eModel$shift.configuration] <- 3
#' plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#'@export
#'
plot.l1ou <- function (model, palette = NA, 
                       edge.shift.ann=TRUE,  edge.shift.adj=c(0.5,-.025),
                       edge.label=c(), asterisk = TRUE,
                       edge.label.ann=FALSE, edge.label.adj=c(0.5,    1), 
                       edge.label.pos=NA,
                       edge.ann.cex = 1, 
                       plot.bar = TRUE, bar.axis = TRUE, ...) 
{
    tree = model$tree
    s.c = model$shift.configuration
    stopifnot(identical(tree$edge, reorder(tree, "postorder")$edge))
    nShifts = model$nShifts
    nEdges = Nedge(tree)
    if (bar.axis) 
        par(oma = c(3, 0, 0, 3))
    Y = as.matrix(model$Y)
    stopifnot(identical(rownames(Y), tree$tip.label))

    if (bar.axis) 
        par(oma = c(3, 0, 0, 3))

    if (plot.bar) {
        layout(matrix(c(1+ncol(Y),1:ncol(Y)), nrow=1), 
               widths = c(2,rep(1, ncol(Y)))
               )
    }

    #NOTE: assiging colors the shifts
    if (all(is.na(palette))) {
        palette = c(sample(rainbow(nShifts)), "gray")
        if( !is.null(names(s.c)) ){
            ids = unique(names(s.c))
            tmp = sample(rainbow(length(ids)))
            for( id in ids )
                palette[which(names(s.c)==id)] = tmp[which(ids==id)]
        }
    }

    stopifnot(length(palette) == model$nShifts + 1)

    edgecol = rep(palette[nShifts + 1], nEdges)
    counter = 1
    Z = model$l1ou.options$Z
    if(length(s.c) > 0)
        for (shift in sort(s.c, decreasing = T)) {
            edgecol[[shift]] = palette[[which(s.c == shift)]]
            tips = which(Z[, shift] > 0)
            for (tip in tips) {
                edgecol[which(Z[tip, 1:shift] > 0)] = palette[[which(s.c == 
                    shift)]]
            }
            counter = counter + 1
        }



    ##A dummy plot just to get the plotting order
    plot.phylo(tree, plot=FALSE)
    lastPP = get("last_plot.phylo", envir = .PlotPhyloEnv)
    o = order(lastPP$yy[1:length(tree$tip.label)])
    par.new.default <- par()$new ##just to be careful with the global variable
    par(new=TRUE)

    #NOTE: plotting bar plot .....
    if (plot.bar) {
        nTips = length(tree$tip.label)
        barcol = rep("gray", nTips)
        for (i in 1:nTips) {
            barcol[[i]] = edgecol[which(tree$edge[, 2] == i)]
        }
        if (bar.axis) 
            par(mar = c(0, 0, 0, 3))
        for (i in 1:ncol(Y)) {
            normy = (Y[, i] - mean(Y[, i], na.rm=TRUE))/sd(Y[, i], na.rm=TRUE)
            barplot(as.vector(normy[o]), border = FALSE, col = barcol[o], 
                    horiz = TRUE, names.arg = "", xaxt = "n")
            if (bar.axis){
                axis(1, at = range(normy, na.rm=TRUE), 
                     labels = round(range(normy, na.rm=TRUE), 
                                    digits = 2))
            }
            if (!is.null(colnames(Y)) && length(colnames(Y)) > 
                (i - 1)) 
                mtext(colnames(Y)[[i]], cex = 1, line = +1, side = 1)
        }
    }

    #NOTE: plotting the tree etc etc
    plot.phylo(tree, edge.color = edgecol, no.margin = TRUE, 
        ...)

    if (length(s.c) > 0) {
        if (asterisk) {
            Z = l1ou:::generate_design_matrix(tree, type = "apprX")
            for (idx in 1:length(s.c)) {
                sP = s.c[[idx]]
                pos = max(Z[, sP])
                edge.labels = rep(NA, length(tree$edge[, 1]))
                edge.labels[sP] = "*"
                edgelabels(edge.labels, cex = 3 * edge.ann.cex, 
                  adj = c(0.5, 0.8), frame = "none", date = pos)
            }
        }
    }
    if (edge.shift.ann) {
        eLabels = rep(NA, nEdges)
        for (shift in s.c) {
            eLabels[shift] = paste(round(model$shift.values[which(s.c == 
                shift), ], digits = 2), collapse = ",")
        }
        edgelabels(eLabels, cex = edge.ann.cex, adj = edge.shift.adj, 
            frame = "none")
    }
    if (edge.label.ann) {
        if (length(tree$edge.label) == 0) {
            if (length(edge.label) == 0) {
                stop("no edge labels are provided via tree$edge.label or edge.label!")
            }
            tree$edge.label = edge.label
        }
        Z = l1ou:::generate_design_matrix(tree, type = "apprX")
        if (!is.na(edge.label.pos)) 
            if (edge.label.pos < 0 || edge.label.pos > 1) 
                stop("edge.label.pos should be between 0 and 1")
        for (idx in 1:length(tree$edge.label)) {
            if (is.na(tree$edge.label[[idx]])) 
                next
            pos = max(Z[, idx])
            if (!is.na(edge.label.pos)) {
                pos = pos - edge.label.pos * tree$edge.length[[idx]]
            }
            edge.labels = rep(NA, length(tree$edge[, 1]))
            edge.labels[[idx]] = tree$edge.label[[idx]]
            edgelabels(edge.labels, cex = edge.ann.cex, adj = edge.label.adj, 
                frame = "none", date = pos)
        }
    }
    par(new=par.new.default)
}

#'
#' Prints out a summary of the shift configurations investigated by \code{\link{estimate_shift_configuration}}  
#'
#' prints the list of the shift configurations sorted by number of shifts and corresponding ic scores.
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param ... further arguments. 
#'
#'@return 
#'\item{shift.configurations}{list of shift configurations sorted by number of shifts.}
#'\item{scores}{list of scores corresponding to shift.configurations.}
#'\item{nShifts}{number of shifts corresponding to the shift configurations.}
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1]
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' model.profile  <- profile(eModel)
#' plot(model.profile$nShifts, model.profile$scores)
#'
#'@export
#'
profile.l1ou <- function(model, ...)
{

    profile.data = model$profile
    p.d = list()
    profile.data$scores = profile.data$scores[order(profile.data$scores)]
    profile.data$configurations = profile.data$configurations[order(profile.data$scores)]
    lens = unlist(lapply(profile.data$configurations, length))
    profile.data$scores = profile.data$scores[order(lens)]
    profile.data$configurations = profile.data$configurations[order(lens)]
    min.score = min(profile.data$scores)
    clength = -1
    counter = 1
    for (i in 1:length(profile.data$scores)) {
        if (clength == length(profile.data$configurations[[i]])) {
            next
        }
        clength = length(profile.data$configurations[[i]])
        p.d$shift.configurations[[counter]] = profile.data$configurations[[i]]

        p.d$nShifts[[counter]] = length(profile.data$configurations[[i]])
        p.d$scores [[counter]] = profile.data$score[[i]]
        #p.d$gamma  [[counter]] = profile.data$moreInfo[[i]][[1]] ##the stationary variance
        #p.d$logLik [[counter]] = profile.data$moreInfo[[i]][[2]] ##the log likelihood

        counter = counter + 1
    }
    return(p.d)
}

#'
#' Returns the best shift configuration with a given number of shifts among the shift configurations that have been evaluated.
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param nShifts number of shifts.
#'
#'@return indices of the edges with shifts
#'
#'@export
get_shift_configuration <- function(model, nShifts){
    p.d = profile(model) 
    if( nShifts > length(p.d$shift.configurations)+1) # starts at 0 shifts
        stop("There is no configuration with the given number of shifts")

    for( i in 1:length(p.d$shift.configurations)){
        if( length(p.d$shift.configurations[[i]]) == nShifts)
            return(p.d$shift.configurations[[i]])
    }
    stop("There is no configuration with the given number of shifts")
}

#'
#' Prints out a summary of the model 
#'
#' prints out a summary of the model 
#'
#'@param model object of class l1ou returned by \code{\link{estimate_shift_configuration}}.
#'@param nTop.scores number of top scores and shift configuration to print out.
#'@param ... further arguments. 
#'
#'@return none.
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' Y <- lizard.traits[,1]
#' eModel <- estimate_shift_configuration(lizard.tree, Y)
#' summary(eModel)
#'
#'@export
#'
summary.l1ou <- function(model, nTop.scores=5, ...){
    cat("number of shifts: ")
    cat(model$nShifts)
    cat("\n")

    cat("edge indices of the shift configuration (column names) and the corresponding shift values:\n")
    #cat(model$shift.configuration)
    #cat(model$shift.values)
    #cat("\n")
    tmp.mat = t(as.matrix(model$shift.values))
    if(length(model$shift.configuration)>0)
        colnames(tmp.mat) = model$shift.configuration
    if(!all(is.null(colnames(model$Y)))){
        rownames(tmp.mat) = colnames(model$Y)
    }
    print(tmp.mat)

    cat("\n")
    cat(paste0(model$l1ou.options$criterion, " score: "))
    cat(model$score)
    cat("\n")

    tmp.mat = rbind(model$alpha, 
                    model$sigma2, 
                    model$sigma2/(2 * model$alpha),
                    model$logLik
                    )
    rownames(tmp.mat) = c("adaptation rate (alpha)", 
                          "variance (sigma2)", 
                          "stationary variance (gamma)",
                          "logLik"
                          )
    if(!all(is.null(colnames(model$Y)))){
        colnames(tmp.mat) = colnames(model$Y)
    }
    print(tmp.mat)
    cat("\n")

    cat("\n")
    cat("optimum values at tips: \n")
    print(model$optimums)

    top.scores = min(nTop.scores, length(model$profile$scores))
    cat(paste0(c("\ntop", top.scores, "best scores among candidate models evaluated during the search:\n")))
    cat("scores\t\tshift.configurations\n")
    for (i in 1:top.scores){
        cat(model$profile$scores[[i]])
        cat("\t")
        cat(model$profile$configurations[[i]])
        cat("\n")
    }
}

#'@export
print.l1ou <- function(model, ...){
    cat("number of shifts: ")
    cat(model$nShifts)
    cat("\n")

    cat(paste0(model$l1ou.options$criterion, " score: "))
    cat(model$score)
    if(!is.null(model$cr.score)){
        cat("\n")
        cat(paste0(model$l1ou.options$criterion, " CR score: "))
        cat(model$cr.score)
    }
    cat("\n")

    cat("edge indices of the shift configuration (column names) and the corresponding shift values:\n")

    if( length(model$shift.configuration) > 0){
        tmp.mat = t(as.matrix(model$shift.values))
        if(length(model$shift.configuration)>0)
            colnames(tmp.mat) = model$shift.configuration
        if(!all(is.null(colnames(model$Y)))){
            rownames(tmp.mat) = colnames(model$Y)
        }
        print(tmp.mat)
        cat("\n")
    }


    sc <- model$shift.configuration
    if( !is.null(names(sc)) ){
        cat("convergent regimes and edge indices of the shift configuration\n")
        for( reg in sort(unique(names(sc))) )
            cat( paste0( "regime ", reg, "-> ", paste0(sc[which(names(sc)==reg)], collapse=", "), "\n" ) )
        cat("\n")
    }


    tmp.mat <- rbind(model$alpha, 
                    model$sigma2, 
                    model$sigma2/(2 * model$alpha),
                    model$logLik
                    )
    rownames(tmp.mat) <- c("adaptation rate (alpha)", 
                          "variance (sigma2)", 
                          "stationary variance (gamma)",
                          "logLik"
                          )
    if(!all(is.null(colnames(model$Y)))){
        colnames(tmp.mat) <- colnames(model$Y)
    }
    print(tmp.mat)
    cat("\n")
}

