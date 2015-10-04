


lnorm          <- function(v,l=1)   { return( (sum(abs(v)^l))^(1/l) ) };
my.time.format <-function()         { return(format(Sys.time(),"%y_%m_%d_%H_%M")); }


add.placement.score.to.list  <- function(shift.placement, score){
    shift.placement = sort(shift.placement);
    add_placement_score_to_db(paste0(shift.placement, collapse=" "), score);
}

get.placement.score.from.list <- function(shift.placement){
    shift.placement = sort(shift.placement);
    res = get_score_of_placement(paste0(shift.placement, collapse=" "));
    if( res$valid == FALSE){
        return(NA);
    }
    return(res$value);
}



print.out.param <- function(eModel, silence){
    if ( silence == FALSE)
        print( paste0( "EST: alpha: ", eModel$alpha, " sigma2: ",  
                 eModel$sigma2, " gamma: ", eModel$sigma2/(2*eModel$alpha),
                 " score: ", eModel$score ) );
}





standardize.matrix <- function(Y){
    for(i in 1:ncol(Y)){
        Y[,i] = Y[,i] - mean(Y[,i]);
    }
    Y   = Y%*%(0.1*nrow(Y)*diag(apply(Y,2,lnorm,l=2)^-1));
    return(Y);
}



correct.unidentifiability <- function(tr, shift.placement, opt){

    if( length(shift.placement) < 2)  { return(shift.placement); }
    shift.placement = sort(shift.placement);
    nTips    = length(tr$tip.label);
    nN       = nrow(opt$Z);

    all.covered.tips = numeric();
    for(sp in shift.placement){
        covered.tips = which( opt$Z[,sp] > 0 );
        nUniqueTips = length( setdiff(covered.tips, all.covered.tips) );
        if ( nUniqueTips == 0 )
            shift.placement = setdiff(shift.placement, sp);
        all.covered.tips = union(covered.tips, all.covered.tips);
    }

    while ( length(shift.placement) > 1 ) {
        coverage = c();
        for(sp in shift.placement)
            coverage = union( coverage, which(opt$Z[,sp] > 0) );

        if( length( setdiff(1:nN, coverage) ) == 0 ){
            shift.placement = shift.placement[-which.max(shift.placement)];
        } else { break; }
    }
    return(shift.placement);
}



alpha.upper.bound <- function(tr){
    nTips       = length(tr$tip.label);
    #eLenSorted  = sort(tr$edge.length[which(tr$edge[,2] < nTips)]); 
    #topMinLen   = ceiling( length(eLenSorted)*(5/100) );
    #return( log(2)/median(eLenSorted[1:topMinLen]) );
    return( log(2)/min(tr$edge.length[which(tr$edge[,2] < nTips)]) );
}



get.num.solutions <- function(sol.path){
    if ( grepl("lars",sol.path$call)[[1]]  ){
        return ( length(sol.path$beta[,1]) );
    } 

    if ( any( grepl("grplasso",sol.path$call) ) ){
        return ( length(sol.path$coefficients[1,]) );
    }
    stop(paste0(match.call(), ":undefined solver!"));
}

get.shift.placement <- function(sol.path, index, Y, tidx=1){
    if ( grepl("lars",sol.path$call)[[1]]  ){
        beta      = sol.path$beta[index,];
        shift.placement  = which( abs(beta) > 0 );
    } else if( any( grepl("grplasso",sol.path$call) ) ){
        beta  = sol.path$coefficients[, index];
        lIdx  = length(beta)/ncol(Y);
        #shift.placement  = which( abs(beta[ (1+(tidx-1)*lIdx):(tidx*lIdx) ]) > 0 );
        ##NOTE: in case, in a group of variables some are zero and some non-zero i consider all as non-zero
        shift.placement  = which( rowSums(matrix(abs(beta),nrow=lIdx)) > 0 ) ;
    } else {  
        stop(paste0(match.call(), ":undefined solver!"));
    }

    return(shift.placement);
}



convert.shifts2regions <-function(tr, shift.placement, shift.values){

    stopifnot( length(shift.placement) == length(shift.values) );

    nTips   = length(tr$tip.label);
    nEdges  = length(tr$edge.length);
    g       = graph.edgelist(tr$edge, directed = TRUE);
    o.vec = rep(0, nEdges);

    options(warn = -1);
    if( length(shift.placement) > 0)
    for(itr in 1:length(shift.placement) ){
        eIdx     = shift.placement[[itr]];
        vIdx     = tr$edge[eIdx, 2];

        o.vec.tmp = rep(0, nEdges);

        path2tips = get.shortest.paths(g, vIdx, to=1:nTips, mode ="out", output="epath")$epath;
        o.vec.tmp[eIdx] = shift.values[[itr]];
        for(i in 1:nTips){
            o.vec.tmp[path2tips[[i]]] =  shift.values[[itr]];
        }
        o.vec = o.vec + o.vec.tmp; 
    }
    options(warn = 0);

    return( o.vec );
}

#' Normalizes the branch lengths so that the distance from the root to all tips are equal to one. 
#'@param tr an ultrametric phylogeny tree.
#'
#'@return returns normalized phylogeny tree.
#'
#'@export
normalize.tree <- function(tr){
    stopifnot(is.ultrametric(tr));

    nTips  = length(tr$tip.label);
    rNode  = nTips + 1; 
    nEdges = Nedge(tr);

    g        = graph.edgelist(tr$edge, directed = TRUE);
    root2tip = get.shortest.paths(g, rNode, to=1:nTips, mode="out", output="epath")$epath;

    Tval     = sum(tr$edge.length[root2tip[[1]] ]);
    tr$edge.length = tr$edge.length / Tval;
    return(tr);
}


my.plot.tree <-function(tr, opt.val=numeric(), plot.title="", colvec=c(),
        out.fn="unnamed", show.el=FALSE, intercept = 0,
        edge.labels=numeric(), plotme=T, 
        el.center=FALSE, nomargins = TRUE, el.cex=0.5, ...){

    if( length(opt.val) > 0 ){
        nConvReg  = length(unique(opt.val));

        if( length(colvec) == 0)
            colvec    = rainbow( length(unique(opt.val)) );

        stopifnot( length(colvec) == length( unique(opt.val) ) );
        edgecol = rep(0, length(opt.val));
        i = 1;
        for (itr in unique(opt.val)){
            if(itr == intercept){
                edgecol[which(opt.val==itr)] = "gray";
            }else{
                edgecol[which(opt.val==itr)] = colvec[[i]];
            }
            i = i + 1;
        }

    }else{edgecol="black";}

    if(!plotme){
    return(edgecol);
    }


    plot.phylo(tr, edge.color=edgecol, no.margin=nomargins, ...);
    title(plot.title, cex=0.8 );
    if(show.el==TRUE){
        if ( length( edge.labels) > 0  ){
            if ( el.center == TRUE){ 
                edgelabels(edge.labels, cex=el.cex, 
                           frame = "none", bg="lightblue");
            } else {
                edgelabels(edge.labels, adj = c(0.5, -0.25), cex=el.cex, 
                           frame = "none", bg="lightblue");
            }
        } else{ edgelabels(cex=el.cex, frame = "circle");}
    }

    return(edgecol);
}

my.plot.phylo <- function(tr, Y, eModel, title.str="", enable.cross=FALSE, ...){

    Y = as.matrix(Y);
    stopifnot(identical(rownames(Y), tr$tip.label));

    layout(matrix(1:(1+ncol(Y)), 1, (1+ncol(Y))));

    if ( ncol(Y) == 1){
        edge.labels = rep(NA, length(eModel$opt.val));
        if ( enable.cross == TRUE){
            if( eModel$nShifts > 0 )
                edge.labels[ eModel$shift.placement ] = "X";
        } else{
            edge.labels[ eModel$shift.placement ] = round(eModel$shift.values, digits = 2);
        }
        my.plot.tree(tr, eModel$opt.val, show.el=TRUE, edge.labels = edge.labels, nomargins=FALSE, edge.width=3, ...);
    } else {
        edge.labels = rep(NA, length(eModel$opt.val));
        if ( enable.cross == TRUE){
            if( eModel$nShifts > 0 )
                edge.labels[ eModel$shift.placement ] = "X";
                my.plot.tree(tr, eModel$opt.val[,1], show.el=TRUE, 
                        edge.labels = edge.labels, nomargins = FALSE, 
                        el.center=TRUE, ...);
        } else {
          if( eModel$nShifts > 0 )
              edge.labels[ eModel$shift.placement ] = 
                  apply( round(eModel$shift.values,2), 1, function(x) paste0(x, collapse = ", "));
          my.plot.tree(tr, eModel$opt.val[,1], show.el=TRUE, 
                       edge.labels = edge.labels, nomargins = FALSE, ...);
          
        }
    }

    for(i in 1:ncol(Y)){
        normy = (Y[,i] - mean(Y[,i]))/sd(Y[,i]);
        barplot(as.vector(normy), horiz = TRUE, names.arg = "", xaxt = "n");
        axis(1, at = range(normy), labels = round(range(normy), digits = 2));
        if ( length( title.str) > 0){
            title(title.str[[i]], cex.main = 0.8);
        }
    }
}


assign.color <- function(theRef){
    colvec    = rainbow( length(unique(theRef)) )
        col = rep(0, length(unique(theRef)));
    i = 1;
    for (itr in unique(theRef)){
        col[which(theRef==itr)] = colvec[[i]];
        i = i + 1;
    }
    return(col);
}

