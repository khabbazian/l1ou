#TODO implement me in cpp
cmp.OU.covariance <- function(tre0, alpha=0, root.model = 
                              c("OUrandomRoot", "OUfixedRoot") ){
    
  tre0 = multi2di(tre0, random=FALSE);
  root.model  = match.arg(root.model);
  tre0 = reorder(tre0, "prun");
  
  initialize = function(tre){
    # external variables are used and modified: F,G,D,B,tips
    n = length(tre$tip.label)
    F <<- cbind(diag(1,n), matrix(0,n,n-1))
    G <<- cbind(diag(1,n), matrix(0,n,n-1))
    D <<- NULL; B <<-NULL
    tips <<- 1:n
  }
  
  onestep = function(i1,i2){
    # external variables are used and modified: F,G,D,B,tre,tips
    e1 = which(tre$edge[,2]==i1)
    e2 = which(tre$edge[,2]==i2)
    i3  = tre$edge[e1,1];
    if (i3 != tre$edge[e2,1]){ # check e1 and e2 have same parent
      stop("Error: the 2 tips do not form a cherry.")
    }
    e3 = which(tre$edge[,2]==i3)
    t1 = tre$edge.length[e1]
    t2 = tre$edge.length[e2]
    if (length(e3)==0){
      if (is.null(tre$root.edge)){ t3 = 0 }
      else {t3 = tre$root.edge}
    } else {t3 = tre$edge.length[e3]}
    
    u = t1+t2
    us = sqrt(u)
    D <<- cbind(D, (F[,i1]-F[,i2])/us)
    B <<- cbind(B, G[,i1]*t1/us - G[,i2]*t2/us)
    F[,i3] <<- F[,i1]*t2/u  + F[,i2]*t1/u
    G[,i3] <<- G[,i1]  + G[,i2]
    if (length(e3)>0){
      tre$edge.length[e3] <<- t3 + 1/(1/t1+1/t2)
    } else {
      tre$root.edge       <<- t3 + 1/(1/t1+1/t2)
    }
    tips <<- c(setdiff(tips, c(i1,i2)), i3)
    
  }
  
  laststep = function(iroot){
    D <<- cbind(D, F[,iroot]/sqrt(tre$root.edge))
    B <<- cbind(B, G[,iroot]*sqrt(tre$root.edge))
    tips <<- setdiff(tips,iroot)
  }
  

  F=NULL; G=NULL; D=NULL; B=NULL; tips=NULL;
  if ( alpha > 0){
    tre <- transf.branch.lengths(tre0, model = root.model, parameters = list(alpha=alpha))$tree;
  } else{
    tre = tre0;
  }
  
  initialize(tre);
  
  #repeat{
  #  if( length(tips) < 2 ){ break;}

  #  ancNode = sapply(tips, function(t){ tre$edge[ which(tre$edge[,2]==t), 1] } );
  #  a  = ancNode[duplicated(ancNode)][[1]];

  #  tp      = c(1,2);
  #  tp[[1]] = tre$edge[ which( tre$edge[,1] == a )[[1]], ][[2]];
  #  tp[[2]] = tre$edge[ which( tre$edge[,1] == a )[[2]], ][[2]];
  #  
  #  onestep(tp[[1]],tp[[2]]);
  #}

  for(i in seq(1,length(tre$edge.length),2) ){
    if( length(tips) < 2 ){ break;}

    tp      = c(1,2);
    tp[[1]] = tre$edge[ i,   ][[2]];
    tp[[2]] = tre$edge[ i+1, ][[2]];
    
    onestep(tp[[1]],tp[[2]]);
  }




  laststep( tips[[1]] )
  return( list( F=F, G=G, D=D, B=B) );
}
