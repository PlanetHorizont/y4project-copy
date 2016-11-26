carriers<-function(a) {
  anc <- tree$edge[a, 1]
  son <- tree$edge[a, 2]
  ind <- which(tree$edge[, 1] == anc)  #indices of edges whose direct ancester is anc
  bro <- which(ind == a) + 1
  
  while(bro > length(ind)){
    b<-which(tree$edge[,2]==anc)
    if(length(b)==0){     #if anc is the root of the tree
      ind <- which(tree$edge[, 1] == anc)    
      bro=2
    }else{
      anc<-tree$edge[b, 1]
      ind <- which(tree$edge[, 1] == anc)
      bro <- which(ind == b) + 1
    }
  }
  c<-ind[bro]
  if (a+2 > c) {
    clade<-tree$edge[a,2]
  } else {
    clade<-tree$edge[(a+1):(c-1), ]
  }
  des<-unique(sort(clade))
  return(des)
}

