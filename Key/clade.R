#Finding a Clade
# a is the index of an edge with mutation(s)

#This function gives all tips who carries the mutation(s)
clade<-function(a){
  son<-tree$edge[a, 2]
  if (son<(n+1)){
    des<-son        #the direct descendant of the edge is a tip
  } else {              ##the direct descendant of the edge is an internal node
    ind<-son-n
    sub<-subtrees(tree)[[ind]]
    des<-c()
    for ( i in 1:sub$Ntip){
      des<-c(des,which(tree$tip.label==sub$tip.label[i]))
    }
  }
  return(des)
}


# This function gives all the carriers of the mutation(s), including tips and internal nodes
clade2<-function(a){
  son<-tree$edge[a, 2]
  if (son<(n+1)){
    des<-son        #the direct descendant of the edge is a tip
  } else {              ##the direct descendant of the edge is an internal node
    ind<-son-n
    sub<-subtrees(tree)[[ind]]
    des<-c()
    for ( i in 1:sub$Ntip){
      des<-c(des,which(tree$tip.label==sub$tip.label[i]))
    }
    des<-c(des,sub$node.label)
  }
  return(des)
}

