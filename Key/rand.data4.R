#This function simulate mutations on a random coalescent tree with given number of tips.
#Input:
#n= : number of tips in the coalescent tree
#plot.tree=TRUE: plot the tree
#theta= : the parameter which controls the frequency of mutations
#s= : the number of sites
library(ape)

rand.data4<-function(n=10,theta=5,s=5,plot.tree=TRUE){
  tree<-rcoal(n)
  if (plot.tree==TRUE){
    plot(tree,show.tip.label = FALSE,edge.width = 2,
         direction='downwards',root.edge=TRUE,label.offset = 0.1)
    title('Random Coalescent Tree')
  }
  
  k<-length(tree$edge.length)    #number of edges
  Mtimes<-c()
  Mnumber<-c()
  
  for (i in 1:k) {
    l<-tree$edge.length[i]
    T<-c()
    A<-c()
    while (sum(T)<l) {
      A<-T
      v<-rexp(1,theta/2)
      T<-c(T,v)
      if (sum(T)>l) {
        Mtimes[i]<-list(A)
        Mnumber[i]<-length(A)}
    }
  }
#  Mtimes     #On each edge, how long is the time until next mutation occur?
#  Mnumber     #On each edge, how many mutations occur?
  S=sum(Mnumber)                #total number of mutations
  m<-tree$Nnode    #number of internal nodes
  edgeMut<-which(Mnumber>0)  #which edges have mutations?
  
  
  states<-c('A','G','C','T')
  t0<-sample(states,s,replace=TRUE); t0  #the sequence of the ancestor
  site.mut<-c()
  for ( i in 1:S){
    site.mut[i]<-sample(1:s,1)
  }
  site.mut     #on which site does the mutation occur?
  
  a<-site.mut
  mut.on.edge<-c()
  for (i in 1:k) {
    c<-Mnumber[i]
    if (c>0) {
      mut.on.edge[i]<-list(a[1:c])
      a<-a[-(1:c)]
    }else{
      mut.on.edge[i]<-NULL
    }
  }
  mut.on.edge
  
  Mutations<-matrix(NA,nrow=S,ncol=7)
  Mutations[,1]<-c(1:S)
  Mutations[,2]<-site.mut
  edge<-c()
  Mnumber.2<-Mnumber[-which(Mnumber==0)]
  for (i in 1:length(Mnumber.2)){
    edge<-c(edge,edgeMut[i])
    while(Mnumber.2[i]>1){
      edge<-c(edge,edgeMut[i])
      Mnumber.2[i]<-Mnumber.2[i]-1
    }
  }
  Mutations[,3]<-edge
  Mutations
  
  #RUN clade!!!
  #RUN clade!!!
  #RUN clade!!!
  #Make sure your have run clade!!!
  Matrix<-matrix(NA,nrow=(n+m),ncol=s)
  for (i in 1:nrow(Matrix)){
    Matrix[i,]<-t0
  }
  original<-c()
  new<-c()
  for (i in 1:S){
    site1<-Mutations[i,2]   #the site on which the mutation takes place
    edge1<-Mutations[i,3]   #the edge on which the mutation takes place
    des1<-clade2(edge1)     #the carriers of this mutation including tips and internal nodes
    anc1<-tree$edge[edge1,1]  #the immediate ancester just before the mutation
    same.edge<-which(edge==edge1)   #which mutations occurs on the same edge as this one
    if (length(which(site.mut[min(same.edge):(i-1)]==site1))==0){  #if no preceding mutations on this edge occurs on the same site as this one does
      original1<-Matrix[anc1,site1]   #then the original state of the site is the same as the state of site of the immediate ancestor
    } else {   #otherwise, the preceding mutations had already change the state of the site for the immediate descendant
      son1<-tree$edge[edge1,2]   #identify the immediate descendant
      original1<-Matrix[son1,site1]   #the origianl state of the site is the state of site of the immediate decendant
    }
    
    original[i]<-original1   #the origianl state of the site on which the mutation occurs
    new1<-sample(states[-which(states==original1)],1)  #choose a new state randomly, besides the original one
    new[i]<-new1    #the new state
    Matrix[des1,site1]<-new1   #replace the state of this site for all decendants including tips and nodes
  }
  
  Matrix   #the matrix gives all sequences of nodes, including the internal nodes (n+1) to (n+m) and tips 1 to n.
  
  alleles<-c()
  descendants<-c()
  for (i in 1:S){
    site1<-Mutations[i,2]    #the site on which the mutation takes place
    edge1<-Mutations[i,3]    #the edge on which the mutation takes place
    des2<-clade(edge1)       #the carriers of this mutation including only tips, no internal nodes
    descendants[i]<-length(des2)
    anc1<-tree$edge[edge1,1]  #the immediate ancester just before the mutation
    nt<-tree$edge[,2]      #all entries of the second column of "edge"
    last.edge<-which((nt)==max(des2))    #according to the order in "edge", what is the index of the last edge that is affected by the mutation
    suc.mutations<-which(edge[(i+1):length(edge)]<(last.edge+1))+i  #which succeeding mutations may revert this mutation
    rev.mutations<-suc.mutations[which(site.mut[suc.mutations]==site1)]  #succeeding mutations which occurs on the same site (i.e. actually revert the mutation)
    if (length(rev.mutations)==0){ #if no succeeding mutations on this edge occurs on the same site as this one does
      alleles[i]<-length(des2)  #the number of tips that carry this mutation 
    } else {
      des3<-c()
      for (j in 1:length(rev.mutations)){
        index<-rev.mutations[j]
        edge2<-Mutations[index,3]
        des3<-c(des3,clade(edge2))
      }
      des4<-unique(des3)
      alleles[i]<-length(des2)-length(des4)
    }
  }
  descendants
  alleles   #the number of tips that actually carries the mutation(not reverted by other mutations)
  
  Mutations[,4]<-original
  Mutations[,5]<-new
  Mutations[,6]<-descendants
  Mutations[,7]<-alleles
  colnames(Mutations)<-c('order','site','edge','original','new','descendants','alleles')
  Mutations
  Tip.matrix<-Matrix[1:n,]   #This matrix only gives all sequences of tips.
  Tip.matrix
  
  plot(descendants,pch=15,ylim=c(0,max(descendants)))
  points(alleles,col='red',pch=16)
  
  
}