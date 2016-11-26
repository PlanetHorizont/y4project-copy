#####################################################
#This function with input n(number of tips) and theta, gives the matrix in which rows represent sequences 
#and columns represent sites whose states are in nucleotides base.


#In the infinite sites model, each mutation affects a new site.
#The number of polymorphic sites = the number of mutations that occured
#Hence, on each sites, at most two states can occur. Otherwise, it violates the infinite sites assumption.

library(ape)

rand.data<-function(n=20,theta=5){
  tree<-rcoal(n)
  plot(tree)
  
  k<-length(tree$edge.length)
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
  Mtimes     #On each edge, how long is the time until next mutation occur?
  Mnumber     #On each edge, how many mutations occur?
  S=sum(Mnumber)                 #total number of mutations      
  m<-tree$Nnode     #number of internal nodes
  edgeMut<-which(Mnumber>0)
  
  label.mutations<-runif(S)  #label mutations using independent uniform random variables in [0,1]
  order.mutations<-order(label.mutations)   #order mutations
  a<-order.mutations
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
  
  Incidence.matrix<-matrix(0,nrow=n,ncol=S)
  for (i in 1:length(edgeMut)){
    a<-edgeMut[i]
    des <- clade(a)
    des.tip<-des[which(des<(n+1))]
    sites <- unlist(mut.on.edge[a])
    Incidence.matrix[des.tip,sites]<-1   
  }
  
  #in terms of nucleotides base
  IM2<-matrix(NA,nrow=nrow(Incidence.matrix),ncol=ncol(Incidence.matrix))
  row.names(IM2)<-c(paste(c('tip'),c(1:n),sep=''))
  colnames(IM2)<-c(paste(c('sites'),c(1:S),sep=''))
  s0<-sample(states,ncol(IM2),replace=TRUE);s0
  for (j in 1:ncol(IM2)) {
    a<-which(states==s0[j])
    k<-sample(states[-a],1)
    for (i in 1:nrow(IM2)){
      if (Incidence.matrix[i,j]==0) {IM2[i,j]<-s0[j]}
      else{ IM2[i,j]<-k }
    }
  }
  return(IM2)
}

data<-rand.data(theta=7,n=5)
df<-data.frame(data)
df 
