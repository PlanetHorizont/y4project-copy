library(ape)

#simulate a tree
tree<-rcoal(20)
plot(tree)
str(tree)
tree$edge
tree$edge.length


#simulate mutations
theta<-5

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
S=sum(Mnumber);S                 #total number of mutations
n<-length(tree$tip.label); n #number of genes
m<-tree$Nnode;m     #number of internal nodes

meanS<-function(n){
  s<-0
  for (i in 1:(n-1)) {
    s<-s+1/i
  }
  return(s*theta)
}
meanS(n)    #expectation of the number of mutations on a coalescent tree of n genes

varS<-function(n){
  s<-0
  v<-0
  for (i in 1:(n-1)) {
    s<-s+1/i
    v<-v+1/(i^2)
  }
  return(s*theta+(theta^2)*v)
}
varS(n)     #variation of the number of mutations on a coalescent tree of n genes


sapply(Mtimes,function(x)sum(x))   #compare to
tree$edge.length
sapply(Mtimes,function(x)sum(x))<tree$edge.length  #double check

edgeMut<-which(Mnumber>0)
edgeMut       #which edges have mutation(s)
tree$edge[,2][edgeMut] #the immediate descendant after the mutation


label.mutations<-runif(S)  #label mutations using independent uniform random variables in [0,1]
sort(label.mutations)    
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
mut.on.edge  #this list tells which mutation(s) are on which edge


#Finding a Clade
# a is the index of an edge with mutation(s), and this function gives all the carriers of the mutation(s)

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

Incidence.matrix<-matrix(0,nrow=n,ncol=S)
row.names(Incidence.matrix)<-c(paste(c('tip'),c(1:n),sep=''))
colnames(Incidence.matrix)<-c(paste(c('sites'),c(1:S),sep=''))

for (i in 1:length(edgeMut)){
  a<-edgeMut[i]
  des <- clade(a)
  des.tip<-des[which(des<(n+1))]
  sites <- unlist(mut.on.edge[a])
  Incidence.matrix[des.tip,sites]<-1   
}

Incidence.matrix
#This gives the incidence matrix in which each column representsa segregating site and each row represents a haplotype
#the row names of the matrix corresponds to the relative position of the tips (lower to upper, 1 to n), not the tip labels(t1,t2,etc).

#try compatibility matrix
compatibility.matrix(Incidence.matrix)   #all compatible as expected
