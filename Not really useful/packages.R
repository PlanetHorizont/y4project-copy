
install.packages("TreeSim")
library(TreeSim)
trees<-sim.bd.taxa.age(10,2,2,0.5,0.6,2,mrca=FALSE)
plot(trees[[1]])
LTT.plot(trees)

install.packages('ctv')
library('ctv')
install.views('Phylogenetics')
update.views('Phylogenetics')
library(ape)
tree<-rtree(n=20)
plot(tree,edge.width=2)
tree
str(tree)
tree$tip.label
tree$edge.length
tree$Nnode
tree<-read.tree(text="(((A,B),(C,D)),E);")
plot(tree,type="cladogram",edge.width = 2)


install.packages("coalesceR")

n<-5
T<-c(0)
for ( i  in 2:n ) {
  t<-rexp(1,choose(i,2))
  T<-c(T,t)
}
T

A<-seq(1:n)
if (length(A)>1) {
  sample(A)
}
