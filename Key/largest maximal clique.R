#This function takes an adjacency matrix as an input and works out the largest maximal cliques for the induced graph.

library(igraph)

largest.maximal.clique<-function(X){
  g<-graph_from_adjacency_matrix(X,mode="undirected",diag=F)
  a<-largest.cliques(g)
  clique1<-a[[1]]
  return(clique1)
}

largest.maximal.clique(CM)

plot(graph_from_adjacency_matrix(CM,mode="undirected",diag=F))
#useful functions:
#(g is graph)
maximal.cliques(g)
largest_cliques(g)
plot(g)
g2<-induced.subgraph(graph=g,vids=clique1)
plot(g2)
