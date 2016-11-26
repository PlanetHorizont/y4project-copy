library(igraph)
D<-read.table(header=TRUE, text=
'from to
A B
A C
C D
C F
C E
D E
D F
E F'
)

g1<-graph.data.frame(D,directed = F)
plot(g1)
a<-largest.cliques(g1)
a
clique1<-a[[1]]
clique1

g2<-induced.subgraph(graph=g1,vids=clique1)
plot(g2)
g2
