rm(list=ls())

#simulate the times in the coalescent tree
n<-10
T<-c()
T[1]<-0
for (k in 2:n) {
  T[k]<-rexp(1,rate=choose(k,2))
}
T

#Example 1

c1<-c(1,0,1,1,0,0)
c2<-c(0,0,1,1,0,0)
c3<-c(0,1,0,0,1,0)
c4<-c(1,0,0,1,1,0)
c5<-c(1,0,0,1,1,0)
c6<-c(0,0,0,1,0,0)
Y<-data.frame(c1,c2,c3,c4,c5,c6)
rownames(Y)<-c("Alpha","Beta","Gamma","Delta","Epsilon","Omega")

X<-matrix(data=c(c1,c2,c3,c4,c5,c6),nrow=6,ncol=6)
colnames(X)<-c("c1","c2","c3","c4","c5","c6")
row.names(X)<-c("Alpha","Beta","Gamma","Delta","Epsilon","Omega")
X
dim(X)
ncol(X)

CM<-compatibility.matrix(X)
lmc<-largest.maximal.clique(CM)  #vertices of the largest maximal clique
X[,c(data.matrix(lmc)[,1])]   #a sub-matrix of original data that only includes vertices of the largest maximal clique



i<-1

   A<-list()
     v1<-c()
     v2<-c()
  for (j in 1:nrow(Y)) {

     if (Y[j,i]==0) {
        v1<-c(v1,rownames(Y)[j])
     } 
     if (Y[j,i]==1) {
        v2<-c(v2,rownames(Y)[j])
     }   
  }
 A<-c(A,list(v1,v2))
A



#generate an object of class "phylo"
tr <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "a", Nnode = 1L)  #simplest tree
class(tr) <- "phylo"
str(tr)


