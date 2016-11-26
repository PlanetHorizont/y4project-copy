A<-matrix(0,nrow=1000,ncol=6)
A[1,1]<-1
for ( i in 1:999) {
  U1<-runif(1)
  U2<-runif(1)
  A[i,2]=-log(U1)/(choose(A[i,1]+1,2)+2.5*(A[i,1]+1)+1.5*(A[i,1]+1))
  A[i,3]=U2<(choose(A[i,1]+1,2)+2.5*(A[i,1]+1))/(choose(A[i,1]+1,2)+2.5*(A[i,1]+1)+1.5*(A[i,1]+1))
  A[i,4]=U2<(2.5*(A[i,1]+1))/(choose(A[i,1]+1,2)+2.5*(A[i,1]+1)+1.5*(A[i,1]+1))
  A[i,5]=-1+A[i,3]+A[i,4]
  A[i,6]=floor(runif(1)*(A[i,1]+A[i,5]))
  A[i+1,1]=A[i,1]+A[i,5]
}
A

i<-1
