states<-c('A','G','C','T')
s1<-sample(states,20,replace=TRUE)
data<-data.frame(s1,s2=NA,s3=NA,s4=NA,s5=NA,s6=NA)
data

matrix<-data.matrix(data)
matrix
tmatrix<-t(matrix)
tmatrix
colnames(tmatrix)<-c(paste(c('c'),c(1:ncol(tmatrix)),sep=''))



for ( i in 2:nrow(tmatrix)){
  tmatrix[i,]<-tmatrix[i-1,]
  j<-sample(1:ncol(tmatrix),1)
  a<-tmatrix[i,j]
  x<-c(1:4)[-a]
  tmatrix[i,j]<-sample(x,1)
}

m<-matrix(NA,nrow=20,ncol=20)
rownames(m)<-c(paste(c('s'),c(1:20),sep=''))
colnames(m)<-c(paste(c('c'),c(1:20),sep=''))
m[1,]<-s1
for ( i in 2:nrow(m)){
  m[i,]<-m[i-1,]
  j<-sample(1:ncol(m),1)
  a<-which(states==m[i,j])
  m[i,j]<-sample(states[-a],1)
}
m

 

