# The argument X is a matrix with species in rows and characters in columns
# compatibility.matrix gives the compatibility matrix of X, with 1 indicating the two characters are compatible and 0 indicating they are not compatible.

compatibility.matrix<-function(X,names=TRUE) {
  cm<-matrix(NA,nrow=ncol(X),ncol=ncol(X))
  for ( i in 1:nrow(cm)) {
    for (j in 1:ncol(cm)) {
      if (compatibility.test(X[,i],X[,j])==TRUE) {
        cm[i,j]<-1
      } else { cm[i,j]<-0}
   }
  }
  if (names==TRUE) {
    colnames(cm)<-colnames(X)
    row.names(cm)<-colnames(X)
  }
  cm
}

CM<-compatibility.matrix(X)
