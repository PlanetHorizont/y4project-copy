# this function tells whether two characters are compatible or not. 
#It returns TRUE for compatible characters and FALSE for incompatible characters.

compatibility.test<-function(a,b){
  result=c()
  ab<-matrix(data=c(a,b),nrow=length(a),ncol=2)
  c<-matrix(data=0,nrow=2,ncol=2)
  for (i in 1:length(a)) {
    if (identical(ab[i,],c(0,0))==T) { c[1,1]<-1 }
    if (identical(ab[i,],c(0,1))==T) { c[1,2]<-1 }
    if (identical(ab[i,],c(1,0))==T) { c[2,1]<-1 }
    if (identical(ab[i,],c(1,1))==T) { c[2,2]<-1 }
  }
  result=(sum(c)<4)
  return(result)
}


a<-c(1,1,0,0)
b<-c(0,1,0,1)
compatibility.test(a,b)


a<-c(1,1,0,0)
b<-c(1,1,0,0)
compatibility.test(a,b)

