random.sample<-sample(c(0,1), 200, replace = TRUE)
X<-matrix(random.sample,nrow=10,ncol=20)
colnames(X)<-c(paste(c('c'),c(1:ncol(X)),sep=''))  #20 characters
row.names(X)<-c(paste(c('s'),c(1:nrow(X)),sep=''))  #10 species
X
CM<-compatibility.matrix(X)
CM
largest.maximal.clique(CM)


