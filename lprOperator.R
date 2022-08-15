lprOperator<- function(h,x,z,p){
# Construct a matrix based on local polynomial regression at a corresponding sequence of x data and sequence of gridpoint z. \\
# Output: A local polynomial regression operator. \\
# Input: The bandwidth $h$; numeric vector of $x$ data; numeric vector of gridpoint $z$ data; degree of local polynomial $p$. \\
  A<-matrix(0,ncol=length(z),nrow=length(x))
  for(j in 1:length(z))
    A[,j]<-getaa(h,x,z[j],p)
  A
}

getX<-function(p,xx,x){	
  X<-matrix(0,ncol=p+1,nrow=length(xx))
  for(k in 1:(p+1)) X[,k]<-(xx-x)^(k-1)
  X
}

getWi<-function(kernel,h,xx,x){ 
  xxtmp<-(xx-x)/h
  Wi<-ifelse(abs(xxtmp)<0.0001,1,0)
  if(kernel=="Gaussian")
    Wi <-dnorm(xxtmp, sd=1)
  if(kernel=="Uniform")
    Wi <-ifelse(abs(xxtmp)<1,1,0)
  if(kernel=="Triangular")
    Wi <- ifelse(abs(xxtmp)<1,(1-abs(xxtmp)),0)
  if(kernel=="Tricube")
    Wi <-ifelse(abs(xxtmp)<1,(1-abs(xxtmp)^3)^3,0)
  Wi
}

getaa<-function(h,xx,x,p){
  X<-getX(p,xx,x)
  W<-diag(getWi(kernel="Gaussian",h,xx,x))
  aa<-(solve(t(X)%*%W%*%X)%*%t(X)%*%W)[1,]
  aa
}


