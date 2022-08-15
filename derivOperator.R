derivOperator<-function(PL,gamma,h,x,z,p){
#Construct a shape constraint matrix at a corresponding sequence of x data and sequence of gridpoint z. \\
#Output: A differential operator corresponding to the penalty. \\
#Input: Penalty label $PL$, which is the type of shape constraint, can be 
#\item \verb!R1! (first derivative constraint), \item \verb!Roughness! (second derivative constraint), 
#\item \verb!Exponential! (global exponential shape constraint), 
#\item \verb!Periodicity! (global periodic shape constraint, use the second and fourth derivative),
#\item \verb!Periodicity2! (global periodic shape constraint, use the first and second derivative);
#the shape constraint parameter $gamma$; the bandwidth $h$; numeric vector of $x$ data; numeric vector of gridpoint $z$ data; degree of local polynomial $p$. \\
  B<-matrix(0,ncol=length(z),nrow=length(x))
  for(j in 1:length(z))
    B[,j]<-getbb(PL,gamma,h,x,z[j],p)
  B
}
getbb<-function(PL,gamma,h, xx,x,p){
  aaFunt<-function(x) {getaa(h,xx,x,p)}
  bb<-rep(0,length(xx))
  if(PL=="R1"){
    bb<-numericalDerivative(x,aaFunt,k=1)
  }
  if(PL=="Roughness"){
    bb<-numericalDerivative(x,aaFunt,k=2)
  }
  if(PL=="Exponential"){
    bb<-numericalDerivative(x,aaFunt,k=2)+gamma*numericalDerivative(x,aaFunt,k=1)
  }
  if(PL=="Periodicity"){
    bb<-numericalDerivative(x,aaFunt,k=4)+gamma*numericalDerivative(x,aaFunt,k=2)
  }
  if(PL=="Periodicity2"){
    bb<-numericalDerivative(x,aaFunt,k=2)+gamma*numericalDerivative(x,aaFunt,k=0)
  }
  bb
}
