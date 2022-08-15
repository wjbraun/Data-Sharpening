ProxV<-function(x,zeta,radius=NULL,family){
#Construct the proximal operator at a corresponding sequence of $x$ data. Return the proximal result.
#Output: The proximal result.
#Input: numeric vector of $x$ data; the proximal operator parameter $zeta$; the norm ball radius $radius$; the proximal $family$:
#    \item \verb!norm1! ($l_1$ norm ball). 
#    \item \verb!norm2! ($l_2$ norm ball).
#    \item \verb!$norminf! ($l_{\infty}$ norm ball).
#    \item \verb!rectangle! 
#    \item \verb!nonnegative! (non negative finite space)
  if(family=="norm1"){
    out<-x-prox_ball(zeta,radius,family="norm1",x)
  }
  if(family=="norm2"){
    out<-x-prox_ball(zeta,radius,family="norm2",x)
  }
  if(family=="norminf"){
    out<-x-prox_ball(zeta,radius,family="norminf",x)
  }
  if(family=="nonnegative"){
    out<-x-prox_box(zeta,family="nonnegative",x)
  }
  if(family=="rectangle"){
    out<-x-prox_box(zeta,family="rectangle",x)
  }
  out
}

prox_ball<-function(sigma,radius,family=c("norm2","norm1","norminf"),input,a=rep(0,1),b=rep(0,1)){
  x<-input/sigma
  out<-numeric()
  if(family=="norm2"){
    out<-radius*(x)/max(radius,norm(x,type="f"))
  }
  if (family=="norm1"){
    if(sum(abs(x))<=radius){
      lambda<-0
    }else{
      xs<-sort(abs(x),decreasing=TRUE)
      abssum<-xs[1]-xs[2]
      i<-2
      while(abssum<radius){
        if(i<length(xs)){
          abssum<-abssum+i*(xs[i]-xs[i+1])
          i<-i+1
        }else{
          break
        }
      }
      xss<-xs[1:i]
      lambda<-(sum(xss)-radius)/i
    }
    for (j in 1:length(x)){
      if(abs(x[j])>lambda){
        out[j]<-sign(x[j])*(abs(x[j])-lambda)
      }else{
        out[j]<-0
      }
    }
  }
  if (family=="norminf"){
    out<-x
    if(max(abs(x))<radius){
      out<-x
    }else{
      out[which(x>radius)]<-radius
      out[which(x<(-radius))]<--radius
      
    }
  }
  return(sigma*out)
}

prox_box<-function(sigma,family=c("rectangle","nonnegative"),input,bound=c(-1,1)){
  x<-input/sigma
  out<-numeric()
  if(family=="rectangle"){
    for (j in 1:length(x)){
      if(x[j]<=bound[1]){
        out[j]<-bound[1]
      }else if (x[j]>=bound[2]){
        out[j]<-bound[2]
      }else{
        out[j]<-x[j]
      }
    }

  }
  if(family=="nonnegative"){
    out<-x
    out[which(x<0)]<-0
  }
  return(sigma*out)
}


