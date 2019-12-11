
score=function(tau,x){
  tau-ifelse(x<0,1,0)
}

resamples=function(w,resids,fit){
   return(fit+w*abs(resids))
}

boot.wild=function(x,z,y=NULL, result=NULL,tau=0.5,weights=NULL,B=1000,sim.level,method,eps)
{
    n = nrow(x)
    if(tau>=1|tau<=0)
      stop('tau is out of range!')
    else{
      if(is.null(weights)){
       boot.n=B
       w.boot=matrix(sample(c(-2*tau,2*(1-tau)),size=n*boot.n,prob=c(tau,1-tau),replace=T),n,boot.n)
      }else{
       boot.n=dim(weights)[2]
       w.boot=weights
      }
      if(is.null(result)&is.null(y))
        stop("No responses or residuals!")
      else if(is.null(result))
        result=rq.fit(x,y,tau=tau)
      w=diag(x%*%solve(t(x)%*%x)%*%t(x))
      residuals.b=result$residual
      y.boot=apply(w.boot,2,resamples,resids=residuals.b,fit=result$fitted)
      temp=apply(y.boot,2,my.est,x=x,z=z,tau=tau,sim.level=sim.level,method=method,eps=eps)
      coeffs=sapply(temp,coef)
      return(coeffs)
    }
}