
Pi.n = function(u,x,tau)
{
  weights=tau-ifelse(u<=tau,1,0)
  tmp=x*weights
  Pi.value=max(abs(apply(tmp,2,sum)))
  return(Pi.value)
}

penalty = function(x,tau,level, sim.n)
{
  n=dim(x)[1];
  U=matrix(runif(n*sim.n),n,sim.n) 
  values = apply(U,2,Pi.n,x=x,tau=tau)
  values = sort(values)
  return(values[ceiling(level*sim.n)])
}

ini.est = function(x,y,index,tau,level=0.85,lam=NULL,sim.n=1000){
  if(is.null(lam)){
    lambda=penalty(x[,-index],tau,level=level,sim.n=sim.n)
  }else{lambda=lam}
  p=dim(x)[2]
  lam.vec=numeric(p)
  lam.vec[-index]=lambda
  result = rq.fit.lasso(x,y,tau=tau,lambda=lam.vec)
  return(list(residual=result$residual,coefficient=result$coefficient))
}

eff.est.glasso = function(z,x,family='gaussian',nlam=100,crit="BIC"){
  fit = grpreg(X=z,y=x,penalty="grLasso",family=family,nlambda=nlam)
  if(crit=="AIC"){
    tmp=2*fit$loss+2*fit$df
  }else{tmp=2*fit$loss+log(fit$n)*fit$df}
  ind=which(tmp==min(tmp))
  return(fit$beta[,,ind])
}

eff.est = function(z,x,family='gaussian',nlam=100,crit="BIC"){
  d = ncol(x)
  p = ncol(z)
  H_hat = NULL
  for(k in 1:d){
    cvfit <- cv.glmnet(z,x[,k],family="gaussian")
    coe=coef(cvfit,s="lambda.min")
    H_hat = cbind(H_hat,coe)
  }
  return(t(H_hat))
}

my.est = function(y,x,z,tau,method,pen,eps,sim.level)
{
  d = ncol(x)
  result1 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau,level=sim.level)
  if(pen=="glasso")
    H = eff.est.glasso(z,x)
  else
    H=eff.est(z,x)
  
  if(method=='Iterative'){
    itnum=iter.num
    m1=dim(x)[2]
    m2=dim(z)[2]
    beta1=result1$coef[1:m1];beta2=beta1-1
    y.tilde=y-cbind(1,z)%*%result1$coef[(m1+1):(m1+m2)]
    x.tilde=x-cbind(1,z)%*%t(H)
    y.new=y.tilde
    while(i<itnum&max(abs(beta1-beta2))>eps){
      beta2=beta1
      result2 = rq(x=x.tilde,y=y.new,tau=tau)
      beta1=result2$coef
      y.new=y.tilde-(cbind(1,z)%*%t(H))%*%beta1
    }
  }else{
    x.tilde=x-cbind(1,z)%*%t(H)
    y.tilde=y-cbind(cbind(1,z)%*%t(H),cbind(1,z))%*%result1$coef
    result2=rq.fit(x=x.tilde,y=y.tilde,tau=tau)
    beta1=result2$coef
  }
  return(result2)
}
inferen = function(y,x,z,tau,method="OneStep",pen="glasso",eps=1e-6,sim.level=0.85,iter.num=100,RCV=F,K=1,weights=NULL,B=100)
{
  result=my.est(y,x,z,tau,method,pen,eps,sim.level=sim.level)
  if(RCV){
    d=dim(x)[2]
    q=dim(z)[2]
    n=dim(x)[1]
    Cov=matrix(0,d,d)
    for(i in 1:K){
      ind=sort(sample(n,as.integer(n/2)))
      x0=x[ind,];x1=x[-ind,]
      y0=y[ind];y1=y[-ind,]
      z0=z[ind,];z1=z[-ind,]
      result0 = ini.est(cbind(x0,1,z0),y0,index=c(1:(d+1)),tau,level=sim.level/sqrt(2))
      ind0=abs(result0$coef[-c(1:(d+1))])>eps
      z1.0=z1[,ind0];
      
      result1=rq.fit(x=cbind(x1,1,z1.0),y=y1,tau=tau)
      betas=boot.wild(x=x1,z=z1,result=result1,tau=tau,pen=pen,weights=weights,B=B,sim.level=sim.level/sqrt(2),method,eps)
      Cov0=cov(t(betas))
      result1 = ini.est(cbind(x1,1,z1),y1,index=c(1:(d+1)),tau,level=sim.level/sqrt(2))
      ind1=abs(result1$coef[-c(1:(d+1))])>eps
      z0.1=z0[,ind1];
      
      result0=rq.fit(x=cbind(x0,1,z0.1),y=y0,tau=tau)
      betas=boot.wild(x=x0,z=z0,result=result0,tau=tau,pen=pen,weights=weights,B=B,sim.level=sim.level/sqrt(2),method=method,eps=eps)
      Cov1=cov(t(betas))
      Cov=Cov+(Cov0+Cov1)/2
    }
    Cov=Cov/K
  }
  else{
    result0 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau,level=sim.level)
    ind0=abs(result0$coef[-c(1:(d+1))])>eps
    z1=z[,ind0]
    result1=rq.fit(x=cbind(x,1,z1),y=y,tau=tau)
    betas=boot.wild(x=x,z=z,result=result1,tau=tau,pen=pen,weights=weights,B=B,sim.level=sim.level,method=method,eps=eps)
    Cov=cov(t(betas))
  }
  return(list(est=result,cov=Cov))
}

