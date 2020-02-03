
##--------------Estimation with Penalty by CV----------------------##
mvr_cv <- function(Y,X,Z,ncv,lambda,opts,opts_pen){
  p <- opts$p
  q <- opts$q
  n <- opts$n
  pz = opts$pz
  nlam <- opts_pen$nlam
  len_cv = floor(n/ncv)
  likhd = rep(0,nlam)
  
  if(opts_pen$isPenColumn){
    bic = rep(0,nlam)
    for(jj in 1:ncv){ # start CV
      cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
      if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
      Ytrain = Y[-cv.id,]
      Xtrain = X[-cv.id,]
      Ytest = Y[cv.id,]
      Xtest = X[cv.id,]
      Ztrain = Z[-cv.id,]
      Ztest = Z[cv.id,]
      nt = nrow(Ytest)
      
      fit = EstMVR_colwise(Ytrain,Xtrain,Ztrain,lambda,opts,opts_pen)
      df = colSums(fit$df)
      if(pz){
        for(kk in 1:nlam){
          Bnew = matrix(fit$betapath[,kk],p)
          Cnew = matrix(fit$Cpath[,kk],pz)
          likhd[kk] = sum((Ytest-Ztest%*%Cnew-Xtest%*%Bnew)^2)
        }
      }
      else{
        for(kk in 1:nlam){
          Bnew = matrix(fit$betapath[,kk],p)
          likhd[kk] = sum((Ytest-Xtest%*%Bnew)^2)
        }
      }
      bic = bic + likhd
    }
    selected = which.min(bic)
    lambda_opt = lambda[selected]
  }
  else{
    lambda_opt = rep(0,q)
    selected = rep(0,q)
    bic = matrix(0,q,nlam)
    for(jj in 1:ncv){ # start CV
      cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
      if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
      Ytrain = Y[-cv.id,]
      Xtrain = X[-cv.id,]
      Ytest = Y[cv.id,]
      Xtest = X[cv.id,]
      Ztrain = Z[-cv.id,]
      Ztest = Z[cv.id,]
      nt = nrow(Ytest)
      
      fit = EstMVR_lasso(Ytrain,Xtrain,Ztrain,lambda,opts,opts_pen)
      for(j in 1:q){
        if(pz){
          for(kk in 1:nlam){
            Bnew = fit$betapath[((j-1)*p+1):(j*p),kk]
            Cnew = fit$Cpath[((j-1)*pz+1):(j*pz),kk]
            likhd[kk] = sum((Ytest[,j]-Ztest%*%Cnew-Xtest%*%Bnew)^2)
          }
        }
        else{
          for(kk in 1:nlam){
            Bnew = fit$betapath[((j-1)*p+1):(j*p),kk]
            likhd[kk] = sum((Ytest[,j]-Xtest%*%Bnew)^2)
          }
        }
        bic[j,] = bic[j,] + likhd
      }
    }
    for(j in 1:q){
      selected1 = which.min(bic[j,])
      selected[j] = selected1
      lambda_opt[j] = lambda[selected1,j]
    }
  } # end of CV
  #---------------- The estimation after selection ---------------------#
  if(opts_pen$isPenColumn){
    fit_opt = EstMVR_colwise(Y,X,Z,lambda,opts,opts_pen)
    activeX = fit_opt$df[,selected]
    Bhat = matrix(fit_opt$betapath[,selected],p)
    if(pz) Chat = matrix(fit_opt$Cpath[,selected],pz)
    else Chat = NULL
    likhd_opt = fit_opt$likhd[selected]
  }
  else{
    activeX = matrix(0,q,p)
    Bhat = matrix(0,p,q)
    likhd_opt = rep(0,q)
    if(pz) Chat = matrix(0,pz,q)
    else Chat = NULL
    
    fit_opt = EstMVR_lasso(Y,X,Z,lambda,opts,opts_pen)
    for(j in 1:q){
      activeX[j,] = fit_opt$df[((j-1)*p+1):(j*p),selected[j]]
      Bhat[,j] = fit$betapath[((j-1)*p+1):(j*p),selected[j]]
      if(pz) Chat[,j] = fit$Cpath[((j-1)*pz+1):(j*pz),selected[j]]
      likhd_opt[j] = fit_opt$likhd[j,selected[j]]
    }
  }
  return(list(rss = likhd_opt,
              activeX = activeX,
              lambda = lambda,
              selectedID = selected,
              bic = bic,
              lambda_opt=lambda_opt,
              Bhat = Bhat,
              Chat = Chat,
              Y = Y,
              X = X
  )
  )
}