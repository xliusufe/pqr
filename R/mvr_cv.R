
##--------------Estimation with Penalty by CV----------------------##
mvr_cv <- function(Y,X,ncv,lambda,opts,opts_pen){
  p <- opts$p
  q <- opts$q
  n <- opts$n
  nlam <- opts_pen$nlam
  len_cv = ceiling(n/ncv)
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
      nt = nrow(Ytest)
      
      fit = EstMVR_glasso(Ytrain,Xtrain,lambda,opts,opts_pen)
      df = colSums(fit$df)
      for(kk in 1:nlam){
        Bnew = matrix(fit$betapath[,kk],p)
        likhd[kk] = sum((Ytest-Xtest%*%Bnew)^2)
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
      nt = nrow(Ytest)
      
      fit = EstMVR_lasso(Ytrain,Xtrain,lambda,opts,opts_pen)
      for(j in 1:q){
        for(kk in 1:nlam){
          Bnew = fit$betapath[((j-1)*p+1):(j*p),kk]
          likhd[kk] = sum((Ytest[,j]-Xtest%*%Bnew)^2)
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
    fit_opt = EstMVR_glasso(Y,X,lambda_opt,opts,opts_pen)
    activeX = c(1,fit_opt$df)
    Bhat = matrix(fit_opt$betapath,p)
    bic_opt = fit_opt$likhd
  }
  else{
    fit_opt = EstMVR_lasso(Y,X,matrix(lambda_opt,1,q),opts,opts_pen)
    activeX = rbind(1,matrix(fit_opt$df,p))
    Bhat = matrix(fit_opt$betapath,p)
    bic_opt = fit_opt$likhd
  }
  return(list(rss=fit_opt$likhd,
              activeX = t(activeX),
              lambda = lambda,
              selectedID = selected,
              RSS = bic_opt,
              lambda_opt=lambda_opt,
              Bhat = Bhat,
              Y = Y,
              X = X
  )
  )
}