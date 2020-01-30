
##--------------Estimation with Penalty by BIC----------------------##
mvr_bic <- function(Y,X,Z,method,lambda,opts,opts_pen){    
  p = opts$p
  q = opts$q
  n = opts$n
  pz = opts$pz
  nlam = opts_pen$nlam
  if(opts_pen$isPenColumn){
    fit = EstMVR_colwise(Y,X,Z,lambda,opts,opts_pen)
    df = q*(colSums(fit$df)+pz)
    loglikelih =  n*q * log(fit$likhd/(n*q))
    bic <- switch (method,
                   BIC = loglikelih + log(n*q)*df,
                   AIC = loglikelih + 2*df,
                   GCV = fit$likhd*(n*q)/(n*q-df)^2,
                   EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p-1)/2+1) 
                                                        - lgamma(df+1) - lgamma(q*p*(p-1)/2-df+1))
    )
    selected = which.min(bic)
    lambda_opt = lambda[selected]
    likhd_opt = fit$likhd[selected]
    df_opt = fit$df[,selected]
    Bhat = matrix(fit$betapath[,selected],p)
    if(pz) Chat = matrix(fit$Cpath[,selected],pz)
    else Chat = NULL
  }
  else{
    lambda_opt = rep(0,q)
    likhd_opt = rep(0,q)
    selected = rep(0,q)
    df_opt = matrix(0,q,p)
    mu_opt = rep(0,q)
    bic = matrix(0,q,nlam)
    Bhat = matrix(0,p,q)
    if(pz) Chat = matrix(0,pz,q)
    else Chat = NULL
    fit = EstMVR_lasso(Y,X,Z,lambda,opts,opts_pen)
    for(j in 1:q){
      df0 = fit$df[((j-1)*p+1):(j*p),]
      df = colSums(df0) + pz
      likhd = fit$likhd[j,]
      loglikelih =  n*log(likhd/n)
      bic0 <- switch (method,
                     BIC = loglikelih + log(n)*df,
                     AIC = loglikelih + 2*df,
                     GCV = likhd*n/(n-df)^2,
                     EBIC = loglikelih + log(n)*df + 2*(lgamma(p*(p-1)/2+1) 
                                                          - lgamma(df+1) - lgamma(p*(p-1)/2-df+1))
      )
      bic[j,] = bic0
      selected1 = which.min(bic0)
      lambda_opt[j] = lambda[selected1,j]
      likhd_opt[j] = likhd[selected1]
      df_opt[j,] = df0[,selected1]
      Bhat[,j] = fit$betapath[((j-1)*p+1):(j*p),selected1]
      if(pz) Chat[,j] = fit$Cpath[((j-1)*pz+1):(j*pz),selected1]
      selected[j] = selected1
    }
  }
  return(list(rss = likhd_opt,
              activeX = df_opt,
              lambda = lambda,
              Bhat = Bhat,
              Chat = Chat,
              selectedID = selected,
              lambda_opt=lambda_opt,
              bic = bic,
              Y = Y,
              X = X
              )
         )
}
