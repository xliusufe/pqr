mvr <- 
  function(Y,X,Z=NULL,method="BIC",ncv=10,penalty="LASSO",isPenColumn=TRUE,lambda=NULL,nlam=20,
           intercept=TRUE,lam_min=1e-4,eps=1e-6,max_step=20,gamma_pen=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
   
    X2bar = colMeans(X)
    Ybar = colMeans(Y)
    X2 = X - matrix(rep(X2bar,each=n),n)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    if(is.null(Z)){
      Zbar = 0
      pz = 0
      Z1 = matrix(0,n,2)
    }
    else{
      Zbar = colMeans(Z)
      pz = ncol(Z)
      Z1 = Z - matrix(rep(Zbar,each=n),n)
      L = solve(chol(t(Z1)%*%Z1/n))
      Z1 = Z1%*%L
    }
    p <- ncol(X2)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma_pen <- 3.7
      pen <- 3
    }  
    if (gamma_pen <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma_pen <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    if (is.null(dfmax)) dfmax = p + 1
    
    opts = list(eps=eps,max_step=max_step,n=n,p=p,q=q,pz=pz) 
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma_pen,alpha=alpha,dfmax=dfmax,isPenColumn=isPenColumn) 
    if(isPenColumn){
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        if(pz){
          fitlm = lm(Y1~Z1-1)
          lambda = setuplambdaMVR_colwise(fitlm$residuals,X2,nlam,setlam)
        }
        else lambda = setuplambdaMVR_colwise(Y1,X2,nlam,setlam)
      }
      else  opts_pen$nlam = length(lambda)
    }
    else{
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        if(pz){
          fitlm = lm(Y1~Z1-1)
          lambda = setuplambdaMVR_lasso(fitlm$residuals,X2,nlam,setlam)
        }
        else  lambda = setuplambdaMVR_lasso(Y1,X2,nlam,setlam)
      }
      else  opts_pen$nlam = nrow(lambda)
    }
    #---------------- The selection by CV or BIC  ---------------------#  
    if(method=="CV") fit_mvr = mvr_cv(Y1,X2,Z1,ncv,lambda,opts,opts_pen)
    else fit_mvr = mvr_bic(Y1,X2,Z1,method,lambda,opts,opts_pen)
    if(pz) fit_mvr$Chat = L%*%fit_mvr$Chat
    if(intercept){
      if(pz) fit_mvr$muhat = Ybar-t(fit_mvr$Bhat)%*%X2bar-t(fit_mvr$Chat)%*%Zbar
      else  fit_mvr$muhat = Ybar-t(fit_mvr$Bhat)%*%X2bar
    }
    else fit_mvr$muhat = rep(0,q)
    
    return(fit_mvr)
  }