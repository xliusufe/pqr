mvr <- 
  function(y,x,method="BIC",ncv=10,penalty="LASSO",isPenColumn=TRUE,lambda=NULL,nlam=50,
           intercept=TRUE,lam_min=1e-4,eps=1e-6,maxstep=20,gamma_pen=2,dfmax=NULL,alpha=1){
    n <- dim(y)[1]
    q <- dim(y)[2]
    X2bar = colMeans(x)
    Ybar = colMeans(y)
    X2 = x - matrix(rep(X2bar,each=n),n)
    Y = y - matrix(rep(Ybar,each=n),n)
    p <- ncol(X2)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma_pen <- 3
      pen <- 3
    }  
    if (gamma_pen <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma_pen <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    if (is.null(dfmax)) dfmax = p + 1
    
    opts = list(eps=eps,max_step=maxstep,n=n,p=p,q=q) 
    opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma_pen=gamma_pen,alpha=alpha,dfmax=dfmax,isPenColumn=isPenColumn) 
    if(isPenColumn){
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        lambda = setuplambdaMVR_glasso(Y,X2,nlam,setlam)
      }
      else  opts_pen$nlam = length(lambda)
    }
    else{
      if (is.null(lambda)) {
        if (nlam < 1) stop("nlambda must be at least 1")
        if (n<=p) lam_min = 1e-2
        setlam = c(1,lam_min,alpha,nlam)
        lambda = setuplambdaMVR_lasso(Y,X2,nlam,setlam)
      }
      else  opts_pen$nlam = nrow(lambda)
    }
    #---------------- The selection by CV or BIC  ---------------------#  
    if(method=="CV") fit_mvr = mvr_cv(Y,X2,ncv,lambda,opts,opts_pen)
    else fit_mvr = mvr_bic(Y,X2,method,lambda,opts,opts_pen)
    
    if(intercept)   fit_mvr$muhat = Ybar-t(fit_mvr$Bhat)%*%X2bar
    else fit_mvr$muhat = rep(0,q)
    return(fit_mvr)
  }