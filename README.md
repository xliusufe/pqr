# pqr
  An R package to construct confidence intervals of coefficients via regularized projection score estimation of treatment effects 
  in high-dimensional quantile regression.
 
  A regularized projection score method is proposed for estimating treatment effects in quantile regression 
  in the presence of high-dimensional confounding covariates. This method is based on an estimated projection 
  score function of the low-dimensional treatment parameters in the presence of high-dimensional confounding 
  covariates. We propose one-step algorithm and a reffitted wild bootstrapping approach for variance estimation. 
  This enables us to construct confidence intervals for the treatment effects in the high-dimensional circumstances.
  
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/pqr")

# Usage

   - [x] [pqr-manual](https://github.com/xliusufe/pqr/blob/master/inst/pqr-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(pqr)
	
    # Example 1
	
    n <- 50
	d <- 3
	s <- 3
	p <- 20
	alpha <- 0.95
	beta <- rep(3,d)
	eta <- c(rep(3,s),numeric(p-s))
	x <- matrix(rnorm(n*d),n,d)
	z <- matrix(rnorm(n*(p-1)),n,p-1)
	y <- x%*%beta + cbind(1,z)%*%eta + rnorm(n)
	fit <- inferen(y,x,z,tau=0.5)
	ests <- fit$ests
    est.coef <- ests$coef
	boot.var <- diag(fit$cov)
    lbounds <- ests$coef - qnorm((1+alpha)/2)*sqrt(boot.var)
    ubounds <- ests$coef + qnorm((1+alpha)/2)*sqrt(boot.var)
    counts <- ifelse(lbounds<beta&beta<ubounds,1,0)
	
	# Example 2
	
    n <- 100
	q <- 5
	s <- 3
	p <- 100
	B <- matrix(runif(q*s, 2,3), s)
	x <- matrix(rnorm(n*p),n,p)
	y <- x[,1:s]%*%B + rnorm(n)
	fit <- mvr(y,x)
	fit$activeX
	fit$Bhat
	which(rowSums(fit$Bhat^2)>0)
	fit$muhat 

# References
 
Cheng, C., Feng, X., Huang, J. and Liu, X. (2020). Regularized projection score estimation of treatment effects 
in high-dimensional quantile regression. Manuscript.

# Development

The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn) and Chao Cheng.
