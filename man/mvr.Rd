\name{mvr}
\alias{mqr-function}
\alias{mvr}
\docType{package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provide CI of individual coefficient of high-dimensional quantile regression
}
\description{
This function provides the coefficient matrix estimator of high-dimensional multivariate regression (MVR) with penalty \code{lasso} (\code{mcp} or \code{scad}). The tuning parameter is selected by \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}
\usage{
mvr(y, x, method="BIC", ncv=10, penalty="LASSO", isPenColumn=TRUE, 
          lambda=NULL, nlam=50, intercept=TRUE, lam_min=1e-4, eps=1e-6, 
          maxstep=20, gamma_pen=2, dfmax=NULL, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The response, a vector of size \eqn{n} or  a matrix of size \eqn{n\times q}}
  
  \item{x}{The treatment effects, a matrix with dimension \eqn{n\times p} }

  \item{method}{The method to be applied to select parameters.  Either \code{BIC}
                (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is not \code{CV}, \code{ncv} is useless.}  
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default), "SCAD", or "MCP".}
 
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{y} is penalized.  Default is \code{TRUE}. If \code{isPenColumn} is \code{TRUE}, the coefficients associating with
                    \eqn{X_j} that affects simultaneously whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. 
                    If \code{isPenColumn} is \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response
                    \eqn{y_l} is penalized for each \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.} 
        
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}
  
  \item{nlam}{The number of lambda values. Default is 50.}
  
  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}  
   
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-3.}
                 
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-4}.}
              
  \item{maxstep}{Maximum number of iterations. Default is \code{20}.}
                  
  \item{gamma_pen}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients.
               Default is no upper bound.  However, for large data sets,
               computational burden may be heavy for models with a large number of
               nonzero coefficients.}
               
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
               the relative contributions from the LASSO, MCP/SCAD penalty and the ridge,
               or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty,
               while \code{alpha=0} would be equivalent to ridge regression.
               However, \code{alpha=0} is not supported; \code{alpha} may be
               arbitrarily small, but not exactly 0.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
   \item{Bhat}{Estimator of \code{Bhat}.}
   
   \item{rss}{Residual sum of squares (RSS).}
   
   \item{activeX}{The active set of \eqn{X}. It is a \eqn{p} dimensional vector.}
  
   \item{lambda}{The sequence of regularization parameter values in the path.}
   
   \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}   
                    
   \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}   
                    
   \item{bic}{BIC value used to select variables.}  
   
   \item{muhat}{Estimator of intercept \eqn{\mu}.}
                    
   \item{Y}{Response \eqn{Y}.}
  
   \item{X}{Design matrix \eqn{X}.}                    
}

\references{
Cheng, C., Feng, X., Huang, J. and Liu, X. (2020). Regularized projection score estimation of treatment effects in high-dimensional quantile regression. Manuscript.
}
\author{
Chao Cheng, Xingdong Feng, Jian Huang and Xu Liu (liu.xu@sufe.edu.cn)
}


\examples{
	library(pqr)

	n <- 100
	q <- 5
	s <- 3
	p <- 100
	B <- matrix(runif(q*s, 2,3), s)
	x <- matrix(rnorm(n*p),n,p)
	y <- x[,1:s]\%*\%B + rnorm(n)
	fit <- mvr(y,x)
	fit$activeX
	fit$Bhat
	which(rowSums(fit$Bhat^2)>0)
	fit$muhat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Multivariate regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ Variable selection }% __ONLY ONE__ keyword per line
