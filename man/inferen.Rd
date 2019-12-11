\name{inferen}
\alias{inferen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provide CI of individual coefficient of high-dimensional quantile regression.
}
\description{
This function provides the confidence intevals of individual coefficient of high-dimensional quantile regression by a regularized projection score method for estimating treatment effects. One-step estimation procedure can speed up computation, and the Bootstrap method can narrow the length of CI. 
}
\usage{
inferen(y,x,z,tau,method="OneStep",pen="glasso",eps=1e-6,level=0.85,iter.num=100,RCV=F,K=1,weights=NULL,B=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response, a vector of size \eqn{n}
}
  \item{x}{
The treatment effects, a matrix with dimension \eqn{n\times p} 
}
  \item{z}{
The confounders a matrix with dimension \eqn{n\times q}
}
  \item{tau}{
The given quantile, a scale in the unit inteval
}
  \item{method}{
The method including "OneStep", "Iterative". "OneStep" denotes one-step method ( Feng et al. 2019); "Iterative" denotes that the iteration stops when algorithm conveges. Default  is "OneStep".
  }
  \item{pen}{
The penalty including "glasso" and "lasso". "glasso" denotes the grouped lasso that is used in the regression of 
treatment effect on confounders; "lasso" denotes the lasso. Default is "glasso".
}
  \item{eps}{
The perturbation when the proposed algorithm is used. Default is epsilon=1e-6.
}
  \item{level}{
The length of tuning parameter \eqn{\alpha} which is selected automatically. Default is 50.
}
  \item{iter.num}{
The number of folds for the tuning selection by CV. Default is 5.
}
  \item{RCV}{
Use refitted cross validation method and wild bootstrap to estimate the asymptotic covariance matrix. Default is False.
}
  \item{K}{
The number of repeated RCV. Default is 1.
}
  \item{weights}{
The weights used for wild bootstrap; if not specified (=NULL). Default is NULL.
}
  \item{B}{
The size for bootstrap. Default is 1000.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{est}{Estimator of \eqn{\beta}. It is a list.}
  
  \item{cov}{Covariance matrix of \eqn{\beta}. It is a \eqn{d\times d}-matrix.}
}

\references{
Feng, X., Huang, J. and Liu, X. (2019). Regularized projection score estimation of treatment effects in high-dimensional quantile regression. Manuscript.
}
\author{
Feng Xingdong, Huang, Jian and Xu Liu (liu.xu@sufe.edu.cn)
}


\examples{
n <- 50
d <- 3
s <- 3
p <- 20
beta <- rep(3,d)
eta <- c(rep(3,s),numeric(p-s))
x <- matrix(rnorm(n*d),n,d)
z <- matrix(rnorm(n*(p-1)),n,p-1)
y <- x\%*\%beta + cbind(1,z)\%*\%eta + rnorm(n)
fit <- inferen(y,x,z,tau=0.5)
ests <- fit$result
est.coef <- ests$coef
boot.var <- diag(fit$cov)
lbounds <- ests$coef - qnorm((1+alpha)/2)*sqrt(boot.var)
ubounds <- ests$coef + qnorm((1+alpha)/2)*sqrt(boot.var)
counts <- ifelse(lbounds<beta&beta<ubounds,1,0)
coverage <- cbind(coverage,counts)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Quantile regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ Projection score }% __ONLY ONE__ keyword per line
