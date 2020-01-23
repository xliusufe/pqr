//[[Rcpp::depends(RcppEigen)]]
#include "mvr.h"

//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_glasso(MatrixXd Y, MatrixXd Z, VectorXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();

	int nlam=opts_pen.nlam, p=opts.p, q=opts.q;
	opts.n=Y.rows();
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0);
	MatrixXi df = MatrixXi::Constant(p,nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
   
	betapath = MVR_glasso(Y, Z, df, lambda, likhd, opts.max_step, opts_pen.alpha, 
		   opts_pen.gamma_pen, opts_pen.pen, opts_pen.dfmax, opts.eps);
	return List::create(Named("betapath") = betapath, Named("df") = df, Named("likhd") = likhd);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_lasso(MatrixXd Y, MatrixXd Z1, MatrixXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.rows();

	int j,k,nlam=opts_pen.nlam, p=opts.p, q=opts.q, n=Y.rows();
	opts.n=n;
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0), Z = Z1, beta;
	MatrixXi df = MatrixXi::Constant(p*q,nlam, 0), activeA;
	VectorXd znorm, likhd0;
	MatrixXd likhd = MatrixXd::Constant(q, nlam, 0);
	likhd0.setZero(nlam);
	

    znorm = Z1.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) Z.col(j) = Z1.col(j)/znorm[j];
	for(k=0; k<q; k++){
        activeA.setZero(p,nlam);		
		beta = MVR_lasso(Y.col(k), Z, activeA, lambda.col(k), likhd0, opts.max_step, opts_pen.alpha, 
						opts_pen.gamma_pen, opts_pen.pen, opts_pen.dfmax, opts.eps);
		df.block(k*p, 0, p, nlam) = activeA;
		for(j=0;j<p;j++) beta.row(j) /= znorm[j];
		betapath.block(k*p, 0, p, nlam) = beta;
		likhd.row(k) = likhd0;
	}
	return List::create(Named("betapath") = betapath, Named("df") = df, Named("likhd") = likhd);
}