//[[Rcpp::depends(RcppEigen)]]
#include "mvr.h"

//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_colwise(MatrixXd Y, MatrixXd Z, MatrixXd W, VectorXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();

	int nlam=opts_pen.nlam, p=opts.p;
	opts.n=Y.rows();
	MatrixXi df = MatrixXi::Constant(p,nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	List fit;
   
	fit = MVR_colwise(Y, Z, W, df, lambda, likhd);
	return List::create(Named("betapath") = fit[0], Named("df") = df, Named("likhd") = likhd, Named("Cpath") = fit[1]);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_lasso(MatrixXd Y, MatrixXd Z1, MatrixXd W, MatrixXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.rows();

	int j,k,nlam=opts_pen.nlam, p=opts.p, q=opts.q, n=Y.rows(), pz=opts.pz;
	opts.n=n;
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0), Z = Z1, beta, Cnew, Cpath;
	MatrixXi df = MatrixXi::Constant(p*q,nlam, 0), activeA;
	VectorXd znorm, likhd0;
	MatrixXd likhd = MatrixXd::Constant(q, nlam, 0);
	likhd0.setZero(nlam);
	List fit;
	if(pz)	Cpath.setZero(pz*q,nlam);
	else Cpath.setZero(q,nlam);

    znorm = Z1.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) Z.col(j) = Z1.col(j)/znorm[j];
	for(k=0; k<q; k++){
        activeA.setZero(p,nlam);		
		fit = MVR_lasso(Y.col(k), Z, W, activeA, lambda.col(k), likhd0);
		beta = fit[0];
		Cnew = fit[1];
		df.block(k*p, 0, p, nlam) = activeA;
		for(j=0;j<p;j++) beta.row(j) /= znorm[j];
		betapath.middleRows(k*p, p) = beta;
		likhd.row(k) = likhd0;
		if(pz) Cpath.middleRows(k*pz,pz) = Cnew;
	}
	return List::create(Named("betapath") = betapath, Named("df") = df, Named("likhd") = likhd, Named("Cpath") = Cpath);
}
