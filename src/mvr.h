#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <list>
#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b) 
#define IDEX(a, b) (a>b?1:0) 
using namespace Rcpp;
using namespace Eigen;

//----------------------------------------------------------------**
//***----------------------parameters for penalization---------------**
struct Options
{
	int p;
	int q;
	int n;
	double eps;
	int max_step;
}opts;

struct Options_pen
{
	int pen; 
	int nlam;
	int dfmax;
	int isPenColumn;
	double lam_max;
	double lam_min;
	double alpha;
	double gamma_pen;
}opts_pen;

//----------------------------------------------------------------**
//***--------------------penalty----------------------------------**
double penalties(double z, double v, double lambda, double alpha, double gamma, int penalty) {
	double beta=0,l1,l2;
	l1 = lambda*alpha; 
	l2 = lambda*(1-alpha);
	if (penalty==1){			  
		if (z > l1) beta = (z-l1)/(v*(1+l2));
		if (z < -l1) beta = (z+l1)/(v*(1+l2));
	}
	if (penalty==2){
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
		if (fabs(z) <= l1) beta = 0;
		else if (fabs(z) <= gamma*l1*(1+l2)) beta = s*(fabs(z)-l1)/(v*(1+l2-1/gamma));
		else beta = z/(v*(1+l2));
	}
	if (penalty==3){
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
		if (fabs(z) <= l1) beta = 0;
		else if (fabs(z) <= (l1*(1+l2)+l1)) beta = s*(fabs(z)-l1)/(v*(1+l2));
		else if (fabs(z) <= gamma*l1*(1+l2)) beta = s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2));
		else beta = z/(v*(1+l2));
	}
	return(beta);
}
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
VectorXd setuplambdaMVR_colwise(MatrixXd Y, MatrixXd Z, int nlam, VectorXd setlam)
{
	int n=Y.rows(), p = Z.cols(), q = Y.cols(), j;
	double lam_max, lam_min, alpha, max_lam, max_tmp=0;
	VectorXd lambda, lambda1, znorm;	
	znorm = Z.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) max_tmp = MAX(max_tmp,(Y.transpose()*Z.col(j)).norm()/znorm[j]);
	//max_tmp = ((Y.transpose()*Z).colwise().norm().array()/Z.colwise().norm().array()).max();
	
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];
	max_tmp/=n*sqrt(q);
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}
	return lambda;
}
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
MatrixXd setuplambdaMVR_lasso(MatrixXd Y, MatrixXd Z1, int nlam, VectorXd setlam)
{
	int n=Y.rows(), p = Z1.cols(), q = Y.cols(), j, k;
	double lam_max, lam_min, alpha, max_lam, max_tmp;
	VectorXd lambda, lambda1, znorm;	
	MatrixXd lambda_all = MatrixXd::Constant(nlam, q, 0);
	
	znorm = Z1.colwise().norm()/sqrt(n);	
	for(k=0; k<q; k++){	
    	max_tmp=0;
		for(j=0; j<p; j++) max_tmp = MAX(max_tmp,fabs(Y.col(k).transpose()*Z1.col(j))/znorm[j]);	
		//max_tmp = ((Z1.transpose()*Y.col(k)).cwiseAbs()/znorm).max();
		lam_max = setlam[0];
		lam_min = setlam[1];
		alpha = setlam[2];
		max_tmp/=n;
		max_lam = lam_max * max_tmp / alpha;
		if (lam_min == 0) {
			lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
			lambda.setLinSpaced(nlam, 0, 0);
			lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
		}
		else {
			lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
			lambda = lambda1.array().exp();
		}
		lambda_all.col(k) = lambda;
	}
	return lambda_all;
}
//----------------------------------------------------------------**
//***----update the jth row of matrix A with penalty--------------**
VectorXd updateAj(VectorXd z, double lambda, double alpha, double gamma, int penalty)
{
	double znorm = z.norm();
	znorm = penalties(znorm, 1, lambda, alpha, gamma, penalty) / znorm;
	return znorm * z;
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
MatrixXd MVR_colwise(MatrixXd Y, MatrixXd Z1, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd)
{
/*
	Input:
	Y is n*q matrix
	Z is n*p matrix, which is standardized
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is q*p coefficient metrix
*/  
	int l,j, active, step, nlam = lambda.size();
	int n = Y.rows(), q = Y.cols(), p = Z1.cols();
	double lambda1,diffmax,diffnorm;
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, penalty = opts_pen.pen;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma_pen;

    MatrixXd Anew=MatrixXd::Constant(q, p, 0), Bnew, Z = Z1;
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0);
	VectorXd ajnew,gj, ajnorm, znorm, diff;
	
	MatrixXd r = Y;
	ajnorm.setZero(p);
	znorm = Z1.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) Z.col(j) = Z1.col(j)/znorm[j];	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<max_step) {
			step++;
			active = 0;
			diffmax=0;
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active = active + 1;
			if (active>dfmax) {
				Anew = MatrixXd::Constant(q, p, -9);
				return Anew.transpose();
			}
			for (j = 0; j < p; j++) {
				gj = r.transpose()*Z.col(j)/n + Anew.col(j);			
				ajnew = updateAj(gj, lambda1, alpha, gamma, penalty);
				diff = ajnew-Anew.col(j);
				diffnorm = diff.norm();
				if(diffnorm!=0){
					r -= kroneckerProduct(Z.col(j),diff.transpose());
					Anew.col(j) = ajnew;
					ajnorm[j] = ajnew.norm();
					if(diffnorm>diffmax) diffmax = diffnorm;
				}
			}
			if(diffmax<eps) break;
		}//end while
		for (j = 0; j<p; j++) 	if (ajnorm[j]) activeA(j,l) = 1;
		likhd[l] = (Y - Z * Anew.transpose()).squaredNorm();
		Bnew = Anew.transpose();
		for(j=0; j<p; j++) Bnew.row(j) /= znorm[j];
		Bnew.resize(p*q,1);
		betapath.col(l) = Bnew;
	}//end for	
	return betapath;
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
MatrixXd MVR_lasso(VectorXd Y, MatrixXd Z, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd)
{
/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is p-dimensional coefficient
*/  
	int l,j, active, step, n = Z.rows(), p = Z.cols(), nlam = lambda.size();
    double ajnew,gj, lambda1,diffmax, diff;
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, penalty = opts_pen.pen;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma_pen;
	
	MatrixXd beta = MatrixXd::Constant(p,nlam,0);;
	VectorXd Anew, ajnorm, r = Y, r1=Y;
	Anew.setZero(p);
	ajnorm = Anew;
	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step=0;
		while (step<max_step){
			step++;
			active = 0;
			diffmax = 0;
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active = active + 1;
			if (active>dfmax)	return VectorXd::Constant(p, -9);			
			for (j = 0; j < p; j++) {
				gj = r.dot(Z.col(j))/n + Anew[j];
				ajnew = penalties(gj, 1, lambda1, alpha, gamma, penalty);
				diff = ajnew-Anew[j];
				if(diff!=0){ 
             		r -= Z.col(j)*diff;	
					Anew[j] = ajnew;
					ajnorm[j] = fabs(ajnew);
					diff = fabs(diff); 
					if(diff>diffmax) diffmax = diff;
				}				
			}						
			if(diffmax<eps) break;
		}//end while
		likhd[l] = (Y - Z * Anew).squaredNorm();
		for (j = 0; j<p; j++) if (ajnorm[j]) activeA(j,l) = 1;
		beta.col(l) = Anew;
	}//end for	
	return beta;
}	
