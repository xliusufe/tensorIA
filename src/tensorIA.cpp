//[[Rcpp::depends(RcppEigen)]]
#include "IA.h"
//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
// [[Rcpp::export]]
List EstInteg(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, VectorXd mu, List optsList)
{	
	opts.n = as<int>(optsList["n"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);	
	opts.p = as<int>(optsList["p"]);
	opts.K = as<int>(optsList["g"]);
	opts.q = as<int>(optsList["q"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.intercept = as<int>(optsList["intercept"]);
	
	double  likhd0 = pow(10, 6), likhd1 = 0;
	int n = opts.n,q = opts.q;
	MatrixXd Dnew,Cnew,Anew,Bnew,Snew, Y1=Y;
	VectorXi convergence1;
	VectorXd Ones;
	Ones.setOnes(n);
	List fit;
	int step = 0;
	while (step<opts.max_step) {
		convergence1 = VectorXi::Constant(4, 1);
		step = step + 1;
		Snew = updateS(Y1, Z, A, B, C);
		Dnew = C * Snew*kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;
		fit = updateC(Y1, Z, A, B, C, S);
		Cnew = fit[0];
		Snew = fit[1];
		Dnew = Cnew * Snew *kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			C = Cnew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;
		fit = updateA(Y1, Z, A, B, C, S);
		Anew = fit[0];
		Snew = fit[1];
		Dnew = C * Snew *kroneckerProduct(B.transpose(), Anew.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			A = Anew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[2]=0;
		fit = updateB(Y1, Z, A, B, C, S);
		Bnew = fit[0];
		Snew = fit[1];		
		Dnew = C * Snew *kroneckerProduct(Bnew.transpose(), A.transpose());
		likhd1 = (Y1  - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			B = Bnew;
			S = Snew;
			if(opts.intercept){
				mu = (Y - Z * Dnew.transpose()).colwise().sum()/n;
				Y1 = Y - kroneckerProduct(Ones,mu);
				Y1.resize(n,q);	
			}
			if ((likhd0 - likhd1) / likhd0<opts.eps) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[3]=0;
		if(convergence1.sum()==0) break;
	}
	Dnew = C * S * kroneckerProduct(B.transpose(), A.transpose());
	return List::create(Named("likhd") = likhd0, Named("Dnew") = Dnew, Named("A") = A, Named("B") = B, Named("C") = C,Named("S") = S, Named("mu") = mu);
}

//----------------------------------------------------------------**
//***--------------------EstimationD3 directly--------------------**
// [[Rcpp::export]]
List EstFR(MatrixXd Y, MatrixXd X)
{
  int j,q = Y.cols(),p=X.cols();
  MatrixXd Dnew = MatrixXd::Constant(q, p, 0), Z = produceX2(X);
  HouseholderQR<MatrixXd> qr;
  qr.compute(Z);
  MatrixXd R = qr.matrixQR().triangularView<Upper>();
  MatrixXd Q = qr.householderQ();
  MatrixXd RQ = UpTriangularInv(R) * Q.transpose();  
  MatrixXd temp, IDEN = MatrixXd::Identity(p, p);
  if (pow(condition_numberQRSym(R),2) > 1e10){
    temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
    for (j = 0; j < q; j++) Dnew.row(j) = (temp.colPivHouseholderQr().solve(Z.transpose()*Y.col(j))).transpose();
  }
  else
    for (j = 0; j < q; j++) Dnew.row(j) = (RQ * Y.col(j)).transpose();
  double likhd = pow((Y - Z * Dnew.transpose()).norm(),2);
  return List::create(Named("likhd") = likhd, Named("Dnew") = Dnew);
}