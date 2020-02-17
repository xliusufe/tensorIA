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
	int r1;
	int r2;
	int r3;
	int p;
	int q;
	int K;
	int n;
	double eps;
	double eps1;	
	int max_step;
    int max_step1;	
	int intercept;
}opts;


//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
VectorXi SEQ(int min, int max, int by)
{
	int leng = static_cast<int>(floor((max - min) / by + pow(3, -12)) + 1);
	VectorXi C = VectorXi::Constant(leng, 0);
	for (int i = 0; i < leng; i++)
		C[i] = min + by * i;
	return C;
}
//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
double dotproduct(VectorXd a, VectorXd b, int n){
	double out=0.0;
	for(int i=0;i<n;i++) out += a[i]*b[i];
	return out;
}

double innerProduct(VectorXd A, VectorXd B, int n){
	double out=0.0;
	for(int i=0;i<n;i++) out += A[i]*B[i];
	return out;
}

double innerProductMatrix(MatrixXd A, MatrixXd B){
	return (A.array()*B.array()).sum();
}

//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
VectorXd SEQN(int min, int n, int by)
{
	// by is the increment with by >0	
	VectorXd C = VectorXd::Constant(n,0);
	for (int i = 0; i < n; i++)	C[i] = min + by * i;
	return C;
}

//----------------------------------------------------------------**
//***----------------------uppertriangleNorm1 of AA---------------**
double uppertriangleNorm1AA(MatrixXd A)
{
	int i, j, k, p, ii;
	double norm1 = 0, temp1, temp2;
	p = A.cols();

	for (k = 0; k < p; k++) {
		temp2 = 0;
		for (j = 0; j < p; j++) {
			temp1 = 0;
			ii = k < j ? k : j;
			for (i = 0; i <= ii; i++) temp1 += A(i, j) * A(i, k);
			temp2 += fabs(temp1);
		}
		if (temp2 > norm1) norm1 = temp2;
	}
	return norm1;
}
//----------------------------------------------------------------**
//***----------------------uppertriangleNorm1 of A----------------**
double uppertriangleNorm1(MatrixXd A)
{
	int j, k, p;
	double norm1 = 0, temp1;
	p = A.cols();

	for (k = 0; k < p; k++) {
		temp1 = 0;
		for (j = 0; j <= k; j++) {
			temp1 += fabs(A(j, k));
		}
		if (temp1 > norm1) norm1 = temp1;
	}
	return norm1;
}

//----------------------------------------------------------------**
//***----------------------UpTriangularInv------------------------**
MatrixXd UpTriangularInv(MatrixXd A)
{
	int i, j, k,n=A.cols();
	MatrixXd B = MatrixXd::Constant(n, n, 0);
	for (i = 0; i < n; i++) B(i, i) = 1;
	for (i = n - 1; i >= 0; i--)//rows
	{
		if (A(i, i) != 1)
			for (j = i; j<n; j++)
				B(i, j) = B(i, j) / A(i, i);
		if (i>0)
		{
			for (j = i; j<n; j++)// columns
				for (k = 0; k<i; k++)// rows
					B(k, j) = B(k, j) - A(k, i) * B(i, j);
		}
	}
	return B;
}
//----------------------------------------------------------------**
//***---------------------- Norm1 of a Matrix---------------------**
double MatrixNorm1(MatrixXd A)
{
  int j, k, p;
  double norm1 = 0, temp1;
  p = A.cols();
  
  for (k = 0; k < p; k++) {
    temp1 = 0;
    for (j = 0; j < p; j++) temp1 += fabs(A(j,k));
    if (temp1 > norm1) norm1 = temp1;
  }
  return norm1;
}
//----------------------------------------------------------------**
//***---------------------- Q*R of qr decomposition --------------**
MatrixXd QbyR(MatrixXd Q, MatrixXd R, int isQR)
{
	//isQR=1 denotes Q*R; otherwise R*Q
	int i, j, k, p = R.cols(),n;
	double temp1;
	MatrixXd A = Q;
	if(isQR){
		n = Q.rows();
		for (k = 0; k < p; k++)
			for (j = 0; j < n; j++) {
				temp1 = 0;
				for (i = 0; i <= k; i++) temp1 += Q(j, i)*R(i, k);
				A(j,k) = temp1;
			}
	}
	else{
		n = Q.cols();
		for (k = 0; k < n; k++)
			for (j = 0; j < p; j++) {
				temp1 = 0;
				for (i = j; i < p; i++) temp1 += R(j, i) * Q(i, k);
				A(j,k) = temp1;
			}
	}
	return A;
}
//----------------------------------------------------------------**
//***---------------------- R*R of Upper triangle ----------------**
MatrixXd tRbyR(MatrixXd R)
{
  int i, j, k, ii, p;
  double temp1;
  MatrixXd A = R;
  p = R.cols();
  
  for (k = 0; k < p; k++) {
    for (j = 0; j < p; j++) {
      temp1 = 0;
	  ii = k < j ? k : j;
      for (i = 0; i <= ii; i++) temp1 += R(i, j) * R(i, k);
      A(j,k) = temp1;
    }
  }
  return A;
}
//----------------------------------------------------------------**
//***----------------------cbind----------------------------------**
MatrixXd cbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n = A.rows();
	int p1 = A.cols();
	int p2 = B.cols();
	MatrixXd C = MatrixXd::Constant(n, p1 + p2, 0);
	C.block(0, 0, n, p1) = A;
	C.block(0, p1, n, p2) = B;
	return C;
}
//----------------------------------------------------------------**
//***----------------------rbind----------------------------------**
MatrixXd rbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n1 = A.rows();
	int n2 = B.rows();
	int p = A.cols();
	MatrixXd C = MatrixXd::Constant(n1 + n2, p, 0);
	C.block(0, 0, n1, p) = A;
	C.block(n1, 0, n2, p) = B;
	return C;
}

//----------------------------------------------------------------**
//***----------------------extract columns------------------------**
MatrixXd extractRows0(MatrixXd A, VectorXi b)
{
	int n = A.rows(), p = A.cols(), p2=b.sum(), j, count = 0;
	if(p2>p) stop("The length of index b must be not great than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(p2, p, 0);
	for(j=0; j< n; j++)  if(b[j])   C.row(count++) = A.row(j);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract columns------------------------**
MatrixXd extractRows(MatrixXd Z, int K, VectorXi b)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>n) stop("The length of index b must be not greater than the number of columns of Z!");
	MatrixXd C = MatrixXd::Constant(K, p, 0);
	for(i=0;i<K;i++) C.row(i) = Z.row(b[i]);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractColsZ(MatrixXd Z, int p, int K, VectorXi b)
{
	int n = Z.rows(), p2=b.sum(), i, j, count;
	if(p2>p) stop("The length of index b must be not great than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(n, p2*K, 0);
	for(i=0;i<K;i++){
	  count = 0;
	  for(j=0; j< p; j++) if(b[j]) C.col(i*p2 + (count++)) = Z.col(i*p+j);
	}
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractCols(MatrixXd Z, int K, VectorXi b)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>p) stop("The length of index b must be not greater than the number of columns of Z!");
	MatrixXd C = MatrixXd::Constant(n, K, 0);
	for(i=0;i<K;i++) C.col(i) = Z.col(b[i]);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractCols_by(MatrixXd Z, int min, int K, int by)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>p) stop("The length of index b must be not greater than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(n, K, 0);
	for(i=0;i<K;i++) C.col(i) = Z.col(min+by*i);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract subvector----------------------**
VectorXd subvector_extract(VectorXd U, int min, int n, int by){
	VectorXd V = VectorXd::Constant(n,0);
	for(int i=0; i<n; i++) V[i] = U[min+by*i];
	return V;
}

//----------------------------------------------------------------**
//***----------------------update subvector-----------------------**
void subvector_updated(VectorXd &U, VectorXd V, int min, int n, int by){
	for(int i=0; i<n; i++) U[min+by*i] = V[i];
}

//----------------------------------------------------------------**
//***----------------------update subvector-----------------------**
MatrixXd myQR(MatrixXd U, int k){
	int i, n=U.rows(), flag=0;
	MatrixXd R, Q, IDEN = MatrixXd::Identity(k, k);
	HouseholderQR<MatrixXd> qr;
	qr.compute(U);
	Q = (qr.householderQ());
	Q = Q.block(0,0,n,k);
	R = qr.matrixQR().triangularView<Upper>();
	for(i=0;i<k;i++) if(R(i,i)<0) {IDEN(i,i) = -1;	flag=1;}
	if(flag) Q = Q*IDEN;
	return Q;
}
//----------------------------------------------------------------**
//***--------QRcondition_number for design matrix-----------------**
double condition_numberQR(MatrixXd R)
{
	return uppertriangleNorm1AA(R) * uppertriangleNorm1AA(UpTriangularInv(R));
}
//----------------------------------------------------------------**
//***------QRcondition_number for symetric matrix-----------------**
double condition_numberQRSym(MatrixXd R)
{
  return uppertriangleNorm1(R) * uppertriangleNorm1(UpTriangularInv(R));
}
//----------------------------------------------------------------**
//***----------------solve linear system by QR--------------------**
VectorXd solveEquationQR(MatrixXd A, VectorXd b)
{
  //solve linear system Ax = b, A = A.transpose();
  int p = b.size();
  HouseholderQR<MatrixXd> qr;
  qr.compute(A);
  MatrixXd R, Q;
  Q = qr.householderQ();
  R = qr.matrixQR().triangularView<Upper>();
  VectorXd X;
  
  if (condition_numberQRSym(R) > 1e10){
    MatrixXd temp, IDEN = MatrixXd::Identity(p, p);
    temp = A + (IDEN.array()*1e-4).matrix();
    X = temp.colPivHouseholderQr().solve(b);
  }
  else X = QbyR(Q.transpose(),UpTriangularInv(R),0)*b;
  return X;
}
//----------------------------------------------------------------**
//***--------------------transfer modal of unfoldings-------------**
// [[Rcpp::export]]
MatrixXd TransferModalUnfoldings(MatrixXd S, int d1, int d2, int r1, int r2, int r3)
{
  //From S_(d1) to S_(d2)
  int j;
  MatrixXd S1,S_3;
  if (d1 == 3) {
    if (d2 == 1){
      S1 = S.row(0).transpose();
      S1.resize(r1, r2); 
      for (j = 1; j < r3;j++) {
        S_3 = S.row(j).transpose();
        S_3.resize(r1, r2);
        S1 = cbind_rcpp(S1, S_3);// S3 is r3 *(r2r1) matrix
      }
    }	
    if (d2 == 2) {
      S1 = S.block(0, 0, r3, r1).transpose();
      S1.resize(1,r1*r3);
      for (j = 1; j < r2; j++) {
        S_3 = S.block(0, j*r1, r3, r1).transpose();
        S_3.resize(1,r1*r3);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  
  if (d1 == 2) {
    if (d2 == 1) {
      S1 = S.block(0, 0, r2, r1).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1,S.block(0,j*r1,r2,r1).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r2, r1).transpose();
      S1.resize(1, r2*r1);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r1, r2, r1).transpose();
        S_3.resize(1, r2*r1);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
    
  }
  
  if (d1 == 1) {
    if (d2 == 2) {
      S1 = S.block(0, 0, r1, r2).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1, S.block(0, j*r2, r1, r2).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r1, r2);
      S1.resize(1, r1*r2);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r2, r1, r2);
        S_3.resize(1, r1*r2);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  return S1;
}

//----------------------------------------------------------------**
//***----------------------reassign columns of a matrix-----------**
MatrixXd submatrix_col(MatrixXd A, VectorXi b)
{
	int n = A.rows();
	int p = b.size();
	int j;
	MatrixXd C = MatrixXd::Constant(n, p, 0);
	for (j = 0; j < p; j++)
		C.col(j) = A.col(b[j]-1);
	return C;
}
//----------------------------------------------------------------**
//***----------------------reassign rows of a matrix--------------**
MatrixXd submatrix_row(MatrixXd A, VectorXi b)
{
	int p = A.cols();
	int n = b.size();
	int i;
	MatrixXd C = MatrixXd::Constant(n, p, 0);
	for (i = 0; i < n; i++)
		C.row(i) = A.row(b[i] - 1);
	return C;
}

//----------------------------------------------------------------**
//***--------------------updateS----------------------------------**
MatrixXd updateS(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C)
{
	int r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), d = r1*r2*r3;
	int k,k1,j,j1;
	MatrixXd ztilde = Z * kroneckerProduct(B, A), vectorS;
	VectorXd U;
	U.setZero(d);
	MatrixXd  V = MatrixXd::Constant(d, d, 0);
	for (k = 0; k < r1*r2; k++) {
		for (j = 0; j < r3; j++) {
			U[k*r3 + j] = ztilde.col(k).transpose()*Y*C.col(j);
			for (k1 = 0; k1 < r1*r2; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(k*r3 + j, k1*r3 + j1) = kroneckerProduct(
						ztilde.col(k1).array()*ztilde.col(k).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
			}
		}
	}
	vectorS = V.colPivHouseholderQr().solve(U);
	vectorS.resize(r3, r1*r2);
	return vectorS;
}
//----------------------------------------------------------------**
//***--------------------updateC----------------------------------**
List updateC(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r3 = C.cols(),q = C.rows(), j,kp;
	MatrixXd ztilde = Z * kroneckerProduct(B, A);
	MatrixXd StZ = ztilde * S.transpose();
	MatrixXd Cnew = MatrixXd::Constant(q, r3, 0);
	HouseholderQR<MatrixXd> qr;
	qr.compute(StZ);
	MatrixXd R = qr.matrixQR().triangularView<Upper>();
	MatrixXd Q = qr.householderQ();	
	kp = StZ.cols();
	MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
	if (pow(condition_numberQRSym(R),2) > 1e10){
	  temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
	  for (j = 0; j < q; j++) Cnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
	}
	else
	  for (j = 0; j < q; j++) Cnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();
	
	qr.compute(Cnew);
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	return List::create(Named("Cnew") = VQ.block(0, 0, q, r3), Named("Snew") = QbyR(S, VR.block(0, 0, r3, r3), 0));
}
//----------------------------------------------------------------**
//***--------------------updateA----------------------------------**
List updateA(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), p = A.rows(), K = B.rows();
	int d = r1 * p, t1,t2,t3,t4,j;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zbw1, zbw2, vectorA;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);
	for (t2 = 0; t2<p; t2++) {
		Zt2 = Z.col(t2);
		for (j = 1; j<K; j++)	Zt2 = cbind_rcpp(Zt2, Z.col(j*p + t2));
		for (t1 = 0; t1<r1; t1++) {
			Wt1 = W.col(t1);
			for (j = 1; j<r2; j++) Wt1 = cbind_rcpp(Wt1, W.col(j*r1 + t1));
			zbw1 = Zt2 * B*(Wt1.transpose());
			tU[t2*r1 + t1] = (Y.array()*zbw1.array()).sum();
			for (t4 = 0; t4<p; t4++) {
				Zt4 = Z.col(t4);
				for (j = 1; j<K; j++) Zt4 = cbind_rcpp(Zt4, Z.col(j*p + t4));
				for (t3 = 0; t3<r1; t3++) {
					Wt3 = W.col(t3);
					for (j = 1; j<r2; j++)	Wt3 = cbind_rcpp(Wt3, W.col(j*r1 + t3));
					
					zbw2 = Zt4 * B*(Wt3.transpose());
					tV(t2*r1 + t1, t4*r1 + t3) = (zbw1.array()*zbw2.array()).sum();
				}
			}
		}
	}
	vectorA = solveEquationQR(tV,tU);
	vectorA.resize(r1, p);
	HouseholderQR<MatrixXd> qr;
	qr.compute(vectorA.transpose());
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	MatrixXd S1new, S1 = TransferModalUnfoldings(S,3,1,r1,r2,r3);
	S1new = QbyR(S1, VR.block(0, 0, r1, r1), 0);
	return List::create(Named("Anew") = VQ.block(0, 0, p, r1), Named("Snew") = TransferModalUnfoldings(S1new,1,3,r1,r2,r3));
}
//----------------------------------------------------------------**
//***--------------------updateB----------------------------------**
List updateB(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), p = A.rows(), K = B.rows(), q = C.rows(), n = Z.rows();	
	int d = r2*K, t1,t2,t3,t4;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);
	for (t2 = 0; t2<K; t2++) {
		Zt2 = Z.block(0, t2*p, n, p);  
		for (t1 = 0; t1<r2; t1++) {
			Wt1 = W.block(0, t1*r1, q, r1);  			
			zaw1 = Zt2 * A*(Wt1.transpose()); 
			tU[t2*r2 + t1] = (Y.array()*zaw1.array()).sum();
			for (t4 = 0; t4<K; t4++) {
				Zt4 = Z.block(0, t4*p, n, p); 
				for (t3 = 0; t3<r2; t3++) {
					Wt3 = W.block(0, t3*r1, q, r1); 					
					zaw2 = Zt4 * A*(Wt3.transpose()); 
					tV(t2*r2 + t1, t4*r2 + t3) = (zaw1.array()*zaw2.array()).sum();
				}
			}
		}
	}
	vectorB = solveEquationQR(tV, tU);
	vectorB.resize(r2, K);	
	HouseholderQR<MatrixXd> qr;
	qr.compute(vectorB.transpose());
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	MatrixXd S2new, S2 = TransferModalUnfoldings(S,3,2,r1,r2,r3);
	S2new = QbyR(S2, VR.block(0, 0, r2, r2), 0);
	return List::create(Named("Bnew") = VQ.block(0, 0, K, r2), Named("Snew") = TransferModalUnfoldings(S2new,2,3,r1,r2,r3));
}
	
