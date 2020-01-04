# tensorIA
 Integrative analysis based on tensor modelling.
 
  For a grouped multivariates gression model, with or without aparsity assumptions, 
  treating the coefficients as a third-order tensor and borrowing Tucker decomposition to reduce the number of parameters.
  
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("xliusufe/tensorIA")

# Usage

   - [x] [tensorIA-manual](https://github.com/xliusufe/tensorIA/blob/master/inst/tensorIA-manual.pdf) ------------ Details of the usage of the package.
# Example

    library(tensorIA)

    n <- 200
	p <- 5
	q <- 5
	g <- 5
	r10 <- 2
	r20 <- 2
	r30 <- 2
	S3 <- matrix(runif(r10*r20*r30,3,7),nrow = r30)
	T1 <- matrix(rnorm(p*r10),nrow = p)
	U1 <- qr.Q(qr(T1))
	T1 <- matrix(rnorm(g*r20),nrow = g)
	U2 <- qr.Q(qr(T1))  
	T1 <- matrix(rnorm(q*r30),nrow = q)
	U3 <- qr.Q(qr(T1))
	D3 <- U3%*%S3%*%t(kronecker(U2,U1))
	X <- matrix(rnorm(n*p*g), nrow = n)
	eps <- matrix(rnorm(n*q),n,q)
	Y <- X%*%t(D3)  + eps
  
    fit <- integ_dr(Y, X, g, intercept=FALSE)
	D3hat <- fit$Dnew
	D2hat <- TransferModalUnfoldings(D3hat,3,2,p,g,q)
	opt <- fit$rk_opt	
 
 # References
Integrative analysis based on tensor modelling. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
