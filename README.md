# tensorIA
 Integrative analysis based on tensor modelling.
 
  For a multivariates grouped gression model, with or without aparsity assumptions, 
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
	D3 <- matrix(runif(p*q*g, 0.7, 1), q, p*g) # tensor with size 5*5*5
	X <- matrix(runif(n*p*g), nrow = n)
	eps <- matrix(rnorm(n*q),n,q)
	Y <- X%*%t(D3)  + eps
  
    fit <- integ_dr(mydata$Y, mydata$X)
	D3hat <- fit$Dnew
	D2hat <- TransferModalUnfoldings(D3hat,3,2,p,g,q)
	opt <- fit$rk_opt	
 
 # References
Integrative analysis based on tensor modelling. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
