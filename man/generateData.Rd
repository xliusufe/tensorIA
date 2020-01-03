\name{generateData}
\alias{generateData}
\docType{package}
\title{
  Generate data from multivariate regression model.
}
\description{
  Generate data for a multivariate regression model.
}

\usage{
generateData(n, q, p, g, D3, SigmaX=NULL, SigmaE=NULL, 
             mu=NULL, sigma2=NULL, seed_id=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.}
  
  \item{q}{The number of responses, \eqn{q\geq1}.}
  
  \item{p}{The number of covariates, \eqn{p\geq1}.}
  
  \item{g}{The number of groups.}
  
  \item{D3}{The mode of unfolding \eqn{D_{(3)}}.}
  
  \item{sigma2}{err variance. Default is 0.2.}
  
  \item{SigmaX}{A \eqn{pg\times pg} positive-definition matrix, which is the covariance matrix of covariates \eqn{X}.}
  
  \item{SigmaX}{A \eqn{q\times q} positive-definition matrix, which is the covariance matrix of error \eqn{E}.}
  
  \item{seed_id}{A positive integer, the seed for generating the random numbers.}
}

\details{This function gives coefficients of multivariate regression. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Y}{Response, a \eqn{n\times q}-matrix.}
  
  \item{X}{Design matrix, a \eqn{n\times pg}-matrix.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
 Integrative analysis based on tensor modelling.
}
\keyword{ GWAS; Integrative analysis; Tensor estimation; Tucker decomposition. }
\seealso{
  mam_sparse
}
\examples{
  # Example 1
  
  D3 <- matrix(runif(30, 0.7, 1), 2, 15)
  mydata <- generateData(200, 3, 5, 5, D3)
    
  Y <- mydata$Y  
  X <- mydata$X
    
    
  # Example 2
  n <- 500
  p <- 10
  q <- 10
  g <- 10
  r10=r20=r30=2
  S3 <- matrix(runif(r10*r20*r30,3,7),nrow = r30)
  T1 <- matrix(rnorm(s0*r10),nrow = s0)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(g*r20),nrow = g)
  U2 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r30),nrow = q)
  U3 <- qr.Q(qr(T1))
  D3 <- U3\%*\%S3\%*\%t(kronecker(U2,U1))
  mydata <- generateData(n,q,p,g,D3)
}
