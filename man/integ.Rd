\name{integ}
\alias{integ-function}
\alias{integ}
\docType{package}
\title{
  Integrative analysis for GWAS data.
}
\description{
  Fit a multivariates grouped gression model by treating coefficients as a order 3 tensor, without aparsity assumptions, and given ranks  \eqn{r_1, r_2, r_3}.
}

\usage{integ(Y, X, g = 1, r1 = NULL, r2 = NULL, r3 = NULL, SABC = NULL, 
           intercept = TRUE, mu = NULL, eps = 1e-4, max_step = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{g}{The number of groups. Default is \code{1}.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}

  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}
  
  \item{mu}{A user-specified initial of intercept(s), a \eqn{q}-vector. Default is \code{0}.}
    
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 20.}
}

\details{This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{mu}{Estimator of intercept \eqn{\mu}.}
  
  \item{rss }{Residual sum of squares (RSS).}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
Integrative analysis based on tensor modelling.
}
\keyword{GWAS; Integrative analysis; Tensor estimation; Tucker decomposition. }
\seealso{
  integ_dr
}
\examples{ 
  n <- 200
  p <- 5
  q <- 5
  g <- 5
  D3 <- matrix(runif(p*q*g, 0.7, 1), q, p*g) # tensor with size 5*5*5
  X <- matrix(runif(n*p*g), nrow = n)
  eps <- matrix(rnorm(n*q),n,q)
  Y <- X\%*\%t(D3)  + eps
  
  fit <- integ(mydata$Y, mydata$X)
  D3hat <- fit$Dnew
  D2hat <- TransferModalUnfoldings(D3hat,3,2,p,g,q)
}