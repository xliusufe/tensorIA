\name{integ_dr}
\alias{integ_dr-function}
\alias{integ_dr}
\docType{package}
\title{
  Integrative analysis for GWAS data without sparsity assumption, and with ranks selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}
\description{
  Fit a multivariates grouped gression model by treating coefficients as a order 3 tensor, without aparsity assumptions, and with ranks \eqn{r_1, r_2, r_3} selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
} 

\usage{integ_dr(Y, X, g = 1, method = "BIC", ncv = 10, r1_index = NULL,
       r2_index = NULL, r3_index = NULL, SABC = NULL, intercept = FALSE, mu = NULL,
       eps = 1e-4, max_step = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}

  \item{g}{The number of groups. Default is \code{1}.}

  \item{method}{The method to be applied to select parameters.  Either \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is  
             \code{"BIC"}, \code{ncv} is useless.}
    
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                   \eqn{r_1} is the first dimension of single value matrix of the 
                   tensor. Default is 
                   \code{r1_index}\eqn{=1,\cdots,\min(\log(n)\rceil,p)}.
                   }
   
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where 
                  \eqn{r_2} is the second dimension of single value matrix of the 
                  tensor. Default is \code{r2_index}\eqn{=1,\cdots,\min(\log(n)\rceil,g)}.
                  }
   
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\log(n)\rceil,q)}.
                  }
                   
  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}, which is a list with values \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}
        
  \item{intercept}{Should intercept(s) be fitted (\code{TRUE}) or set to zero (default=\code{FALSE})?}
  
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
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{selected}{Which \eqn{\lambda} is selection.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
Integrative analysis based on tensor modelling.
}
\seealso{
  integ
}
\examples{ 
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
  D3 <- U3\%*\%S3\%*\%t(kronecker(U2,U1))
  X <- matrix(rnorm(n*p*g), nrow = n)
  eps <- matrix(rnorm(n*q),n,q)
  Y <- X\%*\%t(D3)  + eps
  
  fit <- integ_dr(Y, X, g)
  D3hat <- fit$Dnew
  D2hat <- TransferModalUnfoldings(D3hat,3,2,p,g,q)
  opt <- fit$rk_opt
}