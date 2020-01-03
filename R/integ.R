
##--------------without sparsity----------------------##
integ <- function(Y,X,g=1,r1=NULL,r2=NULL,r3=NULL,SABC=NULL,intercept=TRUE,mu=NULL,eps=1e-4,max_step=20){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  nx <- dim(X)[2]
  p <- nx/g
  if(is.null(r1)) r1 <- 2 
  if(is.null(r2)) r2 <- 2
  if(is.null(r3)) r3 <- 2

  # initial A,B,C,S
  if(is.null(SABC)){
    set.seed(1)
    A = rbind(diag(r1), matrix(0,p-r1,r1))
    B = rbind(diag(r2), matrix(0,g-r2,r2))
    C = rbind(diag(r3), matrix(0,q-r3,r3))
    S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
  }
  else{
    A = SABC$A
    B = SABC$B
    C = SABC$C
    S = SABC$S
  }
  if(!intercept | is.null(mu)) mu = rep(0,q)
  opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,g=g)
  fit = EstInteg(Y,X,S,A,B,C,mu,opts)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              mu = fit$mu,
              Y = Y,
              X = X
              )
         )
}