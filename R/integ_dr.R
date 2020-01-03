
##--------------main by BIC without sparsity----------------------##
integ_dr <- function(Y,X,g=1,method="BIC",ncv=10,r1_index=NULL,r2_index=NULL,r3_index=NULL,SABC=NULL,
                   intercept=TRUE,mu=NULL,eps=1e-4,max_step=20){

  n <- dim(Y)[1]
  q <- dim(Y)[2]
  nx <- dim(X)[2]
  p <- nx/g
  if(is.null(r1_index)) r1_index = 1:min(ceiling(log(n)),p)
  if(is.null(r2_index)) r2_index = 1:min(ceiling(log(n)),g)
  if(is.null(r3_index)) r3_index = 1:min(ceiling(log(n)),q)
  #---------------- The selection by BIC  ---------------------#  
  if(is.null(SABC)){
    set.seed(1)
    r1_max = max(r1_index) 
    r2_max = max(r2_index) 
    r3_max = max(r3_index) 
    A = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
    B = rbind(diag(r2_max), matrix(0,g-r2_max,r2_max))
    C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
    S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
  }
  else{
    A = SABC$A
    B = SABC$B
    C = SABC$C
    S = SABC$S
  }
  if(!intercept | is.null(mu)) mu = rep(0,q)
  opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,g=g)
  
  if((max(r1_index)>dim(A)[2])|(max(r2_index)>dim(B)[2])|(max(r3_index)>dim(C)[2]))
    stop("maximum number of index sequence of r1, r2, and r3 must not be larger than A, B, and C, respectively !")

  if(method=="CV") fit_dr = integ_cv(Y,X,ncv,r1_index,r2_index,r3_index,S,A,B,C,mu,opts)
  else fit_dr = integ_bic(Y,X,method,r1_index,r2_index,r3_index,S,A,B,C,mu,opts)
  
  return(fit_dr)
}