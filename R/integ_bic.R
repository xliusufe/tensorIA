
##--------------main by BIC without sparsity----------------------##
integ_bic <- function(Y,X,method,r1_index,r2_index,r3_index,S,A,B,C,mu,opts){
  n <- opts$n
  q <- opts$q
  p <- opts$p
  g <- opts$g
  RSS = NULL
  for(r3 in r3_index){
    opts$r3=r3
    for(r2 in r2_index){
      opts$r2=r2
      for(r1 in r1_index){
        opts$r1=r1
        fit = EstInteg(Y,X,as.matrix(S[1:r3,1:(r1*r2)]),as.matrix(A[,1:r1]),as.matrix(B[,1:r2]),as.matrix(C[,1:r3]),mu,opts)
        df = r1*r2*r3+p*r1+g*r2+q*r3-r1^2-r2^2-r3^2
        loglikelih =  n*q * log(fit$likhd/(n*q))
        bic <- switch (method,
                       BIC = loglikelih + log(n*q)*df,
                       AIC = loglikelih + 2*df,
                       GCV = fit$likhd*n*q/(n*q-df)^2,
                       EBIC = loglikelih + log(n*q)*df + 2*(lgamma(p*g*q+1) 
                                         - lgamma(df+1) - lgamma(p*g*q-df+1))
                       )        
        RSS = c(RSS,bic)
      }
    }
  }
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  #---------------- The estimation after selection ---------------------#
  opts$r1=r1_opt
  opts$r2=r2_opt
  opts$r3=r3_opt
  fit = EstInteg(Y,X,as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),as.matrix(A[,1:r1_opt]),as.matrix(B[,1:r2_opt]),as.matrix(C[,1:r3_opt]),mu,opts)  
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              mu = fit$mu,
              rk_opt=c(r1_opt,r2_opt,r3_opt),
              selected=selected,
              Y = Y,
              X = X
              )
         )
}