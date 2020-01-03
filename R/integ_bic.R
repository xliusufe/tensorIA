
##--------------main by BIC without sparsity----------------------##
integ_bic <- function(Y,X,method,r1_index,r2_index,r3_index,S,A,B,C,mu,opts){
  n <- opts$n
  q <- opts$q
  p <- opts$p
  RSS = NULL
  for(K in K_index){
    for(r3 in r3_index){
      for(r2 in r2_index){
        for(r1 in r1_index){
          fit = EstInteg(Y,X,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),mu,opts)
          df = r1*r2*r3+p*r1+K*r2+q*r3-r1^2-r2^2-r3^2
          loglikelih =  -n*q * (log(2*pi) + log(fit$likhd))
          if(method=="BIC"){
            bic = loglikelih + log(n*q)*df
          }
          if(method=="AIC") bic = loglikelih + 2*df
          if(method=="EBIC"){
            bic = loglikelih + log(n*q)*df
            bic = bic + 2*(lgamma(p+1) - lgamma(df+1) - lgamma(p-df+1))
          }
          if(method=="GCV") bic = loglikelih/(1-df/n)^2
          RSS = c(RSS,bic)
        }
      }
    }
  }
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  #---------------- The estimation after selection ---------------------#
  fit = EstInteg(Y,X,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),mu,opts)  
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