
##--------------Estimation without Penalty----------------------##
integ_cv <- function(Y,X,ncv,r1_index,r2_index,r3_index,S,A,B,C,mu,opts){
  n <- opts$n
  q <- opts$q
  p <- opts$p
  len_cv = ceiling(n/ncv)  
  RSS = rep(0,length(r1_index)*length(r2_index)*length(r3_index))
  for(jj in 1:ncv){ # start CV
    cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
    if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
    Ytrain = Y[-cv.id,]
    Xtrain = X[-cv.id,]
    Ytest = Y[cv.id,]
    Xtest = X[cv.id,]
    nt = nrow(Ytest)
    RSS0 = NULL
    for(r3 in r3_index){
      opts$r3=r3
      for(r2 in r2_index){
        opts$r2=r2
        for(r1 in r1_index){
          opts$r1=r1
          fit = EstInteg(Ytrain,Xtrain,as.matrix(S[1:r3,1:(r1*r2)]),as.matrix(A[,1:r1]),as.matrix(B[,1:r2]),as.matrix(C[,1:r3]),mu,opts)
          Dnew  = fit$Dnew
          if(opts$intercept)
            RSS0 = c(RSS0,sum((Ytest - matrix(rep(fit$mu,each=nt),nt) - Xtest%*%t(Dnew))^2))
          else
            RSS0 = c(RSS0,sum((Ytest - Xtest%*%t(Dnew))^2))
        }
      }
    }
    RSS = RSS + RSS0
  } # end of CV
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
