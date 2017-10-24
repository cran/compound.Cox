uni.score=function(t.vec, d.vec, X.mat,d0=0){
  n=length(t.vec)
  X.mat=as.matrix(X.mat)
  p=ncol(X.mat)
  atr_t=(matrix(t.vec,n,n,byrow=TRUE)>=matrix(t.vec,n,n))
  S0=atr_t%*%matrix(1,n,p)
  S1=atr_t%*%X.mat
  S2=atr_t%*%X.mat^2
  S=d.vec%*%(X.mat-S1/S0)
  V=d.vec%*%(S2/S0-(S1/S0)^2)
  Z=as.vector( S/(sqrt(V)+d0) )
  P=1-pchisq(Z^2,df=1)
  beta_est=as.vector(S/V) ## S/V is the one-step estimator
  names(Z)=names(P)=names(beta_est)=colnames(X.mat)
  list(beta_est=beta_est,Z=Z,P=P) 
}
