uni.selection=function(t.vec, d.vec, X.mat, P.value=0.001,K=5){
  
  ######### univariate Cox ##########
  P=Symbol=beta_est=X.cut=Z=NULL
  p=ncol(X.mat)
  for(j in 1:p){
    res=summary(coxph(Surv(t.vec,d.vec)~X.mat[,j]))$coef
    if( res[5]<P.value ){
      Z=c(Z,res[4]) # Z-value for Wald test
      P=c(P,res[5]) # P-value for Wald test
      beta_est=c(beta_est,res[1])
      X.cut=cbind(X.cut,X.mat[,j])
      Symbol=c(Symbol,colnames(X.mat)[j])
    }
  }
  
  if(is.null(X.cut)){warning("no gene selected; increase P.value")}else{
    
    q=ncol(X.cut)
    n=length(t.vec)
    CC=X.cut%*%beta_est
    c.index=unname( survConcordance(  Surv(t.vec,d.vec)~CC  )$concordance )
    
    ##### Log-partial likelihood ######
    t.ot=t.vec[order(t.vec)]
    atr_t=(matrix(t.vec,n,n,byrow=TRUE)>=matrix(t.vec,n,n))
    l.func=function(g){
      l1=sum( (CC*g)[d.vec==1] )
      S0=sum( (log(atr_t%*%exp(CC*g)))[d.vec==1] )
      l1=l1-S0
      as.numeric( l1 )
    }
   
    ### Cross-validation ###
    CC.CV=NULL
    CVL=CVL_CV=0
    folds=cut(seq(1,n),breaks=K,labels=FALSE)
    for(k in 1:K){
      temp=which(folds==k)
      t_k=t.vec[-temp]
      d_k=d.vec[-temp]
      CC_k=CC[-temp]
      g_k=coxph(Surv(t_k,d_k)~CC_k)$coef[1]
      n_k=length(t_k)
      atr_t_k=(matrix(t_k,n_k,n_k,byrow=TRUE)>=matrix(t_k,n_k,n_k))
      l_k=sum( (CC_k*g_k)[d_k==1] )
      S0=sum( (log(atr_t_k%*%exp(CC_k*g_k)))[d_k==1] )
      l_k=l_k-S0
      CVL=CVL+l.func(g_k)-as.numeric(l_k)
      
      beta_CV=numeric(q)
      for(j in 1:q){
        res=summary(coxph(Surv(t_k,d_k)~X.cut[-temp,j]))$coef[1]
        beta_CV[j]=res
      }
      CC.CV=c(CC.CV,X.cut[temp,]%*%beta_CV)
      CC.CV_k=X.cut[-temp,]%*%beta_CV
      g_CV_k=coxph(Surv(t_k,d_k)~CC.CV_k)$coef[1]
      l_CV_k=sum( (CC.CV_k*g_CV_k)[d_k==1] )
      S0=sum( (log(atr_t_k%*%exp(CC.CV_k*g_CV_k)))[d_k==1] )
      l_CV_k=l_CV_k-S0
      CVL_CV=CVL_CV+l.func(g_CV_k)-as.numeric(l_CV_k)
    }
    c.index_CV=unname( survConcordance(  Surv(t.vec,d.vec)~CC.CV  )$concordance )

    Gene_list=data.frame(gene=Symbol[order(P)],
                         beta=round(beta_est[order(P)],3),
                         Z=Z[order(P)],P=P[order(P)])
    #print(Gene_list)
    
    list(gene=Symbol[order(P)],
         beta=round(beta_est[order(P)],3),
         Z=Z[order(P)],P=P[order(P)],
         c_index=round(c("CC not cross-validated"=c.index,
                         "CC cross-validated"=c.index_CV),3),
         CVL=c("CC not cross-validated"=CVL,
               "CC cross-validated"=CVL_CV))
  }
  
}
