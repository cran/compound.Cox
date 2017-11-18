uni.selection=function(t.vec, d.vec, X.mat, P.value=0.001,K=10,score=FALSE,d0=0,randomize=FALSE){
  
  n=length(t.vec)
  if(randomize==TRUE){
    rand=sample(1:n)  ### randomize patient IDs before cross-validation ###
    t.vec=t.vec[rand]
    d.vec=d.vec[rand]
    X.mat=X.mat[rand,]
  } 
  
  p=ncol(X.mat)
  if(score==TRUE){ res=uni.score(t.vec, d.vec, X.mat,d0) }else{ 
    res=uni.Wald(t.vec, d.vec, X.mat) 
  }
  temp=res$P<P.value
  if(sum(temp)==0){warning("no gene selected; increase P.value")}else{
    
    beta_est=res$beta_est[temp]
    Z=res$Z[temp]
    P=res$P[temp]
    X.cut=as.matrix(X.mat[,temp])
    q=ncol(X.cut)

    if(score==TRUE){ w=Z }else{ w=beta_est }
    CC=X.cut%*%w
    c.index0=unname( survConcordance(  Surv(t.vec,d.vec)~CC  )$concordance )
    
    ##### Log-partial likelihood ######
    atr_t=(matrix(t.vec,n,n,byrow=TRUE)>=matrix(t.vec,n,n))
    l.func=function(g){
      l=sum( (d.vec*CC*g) )-sum( d.vec*log(atr_t%*%exp(CC*g)) )
      as.numeric( l )
    }
   
    ### Cross-validation ###
    CC.CV=CC.test=NULL
    CVL0=CVL1=CVL2=0
    g_kk_vec=numeric(K)
    folds=cut(seq(1,n),breaks=K,labels=FALSE)
    
    for(k in 1:K){
      temp=which(folds==k)
      t_k=t.vec[-temp]
      d_k=d.vec[-temp]
      CC_k=CC[-temp]
      X_k=X.mat[-temp,]
      n_k=length(t_k)
      
      if(score==TRUE){ res=uni.score(t_k, d_k, X_k,d0) }else{ 
        res=uni.Wald(t_k, d_k, X_k)
      }
      temp_k=res$P<P.value
      
      if(sum(temp_k)==0){CC_kk=rep(0,n)}else{
        if(score==TRUE){ w_k=res$Z[temp_k] }else{ w_k=res$beta_est[temp_k] }
        CC_kk=as.matrix(X.mat[,temp_k])%*%w_k
      }
      CC.test=c(CC.test,CC_kk[temp])

      ##### Not cross-validating #####
      res_k=coxph(Surv(t_k,d_k)~CC_k)
      CVL0=CVL0+l.func(res_k$coef)-res_k$loglik[2] 

      ##### Cross-validating only estimation ### 
      if(score==TRUE){ w_CV=uni.score(t_k, d_k, X.cut[-temp,],d0)$Z }else{ 
        w_CV=uni.Wald(t_k, d_k, X.cut[-temp,])$beta_est
      }
      CC.CV=c(CC.CV,X.cut[temp,]%*%as.matrix(w_CV))
      CC.CV_k=X.cut[-temp,]%*%as.matrix(w_CV)
      res_CV_k=coxph(Surv(t_k,d_k)~CC.CV_k)
      CVL1=CVL1+l.func(res_CV_k$coef)-res_CV_k$loglik[2] 
     
      ##### LCV2 (Cross-validating both selection and estimation) #####
      #if(is.null(w_k)){CC_kk=rep(0,n_k)}else{CC_kk=X_k.cut%*%w_k}
      res_kk=coxph(Surv(t_k,d_k)~CC_kk[-temp])
      CVL2=CVL2-as.numeric(res_kk$loglik[2])
      g_kk_vec[k]=res_kk$coef
    }
    
    ##### LCV2 ######
    l_kk.func=function(g){
      l=sum( (CC.test*g)[d.vec==1] )-sum( (log(atr_t%*%exp(CC.test*g)))[d.vec==1] )
      as.numeric( l )
    }
    for(k in 1:K){ CVL2=CVL2+l_kk.func(g_kk_vec[k]) }
    
    ##### c-index (Cross-validating estimation of CC) 
    c.index1=unname( survConcordance(  Surv(t.vec,d.vec)~CC.CV  )$concordance )
    c.index2=unname( survConcordance(  Surv(t.vec,d.vec)~CC.test  )$concordance )
    
    par(mfrow=c(1,2))
    plot(CC,CC.test,xlab="CC (Not cross-validated)",
         ylab="CC (Fully cross-validated)" )
    plot(CC.CV,CC.test,xlab="CC (Only estimation cross-validated)",
         ylab="CC (Fully cross-validated)" )
    
    list(beta=beta_est[order(P)],Z=Z[order(P)],P=P[order(P)],
         c_index=c("Not cross-validated"=c.index0,"Only estimation cross-validated"=c.index1,
               "Both selection and estimation cross-validated"=c.index2),
         CVL=c("Not cross-validated"=CVL0,"Only estimation cross-validated"=CVL1,
               "Both selection and estimation cross-validated"=CVL2)
    )
  }

}
