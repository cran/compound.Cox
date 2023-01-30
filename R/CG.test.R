CG.test <-
  function(t.vec,d.vec,PI,cutoff=median(PI),
           alpha=2,copula=CG.Clayton,
           S.plot=TRUE,N=10000,mark.time=TRUE){

###### Plot two CG estimators #####
temp=(PI<=cutoff)
t.good=t.vec[temp]
d.good=d.vec[temp]
event.good=sum(d.good)
res.good=copula(t.good,d.good,alpha,S.plot=F)
S.good=res.good$surv
tau.good=res.good$tau

t.poor=t.vec[!temp]
d.poor=d.vec[!temp]
event.poor=sum(d.poor)
res.poor=copula(t.poor,d.poor,alpha,S.plot=F)
S.poor=res.poor$surv
tau.poor=res.poor$tau

if(S.plot==TRUE){
  plot(c(0,sort(t.good)),c(1,S.good),type="s",lwd=2,
     xlim=c(0,max(t.vec)),ylim=c(0,1),
     col="blue",ylab="Survival probability",xlab="Time")
  points(c(0,sort(t.poor)),c(1,S.poor),type="s",lwd=2,col="red")
  if(mark.time==TRUE){
    points(sort(t.good[d.good==0]),
         S.good[d.good[order(t.good)]==0],pch=3,cex=1,col="blue")
    points(sort(t.poor[d.poor==0]),
         S.poor[d.poor[order(t.poor)]==0],pch=3,cex=1,col="red")
  }
}


### mean difference ###
G=1000
t_grid=S.good_grid=S.poor_grid=numeric(G)
Tau=min(max(t.poor),max(t.good)) ## Pepe & Fleming (1989, 1991)
t_grid=Tau*(1:G)/G
for(i in 1:G){
  S.good_grid[i]=c(1,S.good)[sum(c(0,t.good)<=t_grid[i])]
  S.poor_grid[i]=c(1,S.poor)[sum(c(0,t.poor)<=t_grid[i])]
}
MD=mean(S.good_grid-S.poor_grid)
RMST.good=Tau*mean(S.good_grid)
RMST.poor=Tau*mean(S.poor_grid)
polygon(c(t_grid,rev(t_grid)),c(S.poor_grid,rev(S.good_grid)),
        density=40,col="thistle1",border=NA)

#### Permutation test #####
n.good=sum(PI<=cutoff)
n.poor=sum(PI>cutoff)
MD.perm=numeric(N)
set.seed(1)

for(j in 1:N){
  perm=sample(size=n.good,1:(n.good+n.poor),replace=FALSE)
  t.good=t.vec[perm]
  d.good=d.vec[perm]
  S.good=copula(t.good,d.good,alpha,S.plot=F)$surv

  t.poor=t.vec[-perm]
  d.poor=d.vec[-perm]
  S.poor=copula(t.poor,d.poor,alpha,S.plot=F)$surv

  t_grid=S.good_grid=S.poor_grid=numeric(G)
  Tau_perm=min(max(t.poor),max(t.good))
  t_grid=Tau_perm*(1:G)/G
  for(i in 1:G){
    S.good_grid[i]=c(1,S.good)[sum(c(0,t.good)<=t_grid[i])]
    S.poor_grid[i]=c(1,S.poor)[sum(c(0,t.poor)<=t_grid[i])]
  }
  MD.perm[j]=mean(S.good_grid-S.poor_grid)
}
P_MD=mean( abs(MD.perm)>abs(MD) )

abline(v=Tau,col=rgb(0,0,0, alpha=0.1),cex=2)
list(test=c(Survival.diff=MD,RMSTD=Tau*MD,P.value=P_MD),
     Good=c(n=n.good,n.event=event.good,
            RMST=RMST.good,Kendall.tau=tau.good,meanPI=mean(PI[PI<=cutoff])),
     Poor=c(n=n.poor,n.event=event.poor,
            RMST=RMST.poor,Kendall.tau=tau.poor,meanPI=mean(PI[PI>cutoff])))
}
