\name{compound.reg}
\alias{compound.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compound covariate methed under the Cox proportional hazard model. 
}
\description{
This function utilizes the so-called compound covariate methed to calculate the regression coefficients
 under the Cox proportional hazard model, which is proposed by Emura, Chen & Chen (see References).
 The method refines the Cox partial likelihood estimator such that
 the regression coefficients are obtained even when the number of covariates is greater than the number of
 sample. The standard errors are also calculated.    
}
\usage{
compound.reg(t.vec, d.vec, X.mat, K = 5, delta_a = 0.025, a_0 = 0, var = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{Design matrix of n by p, where n is the sample size and p is the number of covariates}
  \item{K}{The number of cross validation folds, K=n corresponds to a leave-one-out cross validation (default=5)}
  \item{delta_a}{The step size for a grid search for the maximum of the cross-validated likelihood (default=0.025)}
  \item{a_0}{The starting value of a grid search for the maximum of the cross-validated likelihood (default=0)}
  \item{var}{If TRUE, the standard errors are given (default=FALSE, to reduce the computational cost)}
}

\details{
K=5 cross validation is recommended for computational efficiency, though the results appear to be robust
against the choice of closs validation folds.
}

\value{
\item{a_hat }{Estimated value of the shrinkage parameter}
\item{beta_hat }{Estimated regression coefficients}
\item{SE }{Standard errors for the estimated regression coefficients}
\item{Lower95CI }{Lower ends of 95 percent confidence intervals}
\item{Upper95CI }{Upper ends of 95 percent confidence intervals}
}

\references{
T. Emura, Y.-H. Chen, H.-Y. Chen, Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models with Microarrays,
submitted to PloS ONE
}

\author{Takeshi Emura, Yi-Hau Chen}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n=100 ### sample size
beta_true=c(rep(1,5),rep(0,5))
p=length(beta_true) ### the number of covariates
t.vec=d.vec=numeric(n)
X.mat=matrix(0,n,p)

for(i in 1:n){
  X.mat[i,]=rnorm(p,mean=0,sd=1)
  eta=sum( as.vector(X.mat[i,])*beta_true )
  T=rexp(1,rate=exp(eta))
  C=runif(1,min=0,max=5)
  t.vec[i]=min(T,C)
  d.vec[i]=(T<=C)
}

mean(1-d.vec) ### censored percentage ### 

res=compound.reg(t.vec,d.vec,X.mat)
res
#### compare the estimate (beta_hat) with the true value ###
beta_true


### PBC data analysis for the first 50 patients ###
data(PBC)
t.vec=PBC[1:50,"T"]
d.vec=PBC[1:50,"d"]
X.mat=as.matrix(PBC[1:50,-c(1,2)])
compound.reg(t.vec,d.vec,X.mat,K=5,delta_a=0.10,a_0=0.5,var=TRUE)
# delta_a=0.10 and a_0=0.5 are used to reduce computational time in this example #
# in general, default values of delta_a=0.025 and a_0=0 are recommended #


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ Compound covariate method }% __ONLY ONE__ keyword per line
