\name{dependCox.reg}
\alias{dependCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate Cox regression under dependent censoring.
}
\description{
This function performs the univariate Cox regression under dependent censoring based on Chen (2010).
 The dependency between the failure and censoring times is modeled via the Clayton copula.
}
\usage{
dependCox.reg(t.vec, d.vec, X.vec, alpha, var = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.vec}{Vector of univariate covariates}
  \item{alpha}{Association parameter; Kendall's tau = alpha/(alpha+2)}
  \item{var}{If TRUE, the standard deviations are given (use FALSE to reduce the computational cost)}
}
\details{
Currently, only the Clayton copula is implemented for association models. The Clayton model yields
positive association between failure and censoring times with the Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0. The independent copula corresponds to alpha = 0. 
}
\value{
\item{beta_hat }{Estimated regression coefficients}
\item{SD }{Standard deviation for the estimated regression coefficients}
\item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
\item{P }{P-value for testing H_0: beta=0 (Wald test)}
}
\references{
Chen Y.-H., 2010. Semiparametric marginal regression analysis for dependent competing risks under an assumed copula,
 Journal of the Royal Statistical Society, Ser B 72, 235-251.
}
\author{Takeshi Emura}


\examples{
### Univariate Cox regression under dependent censoring ###
n=150
beta_true=1.5
t.vec=d.vec=X.vec=numeric(n)
alpha_true=2 ### Kendall's tau=0.5 ##
set.seed(1)

#### Generate data under the Clayton copula ####
for(i in 1:n){
  X.vec[i]=runif(1)
  eta=X.vec[i]*beta_true
  U=runif(1)
  V=runif(1)
  T=-1/exp(eta)*log(1-U)
  W=(1-U)^(-alpha_true)
  C=1/alpha_true/exp(eta)*log(  1-W+W*(1-V)^(-alpha_true/(alpha_true+1))  )
  t.vec[i]=min(T,C)
  d.vec[i]=(T<=C) 
}


dependCox.reg(t.vec,d.vec,X.vec,alpha=0)## same result as "coxph(Surv(t.vec,d.vec)~X.vec)" ##
dependCox.reg(t.vec,d.vec,X.vec,alpha=alpha_true)## analysis under dependent censoring ##
dependCox.reg(t.vec,d.vec,X.vec,alpha=alpha_true,var=FALSE)## faster computation by "var=FALSE"##
beta_true
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line