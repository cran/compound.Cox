\name{dependCox.reg.CV}
\alias{dependCox.reg.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cox regression under dependent censoring.
}
\description{
This function perform estimation and significance testing for survival data under a copula-based dependent censoring model proposed in
Emura and Chen (2014). The dependency between the failure and censoring times is modeled via the Clayton copula.
The method is based on the semiparametric maximum likelihood estimation, where the association parameter is estimated by maximizing
the cross-validated c-index (see Emura and Chen 2004 for details).
}
\usage{
dependCox.reg.CV(t.vec, d.vec, X.mat, K = 5, G = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{K}{The number of cross-validation folds}
  \item{G}{The number of grids in searching for the maximam of the cross-validated c-index}
}
\details{
Currently, only the Clayton copula is implemented for association models. The Clayton model yields
positive association between failure and censoring times with the Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0. The independent copula corresponds to alpha = 0.

If the number of covariates p is large (p>=100), the computational time becomes very long. 
Pre-filtering is recommended to reduce p. 
}
\value{
\item{beta_hat }{Estimated regression coefficients}
\item{SD }{Standard deviation for the estimated regression coefficients}
\item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
\item{P }{P-value for testing H_0: beta=0 (Wald test)}
\item{alpha }{Estimated association parameter for the Clayton copula by 
maximizing the cross-validated c-index}
}
\references{
Emura T and Chen YH (2014) Gene selection for survival data under dependent censoring: 
a copula-based approach, Statistical Methods in Medical Research, doi: 10.1177/0962280214533378.
}
\author{Takeshi Emura}


\examples{
########### Example ##############

n=25 ### sample size 
p=3 ### the number of covariates
set.seed(1)
T=rexp(n) ### survival time
U=rexp(n) ### censoring time
t.vec=pmin(T,U) ### minimum of survival and censoring times
d.vec=as.numeric( c(T<=U) ) ### censoring indicator
X.mat=matrix(runif(n*p),n,p) ### covariate matrix

dependCox.reg.CV(t.vec,d.vec,X.mat,G=10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
