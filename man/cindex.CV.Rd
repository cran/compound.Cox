\name{cindex.CV}
\alias{cindex.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validated c-index for association parameter under a copula-based dependent censoring model.
}
\description{
This function calculates the cross-validated c-index under a copula-based dependent censoring model proposed in
Emura and Chen (2014). The dependency between the failure and censoring times is modeled via the Clayton copula.
}
\usage{
cindex.CV(t.vec, d.vec, X.mat, alpha, K = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{alpha}{Association parameter; Kendall's tau = alpha/(alpha+2)}
  \item{K}{The number of cross-validation folds}
}
\details{
Currently, only the Clayton copula is implemented for association models. The Clayton model yields
positive association between failure and censoring times with the Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0. The independent copula corresponds to alpha = 0.

If the number of covariates p is large (p>=100), the computational time becomes very long. 
Pre-filtering is recommended to reduce p. 
}
\value{
\item{concordant }{Cross-validated c-index}
}

\references{
Emura T and Chen YH (2014) Gene selection for survival data under dependent censoring: 
a copula-based approach, Statistical Methods in Medical Research, in press.
}
\author{Takeshi Emura}

\examples{

n=25 ### sample size ###
p=3  ### the number of covariates ###
set.seed(1)
T=rexp(n) ### survival time
U=rexp(n) ### censoring time
t.vec=pmin(T,U) ### minimum of survival and censoring times
d.vec=as.numeric( c(T<=U) ) ### censoring indicator
X.mat=matrix(runif(n*p),n,p) ### covariates matrix

cindex.CV(t.vec,d.vec,X.mat,alpha=2) ### alpha=2 corresponds to Kendall's tau=0.5

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
