\name{dependCox.reg.CV}
\alias{dependCox.reg.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cox regression under dependent censoring.
}
\description{
This function perform estimation and significance testing for survival data under a copula-based dependent censoring model proposed in
Emura and Chen (2016). The dependency between the failure and censoring times is modeled via the Clayton copula.
The method is based on the semiparametric maximum likelihood estimation, where the association parameter is estimated by maximizing
the cross-validated c-index (see Emura and Chen 2016 for details).
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
\item{SE }{Standard error for the estimated regression coefficients}
\item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
\item{P }{P-value for testing H_0: beta=0 (Wald test)}
\item{alpha }{Estimated association parameter for the Clayton copula by 
maximizing the cross-validated c-index}
\item{c_index }{Maximized value of the cross-validated c_index}
}
\references{
Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57
}
\author{Takeshi Emura}


\examples{
########### Simulated survival data ##############
n=25 ### sample size 
p=3 ### the number of covariates
set.seed(1)
T=rexp(n) ### survival time
U=rexp(n) ### censoring time
t.vec=pmin(T,U) ### minimum of survival and censoring times
d.vec=as.numeric( c(T<=U) ) ### censoring indicator
X.mat=matrix(runif(n*p),n,p) ### covariate matrix

dependCox.reg.CV(t.vec,d.vec,X.mat,G=10)

### Reproduce Section 5 of Emura and Chen (2016) ###
data(Lung)
temp=Lung[,"train"]==TRUE
t.vec=Lung[temp,"t.vec"]
d.vec=Lung[temp,"d.vec"]
X.mat=as.matrix(Lung[temp,-c(1,2,3)])
#dependCox.reg.CV(t.vec,d.vec,X.mat,G=20) ### take time to finish ###
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ univariate Cox regression }% __ONLY ONE__ keyword per line
