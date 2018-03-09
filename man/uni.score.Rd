\name{uni.score}
\alias{uni.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Univariate Cox score test}
\description{
Univariate significance analyses based on the score tests 
}
\usage{
uni.score(t.vec, d.vec, X.mat, d0=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{d0}{A positive constant to stabilize the variance (Witten & Tibshirani 2010)}
}
\details{score test}
\value{
  \item{beta }{Estimated regression coefficients (one-step estimator)}
  \item{Z }{Z-value for testing H_0: beta=0 (score test)}
  \item{P }{P-value for testing H_0: beta=0 (score test)}
}
\references{
Witten DM, Tibshirani R (2010) Survival analysis with high-dimensional covariates. Stat Method Med Res 19:29-51
}
\author{Takeshi Emura and Shigeyuki Matsui}
\examples{
data(Lung)
t.vec=Lung$t.vec[Lung$train==TRUE]
d.vec=Lung$d.vec[Lung$train==TRUE]
X.mat=Lung[Lung$train==TRUE,-c(1,2,3)]
uni.score(t.vec, d.vec, X.mat)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univariate Cox regression }
\keyword{ score test }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
