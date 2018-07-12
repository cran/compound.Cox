\name{uni.Wald}
\alias{uni.Wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Univariate Cox Wald test}
\description{
Univariate significance analyses based on the Wald tests 
}
\usage{
uni.Wald(t.vec, d.vec, X.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
}
\details{Wald test}
\value{
  \item{beta }{Estimated regression coefficients}
  \item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
  \item{P }{P-value for testing H_0: beta=0 (Wald test)}
}
\references{
Witten DM, Tibshirani R (2010) Survival analysis with high-dimensional covariates. Stat Method Med Res 19:29-51

Emura T, Matsui S, Chen HY (2018-). compound.Cox: Univariate Feature Selection and Compound Covariate for Predicting Survival, 
Computer Methods and Programs in Biomedicine, in revision.
}
\author{Takeshi Emura and Shigeyuki Matsui}
\examples{
data(Lung)
t.vec=Lung$t.vec[Lung$train==TRUE]
d.vec=Lung$d.vec[Lung$train==TRUE]
X.mat=Lung[Lung$train==TRUE,-c(1,2,3)]
uni.Wald(t.vec, d.vec, X.mat)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univariate Cox regression }
\keyword{ Wald test }
\keyword{ compound covariate }
\keyword{ feature selection }% __ONLY ONE__ keyword per line
