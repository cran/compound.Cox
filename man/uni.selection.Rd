\name{uni.selection}
\alias{uni.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gene selection based on the univariate Cox regression}
\description{
This function perform gene selection using the univariate Cox regression based on survival data with high-dimensional gene expressions (Matsui 2006; Emura and Chen 2016).
}
\usage{
uni.selection(t.vec, d.vec, X.mat, P.value = 0.001, K = 5,score=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{P.value}{Threshold}
  \item{K}{The number of cross-validation folds}
  \item{score}{If TRUE, the score test is performed instead of the Wald test}
}
\details{Predictive ability of the selected genes are evaluated throught cross-validated log-likelihood (CVL) and c-index are computed.}
\value{
  \item{gene }{Gene symbols}
  \item{beta }{Estimated regression coefficients}
  \item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
  \item{P }{P-value for testing H_0: beta=0 (Wald test)}
  \item{c_index }{c-index}
  \item{CVL }{Cross-validated partial likelihood}
}
\references{
Matsui S (2006). Predicting Survival Outcomes Using Subsets of Significant Genes in Prognostic Marker Studies with Microarrays. BMC Bioinformatics: 7:156.

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57
}
\author{Takeshi Emura}
\examples{
data(Lung)
t.vec=Lung$t.vec[Lung$train==TRUE]
d.vec=Lung$d.vec[Lung$train==TRUE]
X.mat=Lung[Lung$train==TRUE,-c(1,2,3)]
uni.selection(t.vec, d.vec, X.mat, P.value=0.05,K=5)
## the outputs reproduce Table 3 of Emura and Chen (2016) ## 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univariate Cox regression }
\keyword{ c-index }
\keyword{ cross-validated partial likelihood }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
