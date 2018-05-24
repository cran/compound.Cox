\name{uni.selection}
\alias{uni.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Univariate feature selection based on univariate significance tests}
\description{
This function perform univariate feature selection using univariate significance tests, 
where the Wald statistics or score statistics are used to measure significance. 
Features are selected according to whether their P-values are less than a given threshold by user. 
The cross-validated likelihood (CVL) value is computed for selected features (Matsui 2006; Emura et al. 2018).
}
\usage{
uni.selection(t.vec, d.vec, X.mat, P.value = 0.001, K = 10,
score=FALSE,d0=0,randomize=FALSE,c.plot=TRUE,permutation=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{P.value}{A threshold for selecting features}
  \item{K}{The number of cross-validation folds}
  \item{score}{If TRUE, the score test is performed instead of the Wald test}
  \item{d0}{A positive constant to stabilize the variance (Witten & Tibshirani 2010)}
  \item{randomize}{If TRUE, randomize patient ID's before cross-validation}
  \item{c.plot}{If TRUE, the plot of c-index is displayed}
  \item{permutation}{If TRUE, the FDR is computed by randomly permutating the gene expressions}
}
\details{Predictive ability of the selected genes are evaluated throught cross-validated log-likelihood (CVL) and c-index are computed.}
\value{
  \item{gene }{Gene symbols}
  \item{beta }{Estimated regression coefficients}
  \item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
  \item{P }{P-value for testing H_0: beta=0 (Wald test)}
  \item{c_index }{c-index}
  \item{CVL }{Cross-validated partial likelihood}
  \item{Genes }{The number of genes, the number of selected genes, and the number of falsely selected genes}
  \item{FDR }{False Discovery Rate}
}
\references{
Matsui S (2006). Predicting Survival Outcomes Using Subsets of Significant Genes in Prognostic Marker Studies with Microarrays. BMC Bioinformatics: 7:156.

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57

Witten DM, Tibshirani R (2010) Survival analysis with high-dimensional covariates. Stat Method Med Res 19:29-51
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
