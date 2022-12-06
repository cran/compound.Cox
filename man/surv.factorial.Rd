\name{surv.factorial}
\alias{surv.factorial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factorial survival analysis under dependent censoring
}
\description{
Perform factorial survival analysis under dependent censoring under an assumed copula.
}
\usage{
surv.factorial(t.vec,d.vec,group,alpha,copula,t.upper=max(t.vec),C=NULL,S.plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{group}{Vector of group indicators, 1, 2, ..., d}
  \item{alpha}{Copula parameter}
  \item{copula}{Copula function: "CG.Clayton","CG.Gumbel" or "CG.Frank"}
  \item{t.upper}{Follow-up end (default is max(t.vec))}
  \item{C}{Contrast matrix}
  \item{S.plot}{If TRUE, the survival curve is displayed}
}
\details{
Estimates of treatment effects and the test results are shown.
}
\value{
\item{copula.parameter}{Copula parameter}
\item{p}{Estimates of treatment effects}
\item{Var}{Variance estimates}
\item{c_simulated}{Critical value via simulation method}
\item{c_analytical}{Critical value via analytical method}
}
\references{
Emura T, Matsui S, Chen HY (2019). compound.Cox: Univariate Feature Selection and Compound Covariate for Predicting Survival,
Computer Methods and Programs in Biomedicine 168: 21-37.

Emura T, Chen YH (2018). Analysis of Survival Data with Dependent Censoring,
Copula-Based Approaches, JSS Research Series in Statistics, Springer, Singapore.

Rivest LP, Wells MT (2001). A Martingale Approach to the Copula-graphic Estimator for the
Survival Function under Dependent Censoring, J Multivar Anal; 79: 138-55.
}
\author{Takeshi Emura}

\examples{
## to be written ##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ copula-graphic estimator }
