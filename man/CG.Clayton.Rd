\name{CG.Clayton}
\alias{CG.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula-graphic estimator under the Clayton copula.
}
\description{
The copula-graphic estimator calculated based on the formula of Rivest & Wells (2001).
We used this estimator in Emura & Chen (2016). 
}
\usage{
CG.Clayton(t.vec, d.vec, alpha, S.plot = TRUE, S.col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{alpha}{Association parameter that is related to Kendall's tau through "tau= alpha/(alpha+2)"}
  \item{S.plot}{If TRUE, plot the survival curve}
  \item{S.col}{Color of the survival curve}
}
\details{
Estimates for survival probability are calculated at given time points of "t.vec". Association parameter "alpha" of the Clayton copula must be given (alpha>0), where alpha=0 corresponds to the independence copula.
}
\value{
\item{time }{sort(t.vec)}
\item{surv }{survival probability at "time"}
}
\references{
Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.

Rivest LP, Wells MT (2001). A Martingale Approach to the Copula-graphic Estimator for the 
Survival Function under Dependent Censoring, J Multivar Anal; 79: 138-55.
}
\author{Takeshi Emura}

\examples{
t.vec=c(1,3,5,4,7,8,10,13)
d.vec=c(1,0,0,1,1,0,1,0)
CG.Clayton(t.vec,d.vec,alpha=18,S.col="blue")
### CG.Clayton gives identical results with the Kaplan-Meier estimator with alpha=0 ### 
CG.Clayton(t.vec,d.vec,alpha=0.000000001,S.plot=FALSE)$surv
survfit(Surv(t.vec,d.vec)~1)$surv
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
