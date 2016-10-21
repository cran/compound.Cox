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
\item{SE }{Standard error for the estimated regression coefficients}
\item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
\item{P }{P-value for testing H_0: beta=0 (Wald test)}
\item{alpha }{Estimated association parameter for the Clayton copula by 
maximizing the cross-validated c-index}
\item{c_index }{Maximized value of the cross-validated c_index}
}
\references{
Emura T and Chen YH (2014) Gene selection for survival data under dependent censoring: 
a copula-based approach, Statistical Methods in Medical Research, doi: 10.1177/0962280214533378.
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

### 63 lung cancer patients (training samples only) used in Section 5 of Emura and Chen (2014) ###
t.vec=c(49.2739274, 20.0660066, 26.9966997, 35.9380152, 45.2805281, 9.4389439,
41.3201320, 40.8580858, 11.2541254, 31.3531353, 20.0330033, 16.4686469, 36.6996700,
36.2706271, 14.4081767, 30.9900990, 32.1452145, 7.8877888, 32.6402640, 29.0759076, 
4.8466864, 26.5346535, 26.9306931, 20.3300330, 29.6369637, 17.9537954, 18.6468647, 
29.6699670, 27.1287129, 27.0297030, 26.3366337, 27.7887789, 27.7942631, 14.6864687,
20.5610561, 26.2706271, 26.8316832, 22.7392739, 10.8910891, 16.8976898, 22.2112211, 
18.2178218, 21.5181518, 2.5082508, 20.8250825, 20.6270627, 2.9372937, 20.1650165,
19.8349835, 18.4488449, 13.8613861, 18.6798680, 13.8613861, 14.1254125, 16.2046205,
0.2310231, 16.4686469, 12.2112211, 4.7194719, 5.2805281, 7.5577558, 10.8250825, 9.6369637)

d.vec=c(0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1,
 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
 0, 0, 1, 1, 0, 0, 0)

### ZNF264 gene expression ###
ZNF264=c(1, 2, 1, 1, 2, 4, 2, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 4,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 4, 3, 1, 3, 2, 1, 1, 4, 3, 1, 2,
 1, 1, 1, 2, 2, 1)
### HGF gene expression ###
HGF=c(1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1,
 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 4, 4, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1,
 1, 1, 2, 2, 1)

X.mat=cbind(ZNF264,HGF)
#dependCox.reg.CV(t.vec,d.vec,X.mat,G=20) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
