\name{dependCox.reg}
\alias{dependCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate Cox regression under dependent censoring.
}
\description{
This function performs the univariate Cox regression under dependent censoring based on Chen (2010) and Emura and Chen (2014).
 The association between survival time and censoring time is modeled via the Clayton copula.
}
\usage{
dependCox.reg(t.vec, d.vec, X.vec, alpha, var = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.vec}{Vector of univariate covariates}
  \item{alpha}{Association parameter; Kendall's tau = alpha/(alpha+2)}
  \item{var}{If TRUE, the standard deviations are given (use FALSE to reduce the computational cost)}
}
\details{
Currently, only the Clayton copula is implemented for association models. The Clayton model yields
positive association between failure and censoring times with the Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0. The independent copula corresponds to alpha = 0. 
}
\value{
\item{beta_hat }{Estimated regression coefficients}
\item{SE }{Standard error for the estimated regression coefficients}
\item{Z }{Z-value for testing H_0: beta=0 (Wald test)}
\item{P }{P-value for testing H_0: beta=0 (Wald test)}
}
\references{
Chen YH (2010). Semiparametric Marginal Regression Analysis for Dependent Competing Risks under an Assumed Copula,
 JRSSB 72: 235-251.

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.
}
\author{Takeshi Emura}


\examples{
### Univariate Cox regression under dependent censoring ###
n=150
beta_true=1.5
t.vec=d.vec=X.vec=numeric(n)
alpha_true=2 ### Kendall's tau=0.5 ##
set.seed(1)

#### Simulated survival data under the Clayton copula model ####
for(i in 1:n){
  X.vec[i]=runif(1)
  eta=X.vec[i]*beta_true
  U=runif(1)
  V=runif(1)
  T=-1/exp(eta)*log(1-U)
  W=(1-U)^(-alpha_true)
  C=1/alpha_true/exp(eta)*log(  1-W+W*(1-V)^(-alpha_true/(alpha_true+1))  )
  t.vec[i]=min(T,C)
  d.vec[i]=(T<=C) 
}


#dependCox.reg(t.vec,d.vec,X.vec,alpha=0)## same result as "coxph(Surv(t.vec,d.vec)~X.vec)" ##
dependCox.reg(t.vec,d.vec,X.vec,alpha=alpha_true)## analysis under dependent censoring ##
dependCox.reg(t.vec,d.vec,X.vec,alpha=alpha_true,var=FALSE)## faster computation by "var=FALSE"##
beta_true


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

d.vec=c(0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
 0, 0, 0)

### ZNF264 gene expression ###
ZNF264=c(1, 2, 1, 1, 2, 4, 2, 1, 1, 1, 4, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 4, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 4, 3, 1, 3, 2, 1, 1, 4, 3, 1, 2, 1, 1,
 1, 2, 2, 1)
#dependCox.reg(t.vec,d.vec,ZNF264,alpha=18)

### HGF gene expression ###
HGF=c(1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1, 1,
 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 4, 4, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
 2, 2, 1)
#dependCox.reg(t.vec,d.vec,HGF,alpha=18)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
