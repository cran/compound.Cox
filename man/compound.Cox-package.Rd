\name{compound.Cox-package}
\alias{compound.Cox-package}
\alias{compound.Cox}
\docType{package}
\title{
Estimation and Gene Selection Based on the Compound Covariate Method Under the Cox Proportional Hazard Model.  
}
\description{
Calculate regression coefficients and their standard deviations under the Cox proportional hazard model with a large number of covariates.
A shrinkage type estimation based on the compound covariate method (Emura, Chen & Chen, 2012 PLoS ONE) and a univariate gene selection method under dependent 
censoring (Emura & Chen, 2016 Stat Methods Med Res) are implemented in this package. 
Algorithms for simulating gene expressions in these references are also given. A subset of Lung cancer data (Chen et al. 2007) is provided.
}
\details{
\tabular{ll}{
Package: \tab compound.Cox\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2016-11-28\cr
License: \tab GPL-2\cr
}
}

\author{
Takeshi Emura, Hsuan-Yu Chen, Yi-Hau Chen;
Maintainer: Takeshi Emura <takeshiemura@gmail.com>
}
\references{
Chen HY, Yu SL, Chen CH, et al (2007). A Five-gene Signature and Clinical Outcome in Non-small-cell Lung Cancer, N Engl J Med 356: 11-20.

Emura T, Chen YH, Chen HY (2012). Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models.
 PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57
}
\keyword{ package }
