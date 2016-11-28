\name{Lung}
\alias{Lung}
\docType{data}
\title{
Survival data for patients with non-small-cell lung cancer.
}
\description{
A subset of the lung cancer data (Chen et al. 2007) is given. 
The subset consists of 97 gene expressions from 125 patients with non-small-cell lung cancer.
The 97 genes were selected with P-value<0.20 under univariate 
Cox regression analyses (Emura et al. 2016, 2017).
}
\usage{data("Lung")}
\format{
  A data frame with 125 observations on the following 100 variables.
  \describe{
    \item{\code{t.vec}}{survival times (time to either death or censoring) in months}
    \item{\code{d.vec}}{censoring indicators, 1=death, 0=censoring}
    \item{\code{train}}{1=training set, 0=testing set, as defined in Chen et al. (2007)}
    \item{\code{VHL}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{IHPK1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{HMMR}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CMKOR1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PLAU}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{IGF2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{FGB}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MYBL2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ODC1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MTHFD2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{GLIPR1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{EZH2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{HCK}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CCNC}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{XRCC1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CYP1B1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CDC25A}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CD44}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{LCK}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MTHFS}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PON3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PTPN6}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{KIDINS220}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{KLHL22}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RBBP6}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{GABARAPL2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SEH1L}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CITED2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{BARD1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{TLX1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CRMP1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CTNNA1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ANXA5}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PTGS2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SMC4L1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{LOC285086}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ATP11B}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CDK10}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{IRF4}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MYH11}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ME3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CCT6A}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SNCG}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MAK3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{VCPIP1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{JMJD1A}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{STAT2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{DDX6}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ERBB3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PAX2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PCTK2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{NF1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{DLG2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{JMJD1A.1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SUCLA2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MMP16}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{AP3B2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{HGF}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MAP2K3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CPEB4}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ZNF264}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{AXL}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CDC23}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MAST3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{COX11}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PRKAG2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MAN1B1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{F8}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RSU1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MMD}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{AK5}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{IDS}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{BNIP1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{ENG}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PCDHGC3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RALY}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{WDR33}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RNF4}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PRDX1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{FXN}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{PTPRU}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{FRAP1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MMP7}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CST3}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{TIMP2}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{TAL1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{STAT1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{CCND1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{DUSP6}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SNRPF}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{MMP13}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{NR2F6}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{HOXA1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RIPK1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{IL7R}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{SEC13L1}}{gene expression, coded as 1, 2, 3, or 4}
    \item{\code{RPL5}}{gene expression, coded as 1, 2, 3, or 4}
  }
}
\details{
Survival data consisting of 125 patients. 
}
\source{
Chen HY, Yu SL, Chen CH, et al (2007). A Five-gene Signature and Clinical Outcome in Non-small-cell Lung Cancer, N Engl J Med 356: 11-20.
}
\references{
Chen HY, Yu SL, Chen CH, et al (2007). A Five-gene Signature and Clinical Outcome in Non-small-cell Lung Cancer, N Engl J Med 356: 11-20.

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57

Emura T, Chen YH, Rondeau V (2017). Survival Analysis with Dependent Censoring and Correlated Endpoints, Copula-Based Approaches
( to be published in JSS Research Series in Statistics, Springer)
}
\examples{
data(Lung)
Lung[1,]
}
\keyword{datasets}
\keyword{gene expression}