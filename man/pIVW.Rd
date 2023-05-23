\name{mr.pivw}
\alias{mr.pivw}
\title{The penalized inverse-variance weighted estimator for Mendelian randomization}
\usage{
  mr.pivw(data, lambda=1, plei=TRUE, sel.pval=NULL, delta=0, Boot.Fieller=TRUE, sig=0.05)
}
\description{
  The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.}

\arguments{
  \item{data}{A matrix or data frame consists of four columns. The 1st and 2nd columns contain the SNP effects on the outcome and the exposure, respectively. The 3rd and 4th columns contain the standard errors of the SNP effects on the outcome and the exposure, respectively.}
  \item{lambda}{The penalty parameter in the pIVW estimator. The penalty parameter plays a role in the bias-variance trade-off of the estimator. It is recommended to choose lambda=1 to achieve the smallest bias and valid inference. By default, \code{lambda=1}.}
  \item{plei}{If \code{plei=TRUE}, then the horizontal pleiotropy will be taken into account in the pIVW estimator. By default, \code{plei=TRUE}.}
  \item{sel.pval}{A vector containing the P values of the SNP effects on the exposure, which will be used for the IV selection. \code{sel.pval} should be provided when \code{delta} is not zero.}
  \item{delta}{The z-score threshold for IV selection. By default, \code{delta=0} (i.e., no IV selection will be conducted).}
  \item{Boot.Fieller}{If \code{Boot.Fieller=TRUE}, then the P value and the confidence interval of the causal effect based on the bootstrapping Fieller method will be calculated. It is recommended to use the bootstrapping Fieller method when \code{eta} (the estimated effective sample size) is smaller than 10. By default, \code{Boot.Fieller=TRUE}. }
  \item{sig}{The 100(1-\code{sig})\% confidence interval of the causal effect is calculated. By default, \code{sig=0.05}.}
}

\value{
  \item{beta.hat}{The estimated causal effect of the exposure on the outcome}
  \item{beta.se}{The estimated standard error of \code{beta.hat}}
  \item{pval (Normal)}{The P value for testing whether the causal effect is zero, which is derived from the normal distribution.}
  \item{CI (Normal)}{The confidence interval of the causal effect derived from the normal distribution.}
  \item{pval (Bootstrap Fieller)}{The P value for testing whether the causal effect is zero, which is derived from the bootstrapping Fieller method.}
  \item{CI (Bootstrap Fieller)}{The confidence interval of the causal effect derived from the bootstrapping Fieller method.}
  \item{tau2}{The variance of the horizontal pleiotropy. \code{tau2} is calculated by using all IVs in the data before conducting the IV selection.}
  \item{eta}{The estimated effective sample size. It is recommended to be greater than 5 for pIVW to achieve reliable asymptotic properties.}
}

\references{
  Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. <doi:10.1111/biom.13732>
}
\examples{
  mr.pivw(data=example)
}
