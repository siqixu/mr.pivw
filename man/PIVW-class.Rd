\docType{class}
\name{PIVW-class}
\alias{PIVW-class}
\title{PIVW Class}
\description{
An object containing the estimate produced using the penalized inverse-variance weighted (pIVW) method as well as various statistics.
}
\section{Slots}{

\describe{
\item{\code{Over.dispersion}}{Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.}

\item{\code{Boot.Fieller}}{If \code{Boot.Fieller=TRUE}, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. It is recommended to use the bootstrapping Fieller method when \code{Condition} (the estimated effective sample size) is smaller than 10. By default, \code{Boot.Fieller=TRUE}.}

\item{\code{Lambda}}{The penalty parameter in the pIVW estimator. The penalty parameter plays a role in the bias-variance trade-off of the estimator. It is recommended to choose \code{lambda=1} to achieve the smallest bias and valid inference. By default, \code{lambda=1}.}

\item{\code{Delta}}{The z-score threshold for IV selection. By default, \code{delta=0} (i.e., no IV selection will be conducted).}

\item{\code{Estimate}}{The causal point estimate from the pIVW estimator.}

\item{\code{StdError}}{The standard error associated with \code{Estimate}.}

\item{\code{CILower}}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.}

\item{\code{CIUpper}}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.}

\item{\code{Alpha}}{The significance level used in constructing the confidence interval (default is 0.05).}

\item{\code{Pvalue}}{P-value associated with the causal estimate from the pIVW estimator.}

\item{\code{Tau2}}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated by using all IVs in the data before conducting the IV selection.}

\item{\code{SNPs}}{The number of SNPs after IV selection.}

\item{\code{Condition}}{The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.}
}}

\keyword{internal}
