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

\item{\code{Boot.Fieller}}{\code{TRUE} if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, \code{FALSE} otherwise.}

\item{N.boot}{The number of bootstrap samples used in the bootstrapping Fieller method.}

\item{\code{Lambda}}{The penalty parameter in the pIVW estimator. }

\item{\code{Delta}}{The z-score threshold for IV selection.}

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
