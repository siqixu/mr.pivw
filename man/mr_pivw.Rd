\name{mr_pivw}
\alias{mr_pivw}
\title{Penalized Inverse-Variance Weighted (pIVW) Method for Mendelian Randomization}
\usage{
mr_pivw(
  Bx,
  Bxse,
  By,
  Byse,
  lambda = 1,
  over.dispersion = TRUE,
  delta = 0,
  sel.pval = NULL,
  Boot.Fieller = NULL,
  n.boot = 1000,
  alpha = 0.05
)
}
\arguments{
\item{Bx}{A numeric vector of beta-coefficient values for genetic associations with the exposure variable.}

\item{Bxse}{The standard errors associated with the beta-coefficients \code{Bx}.}

\item{By}{A numeric vector of beta-coefficient values for genetic associations with the outcome variable.}

\item{Byse}{The standard errors associated with the beta-coefficients \code{By}.}

\item{lambda}{The penalty parameter in the pIVW estimator. It plays a role in the bias-variance trade-off of the estimator. It is recommended to choose \code{lambda=1} to achieve the smallest bias and valid inference. By default, \code{lambda=1}.}

\item{over.dispersion}{Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.}

\item{delta}{The z-score threshold for IV selection. \code{delta} should be greater than or equal to zero. By default, \code{delta=0} (i.e., no IV selection will be conducted).  See 'Details'.}

\item{sel.pval}{A numeric vector containing the P-values of the SNP effects on the exposure, which will be used for the IV selection. \code{sel.pval} should be provided when \code{delta} is not zero. See 'Details'.}

\item{Boot.Fieller}{If \code{Boot.Fieller=TRUE}, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. By default, \code{Boot.Fieller=TRUE} when \code{Condition} is smaller than 10 (see 'Details'), and \code{Boot.Fieller=FALSE} otherwise .}

\item{n.boot}{The number of bootstrap samples used in the bootstrapping Fieller method. It will be used only when \code{Boot.Fieller=TRUE}. By default, \code{n.boot=1000}. A larger value of \code{n.boot} should be provided when a more precise P-value is needed.}

\item{alpha}{The significance level used to calculate the confidence intervals. The default value is 0.05.}
}
\value{
The output from the function is a \code{PIVW} object containing:

 \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
 \item{Boot.Fieller}{\code{TRUE} if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, \code{FALSE} otherwise.}
 \item{N.boot}{The number of bootstrap samples used in the bootstrapping Fieller method.}
\item{Lambda}{The penalty parameter in the pIVW estimator.}
 \item{Delta}{The z-score threshold for IV selection.}
 \item{Estimate}{The causal point estimate from the pIVW estimator.}
 \item{StdError}{The standard error associated with \code{Estimate}.}
 \item{CILower}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.}
 \item{CIUpper}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.}
 \item{Alpha}{The significance level used in constructing the confidence interval.}
 \item{Pvalue}{P-value associated with \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution.}
 \item{Tau2}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated By using all IVs in the data before conducting the IV selection.}
 \item{SNPs}{The number of SNPs after IV selection.}
 \item{Condition}{The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties. See 'Details'.}
}
\description{
The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.
}
\details{
The penalized inverse-variance weighted (pIVW) estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously
in two-sample MR with summary statistics, i.e., an exposure sample (with IV-exposure effect \code{Bx} and standard error \code{Bxse}) and
an outcome sample (with IV-outcome effect \code{By} and standard error \code{Byse}).

The pIVW estimator also allows for IV selection in three-sample MR, where weak IVs are screened out using
an extra sample (with IV-exposure effect \code{Bx*} and standard error \code{Bxse*}) independent of the exposure sample and outcome sample.
Generally, the P-value for \code{Bx*} can be computed By \code{sel.pval=2*pnorm(abs(Bx*/Bxse*), lower.tail = FALSE)},
Given \code{sel.pval} and a z-score threshold \code{delta}, the variants kept in the analysis will be those
with \code{sel.pval<2*pnorm(delta,lower.tail = FALSE)}.

The \code{mr_pivw} function outputs a measure \code{Condition} that needs to be large for reliable asymptotic properties of the pIVW estimator.
We also refer to \code{Condition} as effective sample size, which is a function of a measure of IV strength and the number of IVs.
When \code{delta} is zero (i.e., no IV selection), \code{Condition = (average F-statistic -1)*sqrt(# snps)}. When \code{delta} is not zero
(i.e., IV selection is conducted), \code{Condition = [(average F-statistic -1)*sqrt(# snps)]/c},
where the numerator is computed using the selected variants, and the denominator \code{c} involves the selection probabilities
of all variants (see more details in the paper \url{ https://doi.org/10.1111/biom.13732}). We suggest that \code{Condition} should be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.
}
\examples{
mr_pivw(Bx = Bx_exp, Bxse = Bxse_exp, By = By_exp, Byse = Byse_exp)

}
\references{
Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. Available at \url{ https://doi.org/10.1111/biom.13732}.
}
