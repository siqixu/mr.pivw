#' PIVW Class
#'
#' @description An object containing the estimate produced using the penalized inverse-variance weighted (pIVW) method as well as various statistics.
#'
#' @slot Over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @slot Boot.Fieller \code{TRUE} if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, \code{FALSE} otherwise.
#' @slot N.boot The number of bootstrap samples used in the bootstrapping Fieller method.
#' @slot Lambda The penalty parameter in the pIVW estimator.
#' @slot Delta The z-score threshold for IV selection.
#' @slot Estimate The causal point estimate from the pIVW estimator.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the pIVW estimator.
#' @slot Tau2 The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated by using all IVs in the data before conducting the IV selection.
#' @slot SNPs The number of SNPs after IV selection.
#' @slot Condition The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.
#'
#' @keywords internal


setClass("PIVW",
         representation(Over.dispersion = "logical",
                        Boot.Fieller = "logical",
                        N.boot = 'numeric',
                        Lambda = "numeric",
                        Delta = "numeric",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        Tau2 = "numeric",
                        SNPs = "numeric",
                        Condition = "numeric")
)















