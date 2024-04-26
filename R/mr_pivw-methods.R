# Required Functions
#' Generate bootstrap samples for the bootstrapping Fieller's confidence interval of the penalized inverse-variance weighted (pIVW) method
#'
#' @description Internal function of the penalized inverse-variance weighted (pIVW) method, which generates bootstrap samples for the bootstrapping Fieller's confidence interval.
#'
#' @param object An data frame containing Bx, By, Bxse and Byse.
#' @param beta_hat The causal effect estimate.
#' @param tau2 The estimated variance of the horizontal pleiotropy.
#' @param lambda The penalty parameter in the pIVW estimator. By default, \code{lambda=1}.
#' @param n.boot The number of bootstrap samples. By default, \code{n.boot=1000}.
#' @param seed_boot The seed for random sampling in the bootstrap method. By default, \code{seed_boot=1}.
#'
#' @keywords internal
#'
#' @return \item{z_b}{A vector containing the bootstrap samples for the bootstrapping Fieller's confidence interval.}
#'
#' @export

BF_dist = function(object,beta_hat=0,tau2=0,lambda=1,n.boot=1000,seed_boot=1){
  i = 0
  set.seed(seed_boot)
  z_b = numeric()

  while(i < n.boot){
    Bx = object$Bx
    Byse  = object$Byse
    Bxse  = object$Bxse
    alpha = rnorm(length(Bx),0, sqrt(tau2))
    By_hat = rnorm(length(Bx),Bx*beta_hat+alpha, Byse)
    Bx_hat = rnorm(length(Bx),Bx, Bxse)

    v1 = sum(Byse^(-4)*(By_hat^2*Bx_hat^2-(By_hat^2-Byse^2-tau2)*(Bx_hat^2-Bxse^2)))
    v2 = sum(Byse^(-4)*(4*(Bx_hat^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
    v12= sum(2*Byse^(-4)*By_hat*Bx_hat*Bxse^2)
    mu1= sum(Byse^(-2)*By_hat*Bx_hat)
    mu2= sum(Byse^(-2)*(Bx_hat^2-Bxse^2))
    mu2p = mu2/2+sign(mu2)*sqrt(mu2^2/4+lambda*v2)
    mu1p = mu1+(v12/v2)*(mu2p-mu2)
    w = mu2p/(2*mu2p-mu2)
    v1p = v1
    v2p = w^2*v2
    v12p = w*v12

    temp = (mu1p-beta_hat*mu2p)^2/(v1p-2*beta_hat*v12p+beta_hat^2*v2p)
    if(temp<0){
      next
    }else{
      z_b = c(temp,z_b)
      i = i+1
    }
  }
  z_b = sort(z_b)
  return(z_b)
}



#' Penalized Inverse-Variance Weighted (pIVW) Method for Mendelian Randomization
#'
#' @description The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.otherwise.
#'
#' @param Bx A numeric vector of beta-coefficient values for genetic associations with the exposure variable.
#' @param Bxse The standard errors associated with the beta-coefficients \code{Bx}.
#' @param By A numeric vector of beta-coefficient values for genetic associations with the outcome variable.
#' @param Byse The standard errors associated with the beta-coefficients \code{By}.
#' @param lambda The penalty parameter in the pIVW estimator. It plays a role in the bias-variance trade-off of the estimator. It is recommended to choose \code{lambda=1} to achieve the smallest bias and valid inference. By default, \code{lambda=1}.
#' @param over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @param delta The z-score threshold for IV selection. \code{delta} should be greater than or equal to zero. By default, \code{delta=0} (i.e., no IV selection will be conducted).  See 'Details'.
#' @param sel.pval A numeric vector containing the P-values of the SNP effects on the exposure, which will be used for the IV selection. \code{sel.pval} should be provided when \code{delta} is not zero. See 'Details'.
#' @param Boot.Fieller If \code{Boot.Fieller=TRUE}, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. By default, \code{Boot.Fieller=TRUE} when \code{Condition} is smaller than 10 (see 'Details'), and \code{Boot.Fieller=FALSE} otherwise. 
#' @param n.boot The number of bootstrap samples used in the bootstrapping Fieller method. It will be used only when \code{Boot.Fieller=TRUE}. By default, \code{n.boot=1000}. A larger value of \code{n.boot} should be provided when a more precise P-value is needed.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#'
#' @details The penalized inverse-variance weighted (pIVW) estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously
#' in two-sample MR with summary statistics, i.e., an exposure sample (with IV-exposure effect \code{Bx} and standard error \code{Bxse}) and
#' an outcome sample (with IV-outcome effect \code{By} and standard error \code{Byse}).
#'
#' The pIVW estimator also allows for IV selection in three-sample MR, where weak IVs are screened out using
#' an extra sample (with IV-exposure effect \code{Bx*} and standard error \code{Bxse*}) independent of the exposure sample and outcome sample.
#' Generally, the P-value for \code{Bx*} can be computed By \code{sel.pval=2*pnorm(abs(Bx*/Bxse*), lower.tail = FALSE)},
#' Given \code{sel.pval} and a z-score threshold \code{delta}, the variants kept in the analysis will be those
#' with \code{sel.pval<2*pnorm(delta,lower.tail = FALSE)}.
#'
#' The \code{mr_pivw} function outputs a measure \code{Condition} that needs to be large for reliable asymptotic properties of the pIVW estimator.
#' We also refer to \code{Condition} as effective sample size, which is a function of a measure of IV strength and the number of IVs.
#' When \code{delta} is zero (i.e., no IV selection), \code{Condition = (average F-statistic -1)*sqrt(# snps)}. When \code{delta} is not zero
#' (i.e., IV selection is conducted), \code{Condition = [(average F-statistic -1)*sqrt(# snps)]/c},
#' where the numerator is computed using the selected variants, and the denominator \code{c} involves the selection probabilities
#' of all variants (see more details in the paper \url{ https://doi.org/10.1111/biom.13732}). We suggest that \code{Condition} should be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.
#'
#'
#' @return The output from the function is a \code{PIVW} object containing:
#'
#'  \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
#'  \item{Boot.Fieller}{\code{TRUE} if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, \code{FALSE} otherwise.}
#'  \item{N.boot}{The number of bootstrap samples used in the bootstrapping Fieller method.}
#'  \item{Lambda}{The penalty parameter in the pIVW estimator.}
#'  \item{Delta}{The z-score threshold for IV selection.}
#'  \item{Estimate}{The causal point estimate from the pIVW estimator.}
#'  \item{StdError}{The standard error associated with \code{Estimate}.}
#'  \item{CILower}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.}
#'  \item{CIUpper}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.}
#'  \item{Alpha}{The significance level used in constructing the confidence interval.}
#'  \item{Pvalue}{P-value associated with \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution.}
#'  \item{Tau2}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated By using all IVs in the data before conducting the IV selection.}
#'  \item{SNPs}{The number of SNPs after IV selection.}
#'  \item{Condition}{The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties. See 'Details'.}
#'
#' @examples mr_pivw(Bx = Bx_exp, Bxse = Bxse_exp, By = By_exp, Byse = Byse_exp)
#'
#' @references Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. Available at \url{ https://doi.org/10.1111/biom.13732}.
#'
#' @export

mr_pivw = function(Bx, Bxse, By, Byse, lambda=1, over.dispersion=TRUE, delta=0, sel.pval=NULL, Boot.Fieller=NULL, n.boot=1000, alpha=0.05)
{

  if(!all(length(Bx) == length(By),
          length(Bxse) == length(Byse),
          length(Bxse) == length(By))){
    cat("Bx, By, Bxse and Byse do not all have the same length.")
    return()
  }

  object = data.frame(Bx, By, Bxse, Byse)
  Bx = object$Bx
  Bxse  = object$Bxse
  p = length(Bx)

  if(lambda<0){
    cat("\'lambda\' cannot be smaller than zero.","\n")
    return()
  }

  if(delta<0){
    cat("\'delta\' cannot be smaller than zero.","\n")
    return()
  }

  if(alpha<=0){
    alpha = 0.05
    cat("\'alpha\' provided is less than or equal to zero. \'alpha\ is set to be 0.05.","\n")
  }
  if(alpha>=1){
    alpha = 0.05
    cat("\'alpha\' provided is greater than or equal to one. \'alpha\ is set to be 0.05.","\n")
  }

  if(delta>0){
    if(is.null(sel.pval)){
      delta=0
      cat("\'sel.pval\' is not provided. \'delta\' is set to be zero and no IV selection conducted.","\n")
    }else if(length(sel.pval)!=p){
      delta=0
      cat("The length of \'sel.pval\' doesn't match the number of snps in the MRInput object. \'delta\' is set to be zero and no IV selection conducted.","\n")
    }
  }



  if(delta==0){
    kappa = mean(Bx^2/Bxse^2)-1
    eta = kappa*sqrt(p)
  }else{
    sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
    if(sum(sel)==0){
      cat("No snps is kept in the anaysis after IV selection. Please try a smaller \'delta\' value.","\n")
      return()
    }
    kappa = mean(Bx[sel]^2/Bxse[sel]^2)-1
    eta = kappa*sqrt(sum(sel))
    sel.z = qnorm(sel.pval/2,lower.tail = FALSE)
    q = pnorm(sel.z-delta,lower.tail = TRUE)+pnorm(-sel.z-delta,lower.tail = TRUE)
    psi2 = sum(((Bx/Bxse)^4-6*(Bx/Bxse)^2+3)*q*(1-q))/sum(sel)
    eta = eta/max(1,sqrt(psi2))
  }

  if(delta!=0 & over.dispersion!=TRUE){
    sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
    object$By = object$By[sel]
    object$Bx = object$Bx[sel]
    object$Byse = object$Byse[sel]
    object$Bxse = object$Bxse[sel]
  }

  By = object$By
  Bx = object$Bx
  Byse  = object$Byse
  Bxse  = object$Bxse

  v2 = sum(Byse^(-4)*(4*(Bx^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
  v12 = sum(2*Byse^(-4)*By*Bx*Bxse^2)
  mu1 = sum(Byse^(-2)*By*Bx)
  mu2 = sum(Byse^(-2)*(Bx^2-Bxse^2))
  mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
  mu1p = mu1 + (v12/v2)*(mu2p-mu2)
  beta_pIVW = mu1p/mu2p

  if(over.dispersion==TRUE){
    tau2 = sum(((By-beta_pIVW*Bx)^2 - Byse^2 -beta_pIVW^2*Bxse^2) * Byse^(-2))/sum(Byse^(-2))
    if (tau2 < 0){tau2 = 0}
  }else{
    tau2 = 0
  }

  if(delta!=0 & over.dispersion==TRUE){
    sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
    object$By = object$By[sel]
    object$Bx = object$Bx[sel]
    object$Byse = object$Byse[sel]
    object$Bxse = object$Bxse[sel]

    By = object$By
    Bx = object$Bx
    Byse = object$Byse
    Bxse = object$Bxse

    v2 = sum(Byse^(-4)*(4*(Bx^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
    v12 = sum(2*Byse^(-4)*By*Bx*Bxse^2)
    mu1 = sum(Byse^(-2)*By*Bx)
    mu2 = sum(Byse^(-2)*(Bx^2-Bxse^2))
    mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
    mu1p = mu1 + (v12/v2)*(mu2p-mu2)
    beta_pIVW = mu1p/mu2p
  }
  se_pIVW = 1/abs(mu2p)*sqrt(sum((Bx^2/Byse^2)*(1+tau2*Byse^(-2))
                                 + beta_pIVW^2*(Bxse^2/Byse^4)*(Bx^2+Bxse^2)))

  if(is.null(Boot.Fieller)){
    Boot.Fieller = ifelse(eta<10, TRUE, FALSE) 
  }

  if(Boot.Fieller==TRUE){

    z_b = BF_dist(object,beta_pIVW,tau2,lambda,n.boot)
    v1p = sum(Byse^(-4)*(By^2*Bx^2-(By^2-Byse^2-tau2)*(Bx^2-Bxse^2)))
    w = mu2p/(2*mu2p-mu2)
    v2p = w^2*v2
    v12p = w*v12
    z0 = mu1p^2/v1p
    pval = sum(z0<z_b)/length(z_b)

    qt = z_b[round(length(z_b)*(1-alpha))]
    A = mu2p^2-qt*v2p
    B = 2*(qt*v12p - mu1p*mu2p)
    C = mu1p^2 -qt*v1p
    D = B^2-4*A*C

    if(A>0){
      r1 = (-B-sqrt(D))/2/A
      r2 = (-B+sqrt(D))/2/A
      CI = cbind(r1,r2)
      CI_ll = r1
      CI_ul = r2
    }else if (D>0){
      r1 = (-B-sqrt(D))/2/A
      r2 = (-B+sqrt(D))/2/A
      CI1 = c(-Inf,r1)
      CI2 = c(r2,Inf)
      CI = rbind(CI1,CI2)
      CI_ll = c(-Inf,r2)
      CI_ul = c(r1,Inf)
    }else{
      CI = cbind(-Inf,Inf)
      CI_ll = -Inf
      CI_ul = Inf
    }

  }else{
    pval = 2*pnorm(abs(beta_pIVW),0,se_pIVW,lower.tail = FALSE)
    CI_ll = beta_pIVW+qnorm(alpha/2)*se_pIVW
    CI_ul = beta_pIVW-qnorm(alpha/2)*se_pIVW
    CI = cbind(CI_ll,CI_ul)
  }


  return(new("PIVW",
             Over.dispersion = over.dispersion,
             Boot.Fieller = Boot.Fieller,
             N.boot = n.boot,
             Lambda = lambda,
             Delta = delta,

             Estimate = beta_pIVW,
             StdError = se_pIVW,
             CILower = CI_ll,
             CIUpper = CI_ul,
             Alpha = alpha,

             Pvalue = as.numeric(pval),
             Tau2 = tau2,
             SNPs = length(By),
             Condition = eta)
  )




}

