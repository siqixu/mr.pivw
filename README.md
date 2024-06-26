# Penalized Inverse-Variance Weighted (pIVW) Method for Mendelian Randomization

The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.

## Setup
Use the following command in R to install the package (The latest version v0.1.3 updated on April 26, 2024):
```
library(devtools)
install_github("siqixu/mr.pivw",ref="main") 
```

## Usage
```
mr_pivw(Bx, Bxse, By, Byse, lambda = 1, over.dispersion = TRUE, delta = 0, sel.pval = NULL, Boot.Fieller = NULL, n.boot = 1000, alpha = 0.05)
```
`Bx`: A numeric vector of beta-coefficient values for genetic associations with the exposure variable.

`Bxse`: The standard errors associated with the beta-coefficients `Bx`.

`By`: A numeric vector of beta-coefficient values for genetic associations with the outcome variable.

`Byse`: The standard errors associated with the beta-coefficients `By`.	

`lambda`:	The penalty parameter in the pIVW estimator. It plays a role in the bias-variance trade-off of the estimator. It is recommended to choose `lambda=1` to achieve the smallest bias and valid inference. By default, `lambda=1`.

`over.dispersion`: Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.

`delta`: The z-score threshold for IV selection. `delta` should be greater than or equal to zero. By default, `delta=0` (i.e., no IV selection will be conducted). 

`sel.pval`: A numeric vector containing the P-values of the SNP effects on the exposure, which will be used for the IV selection. `sel.pval` should be provided when `delta` is not zero. 

`Boot.Fieller`: 	
If `Boot.Fieller=TRUE`, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. By default, `Boot.Fieller=TRUE` when `Condition` is smaller than 10 (see 'Details' in R Help Documentation), and `Boot.Fieller=FALSE` otherwise.

`n.boot`:
The number of bootstrap samples used in the bootstrapping Fieller method. It will be used only when `Boot.Fieller=TRUE`. By default, `n.boot=1000`. A larger value of `n.boot` should be provided when a more precise P-value is needed.

`alpha`: 	
The significance level used to calculate the confidence intervals. The default value is 0.05.

## Value

`Over.dispersion`: `TRUE` if the method has considered balanced horizontal pleiotropy, `FALSE` otherwise.

`Boot.Fieller`: `TRUE` if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, `FALSE` otherwise.

`N.boot`:	The number of bootstrap samples used in the bootstrapping Fieller method.

`Lambda`: The penalty parameter in the pIVW estimator.

`Delta`: The z-score threshold for IV selection.

`Estimate`: The causal point estimate from the pIVW estimator.	

`StdError`: The standard error associated with `Estimate`.	

`CILower`: The lower bound of the confidence interval for `Estimate`.	

`CIUpper`: The upper bound of the confidence interval for `Estimate`. 	

`Pvalue`: P-value associated with `Estimate`.	

`Tau2`:	The variance of the balanced horizontal pleiotropy. `Tau2` is calculated by using all IVs in the data before conducting the IV selection.

`SNPs`: The number of SNPs after IV selection.

`Condition`: The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.	


## Example 
```
library(mr.pivw)  # load the mr.pivw package
mr_pivw(Bx = Bx_exp, Bxse = Bxse_exp, By = By_exp, Byse = Byse_exp)  # analyze the example data with the pIVW method. 

# results 
Penalized inverse-variance weighted method

Account for over-dispersion: TRUE 
CI and P-value: Normal approximation 
Penalty parameter (lambda): 1 
IV selection threshold (delta): 0 
Number of variants: 1000 

-----------------------------------------------------------
Method Estimate Std Error  95% CI       p-value Condition
  pIVW    0.560     0.291 -0.011, 1.131  0.0545    10.071
-----------------------------------------------------------  
```



## Reference
Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. <doi:10.1111/biom.13732>

