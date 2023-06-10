# Penalized Inverse-Variance Weighted (pIVW) Method for Mendelian Randomization

The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.

## Setup
Use the following command in R to install the package (which is the latest version):
```
library(devtools)
install_github("siqixu/mr.pivw",ref="main") 
```
Or install the "mr.pivw" package from R CRAN:
```
install.packages("mr.pivw")
```

## Usage
```
mr_pivw(Bx, Bxse, By, Byse, lambda = 1, over.dispersion = TRUE, delta = 0, sel.pval = NULL, Boot.Fieller = TRUE, alpha = 0.05)
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
If `Boot.Fieller=TRUE`, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. It is recommended to use the bootstrapping Fieller method when `Condition` is smaller than 10. By default, `Boot.Fieller=TRUE`.

`alpha`: 	
The significance level used to calculate the confidence intervals. The default value is 0.05.

## Value
`Estimate`: The causal point estimate from the pIVW estimator.	

`StdError`: The standard error associated with `Estimate`.	

`CILower`: The lower bound of the confidence interval for `Estimate`.	

`CIUpper`: The upper bound of the confidence interval for `Estimate`. 	

`Pvalue`: P-value associated with `Estimate`.	

`Condition`: The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.	


## Example 
```
library(mr.pivw)  # load the mr.pivw package
mr_pivw(Bx = Bx_exp, Bxse = Bxse_exp, By = By_exp, Byse = Byse_exp)  # analyze the example data with the pIVW method. 
```
# results 
Penalized inverse-variance weighted method

Over dispersion: TRUE 
Bootstrapping Fieller: TRUE 
Penalty parameter (lambda): 1 
IV selection threshold (delta): 0 
Number of variants : 1000 
------------------------------------------------------------------
 Method Estimate Std Error  95% CI       p-value Condition
   pIVW    0.560     0.291 -0.025, 1.272   0.059    10.071
------------------------------------------------------------------

## Reference
Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. <doi:10.1111/biom.13732>

