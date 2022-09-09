# A penalized inverse-variance weighted (pIVW) estimator for Mendelian randomization accounting for weak instruments and balanced horizontal pleiotropy

The penalized inverse-variance weighted (pIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The pIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.

## Setup
Use the following command in R to install the package:
```
library(devtools)
install_github("siqixu/mr.pivw",ref="main") # install the "mr.pivw" package
```
Or install the "mr.pivw" package from R CRAN:
```
install.packages("mr.pivw")
```

## Usage
```
mr.pivw(data,lambda=1,plei=TRUE,sel.pval=NULL,delta=0,Boot.Fieller=TRUE,sig=0.05)
```
`data`: A matrix or data frame consists of four columns: the 1st (2nd) column contains the estimated genetic effects on the outcome (exposure); the 3rd (4th) column contains the estimated standard errors of the estimated genetic effects on the outcome (exposure).

`lambda`: The penalty parameter in the pIVW estimator. The penalty parameter plays a role in the bias-variance trade-off of the estimator. We recommend to choose lambda to be 1 to achieve the smallest bias and valid inference. By default, lambda=1.

`plei`: If plei=TRUE, then the horizontal pleiotropy will be taken into account in the pIVW estimator. By default, plei=TRUE.

`sel.pval`:	
A vector containing the P values of the IV-exposure associations, which will be used for the IV selection. "sel.pval" should be provided when "delta" is not zero.

`delta`:	
The z-score threshold for IV selection. By default, delta=0 (i.e., no IV selection will be conducted).

`Boot.Fieller`:
If Boot.Fieller=TRUE, then the P value and the confidence interval of the causal effect based on the bootstrapping Fieller method will be calculated. By default, Boot.Fieller=TRUE.

`sig`:
The 100(1-sig)% confidence interval of the causal effect is calculated. By default, sig=0.05.

## Example 
```
library(mr.pivw)  # load the mr.pivw package
mr.pivw(data=example) # analyze the data with the pIVW method. 

$beta
[1] 0.5602797         # The estimated causal effect of the exposure on the outcome

$se
[1] 0.2913836         # The estimated standard error of estimated causal effect

$pval (Normal)
[1] 0.05450206        # The P value for testing whether the causal effect is zero, which is based on the normal approximation.

$CI (Normal)
     lower bound upper bound
[1,] -0.01082175    1.131381        # The 95% confidence interval of the causal effect based on the normal approximation.

$pval (Bootstrap Fieller)   
[1] 0.059                           # The P value for testing whether the causal effect is zero, which is based on the bootstrapping Fieller method.

$CI (Bootstrap Fieller)          
     lower bound upper bound       
[1,] -0.02538214     1.27153        # The 95% confidence interval of the causal effect based on the bootstrapping Fieller method.

$tau2
[1] 0.0004674883     # The variance of the horizontal pleiotropy
```

## Reference
Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. <doi:10.1111/biom.13732>

