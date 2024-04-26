setMethod("show",
          "PIVW",
          function(object){
            decimals = function(number, places){format(round(number, places), nsmall = places)}

            Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
            Statistic <- c("Method", "Estimate", "Std Error", Interval_type, "", "p-value","Condition")

            blk = rep("", length(object@CILower)-1)
            Value <- c("pIVW", blk,
                       decimals(object@Estimate, 3), blk,
                       decimals(object@StdError,3), blk,
                       paste(decimals(object@CILower, 3), ",", sep = ""),
                       decimals(object@CIUpper, 3),
                       signif(object@Pvalue, 3), blk,
                       decimals(object@Condition, 3), blk)

            output.table <- data.frame(matrix(Value, nrow = length(object@CILower), byrow=FALSE))
            colnames(output.table) <- Statistic

            cat("\nPenalized inverse-variance weighted method\n\n")
            cat("Account for over-dispersion:", object@Over.dispersion, "\n")
            if(object@Boot.Fieller){
              cat("CI and P-value:", "Bootstrapping Fieller", "\n")
              cat("Number of bootstrap samples:", object@N.boot, "\n")
            }else{
              cat("CI and P-value:", "Normal approximation", "\n")
            }
            cat("Penalty parameter (lambda):", object@Lambda, "\n")
            cat("IV selection threshold (delta):", round(object@Delta,3), "\n")
            cat("Number of variants:", object@SNPs, "\n")


            cat("------------------------------------------------------------------\n")
            print(output.table, quote = F, row.names = FALSE, justify= "left")
            cat("------------------------------------------------------------------\n")
          }
)
