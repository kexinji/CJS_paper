source("bivLiu_sc_result_fn.R")

getTable(result.path = "../CJS_revision_results/", 
         result.name = "bivLiu100_sc_cor.txt", 
         s = 16, 
         b = 6)

# calculate corr in example
num <- -0.2996  + 0.3176*exp((-0.9222 + 0.1570*5 -0.0704*5^2-2.2810 + 1.4265*6 -0.5533*6^2)/2)
den <- sqrt(0.6676 + 0.6921 + exp(-0.9222 + 0.1570*5 -0.0704*5^2)) * sqrt(0.6690  + 0.6229 + exp(-2.2810 + 1.4265*6 -0.5533*6^2))

num/den