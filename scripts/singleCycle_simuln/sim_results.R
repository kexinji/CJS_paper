# Get the files names
files <- list.files(pattern="*.csv")
# First apply read.csv, then cbind
myfiles <- do.call(cbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))

res <- rowMeans(myfiles)

theta <- res[1:11]
betaHat <- res[12:13]
f1Hat <- res[14:41]
f2Hat <- res[42:69]

# Record simulation results into csv file
write.csv(res, "sim_results_100", row.names = F)	


# > theta
 # [1]  1.49972814  0.49992872  1.90432531  0.57562808  0.63663826  0.00106279  8.98117184
 # [8] 10.11057326 13.48967809  0.13819028  2.03681944
# > param <- c(1.5, 0.5, 1.2, 0.8, 0.6, 2, 3, 2, 5, 0.15, 2)
# > theta-param
 # [1] -2.718603e-04 -7.127954e-05  7.043253e-01 -2.243719e-01  3.663826e-02 -1.998937e+00
 # [7]  5.981172e+00  8.110573e+00  8.489678e+00 -1.180972e-02  3.681944e-02
 
# > betaHat
# [1] 0.04354703 0.07393497
# #------------------------------------------------------------------------------
# # Initializaiton of fixed effects, beta1
# beta11 <- 0.03
# beta21 <- 0.07


