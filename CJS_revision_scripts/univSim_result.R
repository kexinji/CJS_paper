setwd("~/Downloads/CJS_revision/CJS_revision_scripts")

source("univSim_result_fn.R")

file.path.1 <- "../CJS_revision_results/univSim1/"
file.path.2 <- "../CJS_revision_results/univSim2/"

theta.1 <- c(1, rep(1, 3), 0.2, -0.44, 0.3, -0.2)
theta.2 <- c(0.75, rep(1, 3), 0.15, -1.60, 0.3, -0.1)

result.1 <- getResult(file.path.1, theta.1)
result.2 <- getResult(file.path.2, theta.2)