# Reproduces Figure 7: the explained-variance ratio e(j) for j=1..100,
# j |-> sum(lambda_1..lambda_j) / sum(all computed eigenvalues), for the
# Canada data's estimated covariance operator. This is how d=51 (used by
# the projection test) was chosen: e(51) ~= 0.849.
#
# Usage (from repo root): Rscript canada-weather/scripts/eigenvalue_curve.R

setwd("canada-weather")
source("canada_application.R")

create_sample()
n <- length(years)
CP_est <<- get_cp_estimate(samp = sm)
egs <- get_eigensystem_by_nlam(n, sm)

e_j <- cumsum(Re(egs$vals)) / sum(Re(egs$vals))

dir.create("results", showWarnings = FALSE)
png("results/eigenvalue_curve.png", width = 700, height = 500)
plot(1:length(e_j), e_j, type = "l", xlab = "Number of eigenvalues j",
     ylab = "e(j)", main = "Explained-variance ratio e(j), Canada data",
     ylim = c(0, 1), cex.lab = 2.2, cex.axis = 1.8, cex.main = 1.8)
abline(h = e_j[51], lty = 2, col = "darkgray")
abline(v = 51, lty = 2, col = "darkgray")
text(51, e_j[51], sprintf("e(51) = %.4f", e_j[51]), pos = 4, cex = 1.7)
dev.off()
cat(sprintf("e(51) = %.4f [paper: ~0.849]\n", e_j[51]))
cat("wrote results/eigenvalue_curve.png\n")
