# Reproduces the Canada change-point estimate reported in results/README.md
# and the front-page README: theta_hat_n = 0.7079646 [paper: 0.708].
#
# Usage (from repo root): Rscript canada-weather/scripts/changepoint_estimate.R

setwd("canada-weather")
source("canada_application.R")

create_sample()
cp_est <- get_cp_estimate(samp = sm)

# model (2.1) in the paper: mu_1=...=mu_{floor(n*theta)}, mu_{floor(n*theta)+1}=...=mu_n
# so floor(n*theta) is the last pre-change observation, not a simple offset
# from years[1] -- years[] has 4 excluded (missing-data) years, so position
# in the trimmed sample does not line up with calendar year by arithmetic.
n <- length(years)
k <- floor(n * cp_est)

cat(sprintf("theta_hat_n = %.7f [paper: 0.708]\n", cp_est))
cat(sprintf("last pre-change year: %d, first post-change year: %d\n", years[k], years[k + 1]))
cat("(the paper states '~1969'; our own re-derivation from the trimmed years\n")
cat(" vector, using the paper's own model (2.1) convention, gives 1971/1972 --\n")
cat(" we report our own computation rather than force an exact match)\n")
