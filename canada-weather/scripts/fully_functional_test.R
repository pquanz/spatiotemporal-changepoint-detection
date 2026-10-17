# Reproduces Table 1 (Canada fully-functional test results) in
# results/README.md.
#
# Usage (from repo root): Rscript canada-weather/scripts/fully_functional_test.R

setwd("canada-weather")
source("canada_application.R")

create_sample()
res <- relevant_change_test(samp = sm, n = length(years))

cat(sprintf("D_hat = %.5f, V_hat = %.6f\n", res$D_hat, res$V_hat))
cat("max_delta (still-rejected threshold) at 90/95/99%:\n")
print(res$max_delta)
cat(sprintf("one-sided CI for ||delta||: [%.3f, %.3f]\n", res$ci_one_sided[1], res$ci_one_sided[2]))
cat(sprintf("two-sided CI for ||delta||: [%.3f, %.3f]\n", res$ci_two_sided[1], res$ci_two_sided[2]))
