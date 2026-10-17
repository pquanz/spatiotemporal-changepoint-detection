# Reproduces Figure 3: histogram of the change-point estimator, n=250,
# Delta=0.15, 1000 runs, one panel per error process.
#
# Usage (from repo root): Rscript scripts/changepoint_histogram.R
#
# Much cheaper than the power-curve scripts (one fixed shift magnitude
# instead of a d_mean sweep). mc.cores: see note in
# fully_functional_power_curve.R.

setwd("fully-functional")
source("gen_cp_histogram.R")

seeds <- c(bm = 503917, ns = 672140, ma = 388256)
process_labels <- c(bm = "Brownian motion", ns = "Non-separable", ma = "fMA(1)")
mc.cores <- min(40, parallel::detectCores())

out_dir <- "data/cp-histogram"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()
for (ty in c("bm", "ns", "ma")) {
    type <- ty
    cat("=== cp histogram:", ty, "===\n")
    res <- generate_cp_histogram_data(num_runs = 1000, n = 250, seed = seeds[[ty]],
                                       mc.cores = mc.cores, target_delta = 0.15)
    results[[ty]] <- res
    write.csv(data.frame(cp_est = res), file.path(out_dir, paste0(ty, ".csv")), row.names = FALSE)
    cat(ty, ": mean=", mean(res), " sd=", sd(res), "\n")
}

png("../results/fully-functional/changepoint-histogram.png", width = 1500, height = 500)
par(mfrow = c(1, 3), mar = c(5.5, 6, 4, 2) + 0.1)
for (ty in c("bm", "ns", "ma")) {
    hist(results[[ty]], breaks = 30, xlim = c(0, 1), main = process_labels[[ty]],
         xlab = expression(hat(vartheta)[n]), col = "grey80",
         cex.lab = 2.2, cex.axis = 1.8, cex.main = 2.0)
    abline(v = 0.6, col = "darkgray", lwd = 2, lty = 2)
}
dev.off()
cat("wrote results/fully-functional/changepoint-histogram.png\n")
