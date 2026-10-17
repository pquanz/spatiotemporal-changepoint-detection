# Reproduces Figure 1 (fully-functional test power curves, all 3 error
# processes) and the Fully-functional test section of results/README.md.
#
# Usage (from repo root): Rscript scripts/fully_functional_power_curve.R
#
# Generates 1000 fresh runs per (process, sample size) combination -- same
# scale as the committed results. On a single machine this is slow; the
# committed results were generated on a cluster (see mc.cores below -- set
# it to the number of cores you actually have, never more than 2 on a
# personal laptop you don't own exclusively).

setwd("fully-functional")
source("gen_data_cp.R")
source("../cluster/run_parallel.R")

seeds <- c(bm = 847213, ns = 592084, ma = 736951)
process_labels <- c(bm = "Brownian motion", ns = "Non-separable", ma = "fMA(1)")
mc.cores <- min(40, parallel::detectCores())

ns <- c(150, 250, 500)
for (ty in c("bm", "ns", "ma")) {
    type <- ty
    data_path <- paste(getwd(), "/data/fully functional/cp/", ty, sep = "")
    cat("=== generating:", ty, "===\n")
    generate_data_parallel(num_runs = 1000, mc.cores = mc.cores, seed = seeds[[ty]])
}

source("data_analysis_ff_cp.R")
png("../results/fully-functional/power-curve-3process.png", width = 1500, height = 500)
par(mfrow = c(1, 3))
for (ty in c("bm", "ns", "ma")) {
    type <- ty
    analyze_all(Delta = 0.15)
    title(main = process_labels[[ty]], cex.main = 2.0)
}
dev.off()
cat("wrote results/fully-functional/power-curve-3process.png\n")

for (ty in c("bm", "ns", "ma")) {
    type <- ty
    cat("=== type:", ty, "===\n")
    for (n in c(150, 250, 500)) {
        directory <- paste(getwd(), "/data/fully functional/cp/", ty, "/", n, sep = "")
        set_directory(directory)
        read_data()
        cat("n=", n, ": ", sep = "")
        for (dm in c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) {
            res <- sapply(1:get_num_of_runs(), function(r) test(0.15, r, dm, q_95))
            cat(sprintf("%.3f ", mean(res)))
        }
        cat("\n")
    }
}
