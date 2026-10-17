# Reproduces Figure 5 (projection test power curves, all 3 error processes)
# and the Projection-based test section of results/README.md.
#
# Usage (from repo root): Rscript scripts/projection_power_curve.R
#
# Generates 1000 fresh runs per (process, sample size) combination -- same
# scale as the committed results. On a single machine this is slow; the
# committed results were generated on a cluster (see mc.cores below -- set
# it to the number of cores you actually have, never more than 2 on a
# personal laptop you don't own exclusively).

setwd("projection")
source("gen_data_proj_cp.R")
source("../cluster/run_parallel.R")

seeds <- c(bm = 418670, ns = 905327, ma = 261849)
d_map <- c(bm = 15, ns = 11, ma = 11)
process_labels <- c(bm = "Brownian motion", ns = "Non-separable", ma = "fMA(1)")
mc.cores <- min(40, parallel::detectCores())

ns <- c(150, 250, 500)
for (ty in c("bm", "ns", "ma")) {
    type <- ty
    ds <- c(d_map[[ty]])
    data_path <- paste(getwd(), "/data/projection/cp/", ty, sep = "")
    cat("=== generating:", ty, "(d=", ds, ") ===\n")
    generate_data_parallel(num_runs = 1000, mc.cores = mc.cores, seed = seeds[[ty]])
}

source("data_analysis_proj_cp.R")
png("../results/projection/power-curve-3process.png", width = 1500, height = 500)
par(mfrow = c(1, 3))
for (ty in c("bm", "ns", "ma")) {
    type <- ty
    analyze_all(Delta = 0.15, d = d_map[[ty]])
    title(main = paste0(process_labels[[ty]], " (d=", d_map[[ty]], ")"), cex.main = 2.0)
}
dev.off()
cat("wrote results/projection/power-curve-3process.png\n")

for (ty in c("bm", "ns", "ma")) {
    type <- ty
    d <- d_map[[ty]]
    cat("=== type:", ty, "d=", d, "===\n")
    for (n in c(150, 250, 500)) {
        directory <- paste(getwd(), "/data/projection/cp/", ty, "/", n, sep = "")
        set_directory(directory)
        read_data()
        cat("n=", n, ": ", sep = "")
        for (dm in c(0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)) {
            res <- sapply(1:get_num_of_runs(), function(r) test(0.15, r, dm, q_95, d))
            cat(sprintf("%.3f ", mean(res)))
        }
        cat("\n")
    }
}
