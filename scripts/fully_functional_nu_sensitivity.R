# Reproduces Figure 4: sensitivity of the fully-functional test to the
# measure nu in the self-normalizing statistic, fMA(1) process only, three
# choices nu_4/nu_9/nu_19 (num_eigsys = 5/10/20 in this code's notation --
# support_lam has num_eigsys-1 points, matching nu_k with k=num_eigsys-1).
#
# Usage (from repo root): Rscript scripts/fully_functional_nu_sensitivity.R
#
# mc.cores: see note in fully_functional_power_curve.R.

setwd("fully-functional")
source("gen_data_cp.R")
source("../cluster/run_parallel.R")

nu_settings <- list(nu4 = list(num_eigsys = 5, seed = 714823),
                     nu9 = list(num_eigsys = 10, seed = 619472),
                     nu19 = list(num_eigsys = 20, seed = 736951))
nu_labels <- list(nu4 = expression(nu[4]), nu9 = expression(nu[9]), nu19 = expression(nu[19]))
mc.cores <- min(40, parallel::detectCores())

type <- "ma"
ns <- c(150, 250, 500)
for (nm in names(nu_settings)) {
    num_eigsys <- nu_settings[[nm]]$num_eigsys
    data_path <- paste(getwd(), "/data/nu-sensitivity/ff/", nm, sep = "")
    cat("=== generating:", nm, "(num_eigsys=", num_eigsys, ") ===\n")
    generate_data_parallel(num_runs = 1000, mc.cores = mc.cores,
                            seed = nu_settings[[nm]]$seed)
}

source("data_analysis_ff_cp.R")
png("../results/fully-functional/nu-sensitivity.png", width = 1500, height = 500)
par(mfrow = c(1, 3))
for (nm in names(nu_settings)) {
    path_prefix <- paste(getwd(), "/data/nu-sensitivity/ff/", nm, sep = "")
    analyze_all_directory(path_prefix, Delta = 0.15)
    title(main = nu_labels[[nm]], cex.main = 2.0)
}
dev.off()
cat("wrote results/fully-functional/nu-sensitivity.png\n")
