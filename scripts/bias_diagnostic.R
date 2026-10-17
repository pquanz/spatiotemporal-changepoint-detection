# Reproduces Figure 2: (Delta - ||delta||^2) / tau_{delta,theta0} as a
# function of ||delta||^2, one line per error process -- the quantity that
# explains the power differences seen in Figure 1 (a smaller tau relative
# to the Delta - ||delta||^2 gap means higher power).
#
# Usage (from repo root): Rscript scripts/bias_diagnostic.R
#
# Cheap: tau^2 is a quadratic form in delta, so only ONE Monte Carlo
# estimate per process is needed (at a reference ||delta||^2), not one per
# grid point -- see fully-functional/tau_sq_diagnostic.R.
#
# The Monte Carlo estimate itself is cached in data/tau-sq-cache.csv, since
# it's the only genuinely expensive part of this script (the rest is just
# plotting). Delete that file to force a recompute (e.g. after changing the
# seeds or num_reps below).

setwd("fully-functional")
source("tau_sq_diagnostic.R")

seeds <- c(bm = 182340, ns = 182341, ma = 182342)
mc.cores <- min(40, parallel::detectCores())

delta_grid <- seq(0.01, 0.3, by = 0.01)
Delta <- 0.15

cache_path <- "data/tau-sq-cache.csv"
if (file.exists(cache_path)) {
    cached <- read.csv(cache_path)
    tau_sq_per_delta <- setNames(cached$tau_sq_per_delta, cached$process)
    cat("using cached tau_sq_per_delta from", cache_path, "\n")
} else {
    tau_sq_per_delta <- c()
    for (ty in c("bm", "ns", "ma")) {
        type <- ty
        res <- estimate_tau_sq(n = 500, num_reps = 2000, ref_delta = 0.15, cp_at = 0.6,
                                mc.cores = mc.cores, seed = seeds[[ty]])
        tau_sq_per_delta[ty] <- res$tau_sq_per_delta
        cat(ty, ": tau_sq_per_delta =", res$tau_sq_per_delta, "\n")
    }
    dir.create("data", showWarnings = FALSE)
    write.csv(data.frame(process = names(tau_sq_per_delta), tau_sq_per_delta = tau_sq_per_delta),
              cache_path, row.names = FALSE)
}

curves <- list()
for (ty in c("bm", "ns", "ma")) {
    tau <- sqrt(delta_grid * tau_sq_per_delta[[ty]])
    curves[[ty]] <- (Delta - delta_grid) / tau
}

png("../results/fully-functional/bias-diagnostic.png", width = 800, height = 550)
par(mar = c(5.5, 8, 4.5, 2) + 0.1)
plot(delta_grid, curves[["bm"]], type = "l", lty = 1, ylim = range(unlist(curves)),
     xlab = expression(paste("||", delta, "||"^2)),
     ylab = expression((Delta - "||" * delta * "||"^2) / tau[list(delta, vartheta[0])]),
     main = "Size of the bias term (Figure 2)",
     cex.lab = 1.9, cex.axis = 1.8, cex.main = 2.0)
lines(delta_grid, curves[["ns"]], lty = 2)
lines(delta_grid, curves[["ma"]], lty = 3)
abline(v = Delta, lwd = 0.5)
legend("topright", c("bm", "ns", "ma"), lty = 1:3, cex = 1.7)
dev.off()
cat("wrote results/fully-functional/bias-diagnostic.png\n")
