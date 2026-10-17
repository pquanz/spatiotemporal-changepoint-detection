# New functionality (not part of the original get_export_matrix() pipeline):
# generates ONLY the change-point estimate per run, at a fixed mean-shift
# ||delta||^2 = target_delta. This is what Figure 3 of the paper needs
# (histogram of the change-point estimator, n=250, Delta=0.15, 1000 runs,
# one panel per error process) -- much cheaper than the full d_mean sweep
# in get_export_matrix(), since only one shift magnitude is needed.
source("gen_data_cp.R")
suppressMessages(library(parallel))
suppressMessages(library(data.table))

get_cp_estimate_at_delta <- function(samp, target_delta = 0.15){
    gamma <- sqrt(get_gamma_norm(data_mean = target_delta, S = S, type = type))
    shifted_samp <- get_shifted_sample(samp = samp, cp_at = CP_at, factor = gamma)
    get_cp_estimate(samp = shifted_samp)
}

generate_cp_histogram_data <- function(num_runs, n, seed, mc.cores, target_delta = 0.15){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    run_one <- function(run){
        if (type == "bm"){
            samp <- generate_brownian_motion_sample(n = n, S = S)
        } else if (type == "ns"){
            samp <- generate_nonseparable_sample(n = n, S = S)
        } else if (type == "ma"){
            samp <- generate_ma_sample(n = n, S = S)
        }
        get_cp_estimate_at_delta(samp, target_delta = target_delta)
    }

    results <- mclapply(1:num_runs, run_one, mc.cores = mc.cores,
                         mc.preschedule = FALSE, mc.set.seed = TRUE)
    errs <- which(vapply(results, function(x) inherits(x, "try-error"), logical(1)))
    if (length(errs) > 0){
        stop(paste("run(s) failed:", paste(errs, collapse = ", "),
                   "\nfirst error:", attr(results[[errs[1]]], "condition")$message))
    }
    unlist(results)
}
