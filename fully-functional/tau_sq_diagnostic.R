# New functionality: Monte Carlo estimate of tau^2_{delta,theta0} (eq 3.7
# in the paper), the asymptotic variance in Theorem 3.1, needed for Figure
# 2 (the diagnostic plot of (Delta - ||delta||^2)/tau_{delta,theta0}).
#
# tau^2_{delta,theta0} is a quadratic form in delta (eq 3.7 integrates a
# product of two copies of delta against a fixed kernel), so for
# delta = sqrt(gamma)*shape, tau^2 scales linearly in gamma = ||delta||^2.
# We therefore only need ONE Monte Carlo estimate per process, at a
# reference ||delta||^2 (not one per candidate ||delta||^2 value), using
# the asymptotic characterization in Theorem 3.1:
#   tau^2_{delta,theta0} = lim_n n * Var(||D_n(1,theta0)||^2),
# with the TRUE theta0 (no change-point estimation needed -- (3.5) uses
# theta0, not the estimated theta_hat_n).
source("gen_data_cp.R")
suppressMessages(library(parallel))

estimate_tau_sq <- function(n, num_reps, ref_delta, cp_at, mc.cores, seed){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    gamma <- sqrt(get_gamma_norm(data_mean = ref_delta, S = S, type = type))

    run_one <- function(rep){
        if (type == "bm"){
            samp <- generate_brownian_motion_sample(n = n, S = S)
        } else if (type == "ns"){
            samp <- generate_nonseparable_sample(n = n, S = S)
        } else if (type == "ma"){
            samp <- generate_ma_sample(n = n, S = S)
        }
        shifted <- get_shifted_sample(samp = samp, cp_at = cp_at, factor = gamma)
        first_samp <- get_lambda_subsample(lambda = cp_at, samp = shifted, second_half = FALSE)
        second_samp <- get_lambda_subsample(lambda = cp_at, samp = shifted, second_half = TRUE)
        m1 <- get_mean(lambda = 1, samp = first_samp)
        m2 <- get_mean(lambda = 1, samp = second_samp)
        norm_spatial(m1 - m2)^2
    }

    D_vals <- unlist(parallel::mclapply(1:num_reps, run_one, mc.cores = mc.cores,
                                         mc.preschedule = FALSE, mc.set.seed = TRUE))
    list(tau_sq_per_delta = n * var(D_vals) / ref_delta, D_vals = D_vals)
}
