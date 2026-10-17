# This file solely purposes the calculation of specific values
# for gamma in the mean function


hyp = "proj"    # the hypothesis that should be tested


# Projection-test gamma calibration (get_gamma): Sum_{k=1}^d <delta,v_k>^2
# depends on the error process's covariance eigenfunctions v_k, since those
# are process-specific. c_k below are <shape,v_k> for the unscaled mean-shift
# shape(s,t) = s*cos(pi/2*t) (so delta = sqrt(gamma)*shape), estimated
# numerically from each process's population covariance operator (large
# reference sample), rather than hand-derived in closed form -- see
# cluster/calibrate_gamma.R for the full derivation and provenance. This
# avoids a real risk of a subtle hand-derivation error for the non-separable
# process (5.3), which has no simple closed-form eigendecomposition due to
# its spatial-temporal coupling.
#
# get_gamma_norm (fully-functional test) does NOT need this: ||delta||^2 is
# a property of delta alone, independent of the error process, so it's the
# same formula for every type.

.c_k_bm <- c(2.552495, 2.438923, -0.828451, -0.805304, 0.491890, -0.515516,
             -0.320506, 0.345132, -0.300090, 0.314742, -0.230357, 0.272911,
             0.136609, -0.182941, -0.160241, 0.171894, -0.206089, -0.117478,
             0.080644, -0.103690)

.c_k_ns <- c(-0.101895, -0.123858, -0.096063, 0.080866, -0.058125, -0.056002,
             -0.050670, -0.046969, 0.033583, -0.031124, -0.026901, 0.034448,
             -0.008969, -0.011964, 0.023032)

.c_k_ma <- c(1.323133, -0.624128, -0.402587, -0.311098, 0.205433, 0.188408,
             -0.184286, 0.182318, 0.120596, 0.138121, -0.108541, -0.121809,
             -0.099343, -0.100950, 0.035634)


get_gamma <- function(data_mean, d, S, type = "bm"){

    # for fast return
    if (data_mean == 0){
        return (0)
    }

    c_k <- switch(type, "bm" = .c_k_bm, "ns" = .c_k_ns, "ma" = .c_k_ma)
    if (d > length(c_k)){
        stop(paste("get_gamma: d =", d, "exceeds the calibrated range for type", type,
                    "-- rerun cluster/calibrate_gamma.R with a larger max_d."))
    }

    return(data_mean / sum(c_k[1:d]^2))
}



get_gamma_norm <- function(data_mean, S, type = "bm"){

    if (data_mean == 0){
        return (0)
    }

    # ||delta||^2 = gamma * sum(s^2) * integral(cos(pi/2*t)^2) dt
    #             = gamma * sum(s^2) / 2, independent of the error process
    sum_sq_s <- Reduce("+", c(1:S)^2)
    return(2 * data_mean / sum_sq_s)
}
