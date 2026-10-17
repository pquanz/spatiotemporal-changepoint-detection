# One-off calibration script: numerically estimates, for each error process
# used in the paper's simulation study (bm/ns/ma), the population covariance
# operator's eigenfunctions (from a large reference sample), then computes
# c_k = <shape, v_k> for the unscaled mean-shift shape(s,t) = s*cos(pi/2*t)
# and the cumulative explained-variance ratio e(d). These feed directly into
# get_gamma() in shared/gamma.R (c_k) and the paper's "95%-variance" d choice
# (e(d)), without needing a hand-derived closed-form eigendecomposition.
#
# Rationale: get_gamma_norm() (used by the fully-functional test) doesn't
# depend on the error process's eigenfunctions, since ||delta||^2 is a
# property of delta alone. get_gamma() (used by the projection test) does,
# since Sum_k <delta,v_k>^2 depends on which v_k. Rather than deriving each
# process's eigenfunctions by hand (real risk of a subtle error for the new
# non-separable process, which has no simple closed form due to its
# spatial-temporal coupling), we estimate them empirically from a large
# simulated sample -- the same numerical approach already used elsewhere in
# this repo (e.g. the Canada application's eigensystem), and very plausibly
# how the paper itself chose d=15/11/11 (that figure doesn't match a naive
# theoretical rank-1-in-space derivation for the "bm" process, but does match
# an empirically-estimated finite-sample covariance, which picks up small
# nonzero "noise" eigenvalues in directions that are theoretically zero).

setwd("..")
source("shared/gen_sample.R")
library(RSpectra)
library(fda)

S <- 4
n_ref <- 3000      # reference sample size for covariance estimation
max_d <- 80        # max projection dimension we might ever need

shape_fn <- function(s, t) s * cos(pi/2 * t)

covariance_fnc <- function(func1, func2){
    sample_size <- dim(func1$coefs)[2]
    c_func1 <- center.fd(func1)
    c_func2 <- center.fd(func2)
    final_coef_mat <- 0
    for (i in 1:sample_size){
        coef_mat <- c_func1$coefs[,i] %*% t(c_func2$coefs[,i])
        rownames(coef_mat) <- colnames(coef_mat)
        final_coef_mat <- final_coef_mat + coef_mat
    }
    final_coef_mat <- 1/sample_size * final_coef_mat
    coefnames <- list(dimnames(c_func1$coefs)[[1]], dimnames(c_func2$coefs)[[1]])
    bifd(final_coef_mat, c_func1$basis, c_func2$basis, coefnames)
}

evaluate_sample_in_s <- function(s, samp){
    n <- length(samp$coefs[1,])/S
    selected_cols <- ((1:n) - 1) * S + s
    subsample <- samp
    subsample$fdnames$reps <- samp$fdnames$reps[selected_cols]
    subsample$coefs <- samp$coefs[,selected_cols]
    subsample
}

inprod_spatial <- function(f, g){
    inprods <- inprod(f, g)
    s <- 0
    for (i in 1:S) s <- s + inprods[i,i]
    s
}

calibrate <- function(samp, max_d){
    n_basis <- samp$basis$nbasis

    final_mat <- NULL
    for (s1 in 1:S){
        row_block <- NULL
        for (s2 in 1:S){
            cf <- covariance_fnc(evaluate_sample_in_s(s1, samp), evaluate_sample_in_s(s2, samp))
            row_block <- cbind(row_block, cf$coefs)
        }
        final_mat <- rbind(final_mat, row_block)
    }

    all_basis <- fd(coef = diag(1, nrow = n_basis, ncol = n_basis), basisobj = samp$basis)
    gram_mat <- inprod(all_basis, all_basis)
    joint_mat <- final_mat %*% kronecker(diag(S), gram_mat)

    eig <- eigs(joint_mat, k = max_d, which = "LM")
    vals <- Re(eig$values)
    vecs <- Re(eig$vectors)

    # project shape(s,t) onto each eigenvector (same construction as
    # get_eigenvector() elsewhere: stack S blocks of n_basis coefficients,
    # normalize to unit spatial norm)
    t_grid <- seq(0, 1, length.out = 301)
    shape_block <- NULL
    for (s in 1:S){
        y <- shape_fn(s, t_grid)
        fdobj <- Data2fd(argvals = t_grid, y = matrix(y, ncol = 1), basisobj = samp$basis, lambda = 0)
        shape_block <- c(shape_block, fdobj$coefs[,1])
    }

    c_k <- numeric(max_d)
    for (k in 1:max_d){
        vec_coefs <- matrix(vecs[,k], nrow = n_basis)
        eigfd <- fd(coef = vec_coefs, basisobj = samp$basis)
        nrm <- 0
        for (s in 1:S) nrm <- nrm + inprod(fd(matrix(vec_coefs[,s], ncol=1), samp$basis),
                                             fd(matrix(vec_coefs[,s], ncol=1), samp$basis))[1,1]
        if (nrm > 0) vec_coefs <- vec_coefs / sqrt(nrm)
        eigfd <- fd(coef = vec_coefs, basisobj = samp$basis)

        shape_mat <- matrix(shape_block, nrow = n_basis)
        ip <- 0
        for (s in 1:S){
            ip <- ip + inprod(fd(matrix(shape_mat[,s], ncol=1), samp$basis),
                               fd(matrix(vec_coefs[,s], ncol=1), samp$basis))[1,1]
        }
        c_k[k] <- ip
    }

    list(vals = vals, c_k = c_k)
}

set.seed(392018)
samp_bm <- generate_brownian_motion_sample(n_ref, S)
res_bm <- calibrate(samp_bm, max_d)
cat("bm c_k[1:20]:", round(res_bm$c_k[1:20], 5), "\n")
cat("bm vals[1:20]:", round(res_bm$vals[1:20], 5), "\n")
cat("bm e(d) for d=5,10,15,20:",
    sapply(c(5,10,15,20), function(d) round(sum(res_bm$vals[1:d])/sum(res_bm$vals), 4)), "\n\n")

set.seed(583920)
samp_ns <- generate_nonseparable_sample(n_ref, S)
res_ns <- calibrate(samp_ns, max_d)
cat("ns c_k[1:20]:", round(res_ns$c_k[1:20], 5), "\n")
cat("ns vals[1:20]:", round(res_ns$vals[1:20], 5), "\n")
cat("ns e(d) for d=5,10,11,15,20:",
    sapply(c(5,10,11,15,20), function(d) round(sum(res_ns$vals[1:d])/sum(res_ns$vals), 4)), "\n\n")

set.seed(204857)
samp_ma <- generate_ma_sample(n_ref, S)
res_ma <- calibrate(samp_ma, max_d)
cat("ma c_k[1:20]:", round(res_ma$c_k[1:20], 5), "\n")
cat("ma vals[1:20]:", round(res_ma$vals[1:20], 5), "\n")
cat("ma e(d) for d=5,10,11,15,20:",
    sapply(c(5,10,11,15,20), function(d) round(sum(res_ma$vals[1:d])/sum(res_ma$vals), 4)), "\n\n")

saveRDS(list(bm = res_bm, ns = res_ns, ma = res_ma), "cluster/calibration_result.rds")
cat("saved cluster/calibration_result.rds\n")
