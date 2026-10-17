library(fda)            # functional data analysis library


# parameters
N = 300         # number of grid points (x-axis)


# Error process (5.2) in Dette & Quanz (2023): eta_i(s,t) = (s/S) * W_i(t),
# a single standard Brownian motion per replication i, scaled by s/S (reduces
# to the paper's s/4 at S=4, the only S value used in the paper).
generate_brownian_motion_sample <- function(n, S){

    # parameters:
    # n:    sample size
    # S:    number of elements in spatial set S
    #       (e.g. S = 4 corresponds to S = {1,2,3,4})

    data = matrix(rnorm(N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0

    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))

    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)
        bm_path <- 1/sqrt(N) * cumsum(data[,a_k(i)])
        values[2:(N+1),i] <- (arg_s/S) * bm_path
    }

    basis = create.fourier.basis(nbasis = 149)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)

    return(process)
}


# Error process (5.3): eta_i(s,t) = sum_{k=1}^D (N_ki/(2*pi*k)) *
# (sin(2*pi*k*t) + s*cos(2*pi*k*t)), a non-separable (spatial-temporal
# coupled) process. The paper truncates at infinity; it does not state a
# truncation for this process specifically (only for the fMA(1) process
# below), so we use D=40 terms for consistency with that stated choice.
generate_nonseparable_sample <- function(n, S, D = 40){

    N_coef <- matrix(rnorm(D*n), D, n)

    t_grid <- 0:N/N
    sin_mat <- sin(outer(t_grid, 1:D) * 2 * pi)   # (N+1) x D
    cos_mat <- cos(outer(t_grid, 1:D) * 2 * pi)   # (N+1) x D
    scale_k <- 1 / (2 * pi * (1:D))

    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    values <- matrix(0, nrow = N+1, ncol = n*S)

    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)
        rep_idx <- a_k(i)
        amp <- N_coef[, rep_idx] * scale_k
        values[, i] <- sin_mat %*% amp + arg_s * (cos_mat %*% amp)
    }

    basis = create.fourier.basis(nbasis = 2*D + 1)
    process = Data2fd(argvals = t_grid, y = values, basisobj = basis, lambda = 0)

    return(process)
}


# Error process (5.4)/(5.5): eps_i(s,t) = sum_{j=1}^D N_ji*sqrt(lam_j)*v_j(s,t),
# lam_j = 1/(2*pi*j^2), v_j(s,t) = (sqrt(2)/4)*s*sin(2*pi*j*t); fMA(1):
# eta_i = eps_i + kappa*eps_{i-1}. D=40 terms, as stated in the paper.
generate_ma_sample <- function(n, S, kappa = 0.7, D = 40){

    zlag0 <- matrix(rnorm(D*n), D, n)
    zlag1 <- cbind(0, zlag0[, -n, drop = FALSE])
    coef  <- kappa * zlag1 + zlag0          # D x n, combined N-coefficients for eta_i

    lam_j <- 1 / (2 * pi * (1:D)^2)
    sqrt_lam <- sqrt(lam_j)

    t_grid <- 0:N/N
    sin_mat <- sin(outer(t_grid, 1:D) * 2 * pi)   # (N+1) x D

    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    values <- matrix(0, nrow = N+1, ncol = n*S)

    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)
        rep_idx <- a_k(i)
        amp <- coef[, rep_idx] * sqrt_lam * (sqrt(2)/4) * arg_s
        values[, i] <- sin_mat %*% amp
    }

    basis <- create.fourier.basis(nbasis = 2*D + 1)
    process <- Data2fd(argvals = t_grid, y = values, basisobj = basis, lambda = 0)

    return(process)
}
