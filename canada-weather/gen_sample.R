library(fda)            # functional data analysis library


# parameters
N = 300         # number of grid points (x-axis)




gen_s <- function(n, S){
    kappa = 0
    # number of summands in Karhunen-Loeve expansion
    D = 5
    
    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[, i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0   # ranodm coefs
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    # norm_factor <- (function(k) sqrt(S/2 - (S * (S+1)/(pi * (1-2*k)^2) ) * 4 * cos(pi * k)^2 + 10 * sin( (k-1/2) * pi)) ) # 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    # norm_factor <- c((13*pi - 80)/pi, (99*pi - 80)/(9 * pi), (65*pi - 16)/(5*pi), (539*pi - 80)/(49*pi), (1053*pi - 80)/(81*pi) )
    # v_k <- (function(k, t, s) sqrt(15) * sin( (k -0.5) * pi * t) - (t-s) * sin((k-0.5)*pi))
    v_k <- (function(k, t, s) 1/sqrt(17) * (sin( (k ) *2 * pi * t) + s * cos((k) * 2 * pi * t)))
    # v_k <- (function(k, t, s) 1/s * sin(k * pi * t))
    lam_k <- (function(k) 1 / ( (k) * 2 * pi))
    
    # initialize matrix, which will store values of the process
    # values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    final_coefs <- matrix(replicate((n*S * (2*D + 1)), 0), nrow = 2*D + 1, ncol = n*S)
    values <- matrix(replicate(n*S * 301, 0), ncol=n*S)

    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    # loop through all values to calculate and assign them
    for (i in 1:(n*S)){
        
        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)
        
        for (j in 0:N) {
            # evaluation in v_k, in order to get the pre factor of v_k 
            # chosen such that sin(...) = 1
            # print(paste(coef[j, a_k(i)], sqrt(lam_k(j)), v_k(1, 0.25, arg_s)))
            
            to_add <- 0
            for(k in 1:D){
                to_add <- to_add + (coef[k, a_k(i)] * sqrt(lam_k(k)) * v_k(k, j/N, arg_s)) 
            }
            values[j, i] <- to_add
            # to_add <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0.25, arg_s)
            # final_coefs[(2*j), i] <- to_add
            
            # print(to_add)
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    basis = create.fourier.basis(nbasis = 149)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)
    return(process)

}








gen_non_sep <- function(n, S, kappa=0){
    # creates a sample of moving average processes
    
    # parameters:
    # n:        sample size
    # S:        number of elements in spatial set S
    #           (e.g. S = 4 corresponds to S = {1,2,3,4})
    # kappa:    dependence parameter in MA(1)-process 
    
    
    # number of summands in Karhunen-Loeve expansion
    D = 20
    
    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[, i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0   # ranodm coefs
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    norm_factor <- 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    v_k <- (function(k, t, s) 1/sqrt(17) * (sin( (k ) *2 * pi * t) + s * cos((k) * 2 * pi * t)))
    # v_k <- (function(k, t, s) 1/s * sin(k * pi * t))
    lam_k <- (function(k) 1 / ( (k) * 2 * pi)^2)
    
    # initialize matrix, which will store values of the process
    # values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    final_coefs <- matrix(replicate((n*S * (2*D + 1)), 0), nrow = 2*D + 1, ncol = n*S)
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    # loop through all values to calculate and assign them
    for (i in 1:(n*S)){
        
        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)
        
        for (j in 1:D) {
            # evaluation in v_k, in order to get the pre factor of v_k 
            # chosen such that sin(...) = 1
            # print(paste(coef[j, a_k(i)], sqrt(lam_k(j)), v_k(1, 0.25, arg_s)))
            to_add_sine <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0.25, arg_s)
            final_coefs[(2*j), i] <- to_add_sine

            to_add_cosine <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0, arg_s)
            final_coefs[(2*j + 1), i] <- to_add_cosine
            # print(to_add)
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    process <- fd(coef = final_coefs, basisobj = create.fourier.basis(nbasis = 2 * D + 1))
    
    return(process)
}


gen_non_sep_old <- function(n, S, kappa=0){
    # creates a sample of moving average processes
    
    # parameters:
    # n:        sample size
    # S:        number of elements in spatial set S
    #           (e.g. S = 4 corresponds to S = {1,2,3,4})
    # kappa:    dependence parameter in MA(1)-process 
    
    
    # number of summands in Karhunen-Loeve expansion
    D = 20
    
    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[, i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0   # ranodm coefs
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    norm_factor <- 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    v_k <- (function(k, t, s) 1/sqrt(17) * sin( (k ) *2 * pi * t) + s * cos((k) * 2 * pi * t))
    # v_k <- (function(k, t, s) 1/s * sin(k * pi * t))
    lam_k <- (function(k) 1 / ( (k) * 2 * pi))
    
    # initialize matrix, which will store values of the process
    # values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    final_coefs <- matrix(replicate((n*S * (2*D + 1)), 0), nrow = 2*D + 1, ncol = n*S)
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    # loop through all values to calculate and assign them
    for (i in 1:(n*S)){
        
        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)
        
        for (j in 1:D) {
            # evaluation in v_k, in order to get the pre factor of v_k 
            # chosen such that sin(...) = 1
            # print(paste(coef[j, a_k(i)], sqrt(lam_k(j)), v_k(1, 0.25, arg_s)))
            to_add_sine <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0.25, arg_s)
            final_coefs[(2*j), i] <- to_add_sine

            to_add_cosine <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0, arg_s)
            final_coefs[(2*j + 1), i] <- to_add_cosine
            # print(to_add)
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    process <- fd(coef = final_coefs, basisobj = create.fourier.basis(nbasis = 2 * D + 1))
    
    return(process)
}











gen_new_bm <- function(n, S){
    # parameters:
    # n:        sample size
    # S:        number of elements in spatial set S
    #           (e.g. S = 4 corresponds to S = {1,2,3,4})
    # kappa:    dependence parameter in MA(1)-process 
    

    kappa = 0
    
    # number of summands in Karhunen-Loeve expansion
    D = 20
    
    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[, i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0   # ranodm coefs
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    
    # norm_factor <- 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    norm_factor <- sqrt(15)
    v_k <- (function(k, t, s) norm_factor * s * sin((k - 0.5) * pi * t) - (t-s))
    # v_k <- (function(k, t, s) 1/s * sin(k * pi * t))
    lam_k <- (function(k) 1 / ((k - 0.5) * pi))
    
    # initialize matrix, which will store values of the process
    # values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    final_coefs <- matrix(replicate((n*S * (2*D + 1)), 0), nrow = 2*D + 1, ncol = n*S)
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    # loop through all values to calculate and assign them
    for (i in 1:(n*S)){
        
        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)
        
        for (j in 1:D) {
            # evaluation in v_k, in order to get the pre factor of v_k 
            # chosen such that sin(...) = 1
            # print(paste(coef[j, a_k(i)], sqrt(lam_k(j)), v_k(1, 0.25, arg_s)))
            to_add <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0.25, arg_s)
            final_coefs[(2*j), i] <- to_add
            # print(to_add)
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    process <- fd(coef = final_coefs, basisobj = create.fourier.basis(nbasis = 2 * D + 1))
    
    return(process)
}


generate_brownian_motion_sample1 <- function(n, S){
    data = matrix(rnorm(N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))

    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)
        values[2:(N+1),i] <- arg_s/4 * 1/sqrt(N) * cumsum(data[,a_k(i)])
    }
    
    # create bspline basis and fd-object
    # basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
    basis = create.fourier.basis(nbasis = 149)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)

    return(process)
}

generate_brownian_bridge_sample1 <- function(n, S){
    
    # parameters:
    # n:    sample size
    # S:    number of elements in spatial set S
    #       (e.g. S = 4 corresponds to S = {1,2,3,4})
    
    # construction of BB data points
    data = matrix(rnorm(N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)

        cumulative = 1/sqrt(N) * cumsum(data[,a_k(i)])
        values[2:(N+1),i] =  cumulative - ((1:N/N) - sqrt(arg_s)) * cumulative[N]
    }

    # create bspline basis and fd-object
    basis = create.fourier.basis(nbasis = 149)
    # basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)
    
    return(process)
}


gen_bb2 <- function(n, S){
    
    # parameters:
    # n:    sample size
    # S:    number of elements in spatial set S
    #       (e.g. S = 4 corresponds to S = {1,2,3,4})
    
    # construction of BB data points
    data = matrix(rnorm(2*N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    for(i in 1:(n*S)) {
        arg_s <- (if (i %% S == 0) S else i %% S)
        cumulative = 1/sqrt(2*N) * cumsum(data[,a_k(i)])
        values[2:(N+1),i] =  cumulative[2:(N+1)] - arg_s * (1:N/N) * cumulative[2*N]
    }

    # create bspline basis and fd-object
    basis = create.fourier.basis(nbasis = 149)
    # basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)
    
    return(process)
}





generate_brownian_motion_sample <- function(n, S){

    # parameters:
    # n:    sample size
    # S:    number of elements in spatial set S
    #       (e.g. S = 4 corresponds to S = {1,2,3,4})
    
    # construction of BM data points
    data = matrix(rnorm(N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))

    for(i in 1:(n*S)) {
        values[2:(N+1),i] <- 1/sqrt(N) * cumsum(data[,a_k(i)])
    }
    
    # create bspline basis and fd-object
    # basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
    basis = create.fourier.basis(nbasis = 149)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)

    return(process)
}

generate_brownian_bridge_sample <- function(n, S){
    
    # parameters:
    # n:    sample size
    # S:    number of elements in spatial set S
    #       (e.g. S = 4 corresponds to S = {1,2,3,4})
    
    # construction of BB data points
    data = matrix(rnorm(N*n), ncol = n)
    values <- matrix(nrow = N+1, ncol = n*S)
    values[1,] <- 0
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    for(i in 1:(n*S)) {
        cumulative = 1/sqrt(N) * cumsum(data[,a_k(i)])
        values[2:(N+1),i] =  cumulative - (1:N/N) * cumulative[N]
    }

    # create bspline basis and fd-object
    basis = create.fourier.basis(nbasis = 149)
    # basis = create.bspline.basis(norder = 2, breaks = 0:N/N)
    process = Data2fd(argvals = 0:N/N, y = values, basisobj = basis, lambda = 0)
    
    return(process)
}









generate_ma_sample <- function(n, S, kappa=0.7){
    # creates a sample of moving average processes
    
    # parameters:
    # n:        sample size
    # S:        number of elements in spatial set S
    #           (e.g. S = 4 corresponds to S = {1,2,3,4})
    # kappa:    dependence parameter in MA(1)-process 
    
    
    # number of summands in Karhunen-Loeve expansion
    D = 20
    
    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[, i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0   # ranodm coefs
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    norm_factor <- 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    v_k <- (function(k, t, s) norm_factor * s * sin(2 * k * pi * t))
    # v_k <- (function(k, t, s) 1/s * sin(k * pi * t))
    lam_k <- (function(k) 1 / (2 * k^2 * pi))
    
    # initialize matrix, which will store values of the process
    # values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    final_coefs <- matrix(replicate((n*S * (2*D + 1)), 0), nrow = 2*D + 1, ncol = n*S)
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    
    # loop through all values to calculate and assign them
    for (i in 1:(n*S)){
        
        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)
        
        for (j in 1:D) {
            # evaluation in v_k, in order to get the pre factor of v_k 
            # chosen such that sin(...) = 1
            # print(paste(coef[j, a_k(i)], sqrt(lam_k(j)), v_k(1, 0.25, arg_s)))
            to_add <- coef[j, a_k(i)] * sqrt(lam_k(j)) * v_k(1, 0.25, arg_s)
            final_coefs[(2*j), i] <- to_add
            # print(to_add)
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    process <- fd(coef = final_coefs, basisobj = create.fourier.basis(nbasis = 2 * D + 1))
    
    return(process)
}














generate_ma_sample_old <- function(n, S, kappa=0.7){
    # creates a sample of moving average processes
    
    # parameters:
    # n:        sample size
    # S:        number of elements in spatial set S
    #           (e.g. S = 4 corresponds to S = {1,2,3,4})
    # kappa:    dependence parameter in MA(1)-process 
    

    # number of summands in Karhunen-Loeve expansion
    D = 250# 100# 30

    # initialize coefficient and lag matrices as zero matrix
    coef = zlag0 = zlag1 = matrix(0, D, n)
    
    # fill lag matrices with realizations of st. normal random variable
    for (i in 1:n){ zlag0[,i] = rnorm(D, 0, 1) }
    for (i in 2:n){ zlag1[,i] = zlag0[,i-1] }
    
    # calculate coefficients of MA(1)-process (X_t = kappa * Z_{t-1} + Z_t)
    coef = kappa * zlag1 + zlag0
    
    # define eigenfunctions and -values to be used in Karhunen-Loeve expansion
    # v_k <- (function(k, t, s) sqrt(2) * s / S * sin(2 * k * pi * t))
    norm_factor <- 2 * sqrt(3)/sqrt(S * (S+1) * (2*S + 1))
    v_k <- (function(k, t, s) norm_factor * s * cos(2 * k * pi * t))
    lam_k <- (function(k) 1 / (2 * k * pi))
    # lam_k <- (function(k) 1/(2*k^2 * pi))
    
    # initialize matrix, which will store values of the process
    values <- matrix(replicate((n*S), 0), nrow = N+1, ncol = n*S)
    
    # function to help assign values in s
    # 1, ..., 1 (S times), 2, ..., 2 (S times), 3, ..., 3, ...
    a_k <- (function(k) (k - k %% S)/S + (k %% S != 0))
    

    # loop through all values to calculate and assign them
    for (i in 1:length(values[1,])){

        # will take the values 1,2,3,4,1,2,3,4,1,... which will be plugged into v_k
        arg_s <- (if (i %% S == 0) S else i %% S)

        # loop through 1,...,D to sum up all summands in Karhunen-Loeve exp.
        for (k in 1:D){
            values[,i] <- values[,i] + (coef[k,a_k(i)] * lam_k(k) * v_k(k, 0:N/N, arg_s))
        }
    }
    
    # lastly create fd-object with fourier basis
    # twice as many basis elements, since every second is cosine
    # and we want the first 30 sines.
    process <- Data2fd(argvals = 0:N/N, y = values, 
                       basisobj = create.fourier.basis(nbasis = 2 * D + 1))
    
    return(process)
}
