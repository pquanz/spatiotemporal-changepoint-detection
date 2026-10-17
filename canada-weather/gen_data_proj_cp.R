source("gen_sample.R")
source("gamma.R")
library(R.utils)        # ?? 
library(RSpectra)       # eigenvalue computation
library(data.table)     # for exporting data faster



# parameters

S = 4                         # number of spatial components
ds = c(25) # c(5, 10, 15, 20, 25, 30) # number of projections
num_eigsys = 20#20               # number of eigensystems
ns = c(150, 250, 500) # c(150, 250, 500)   # number of sample sizes

# range_data_mean <- c(0, 0.02, 0.03, 0.04, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
#                      0.13, 0.14, 0.15, 0.16, 0.18, 0.2, 0.225, 0.25, 0.5, 1, 2)

range_data_mean <- c(0, 0.025, 0.05, 0.075, 0.09, 0.095, 0.1, 0.105, 0.11, 0.12, 0.13,
                     0.14, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
# range_data_mean <- c(0, 0.02, 0.03, 0.1, 0.2, 0.3, 0.5, 1, 2)
# range_data_mean <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


CP_at <- 0.6
CP_est <- -1

test = "cp"             # test for cp
type = "nonsep"             # type of data (bm, bb, ma)

data_path = paste(getwd(), "/data/projection/", test, "/", type, sep="")

# params end



# mean function

t <- 0:N/N
mu <- (function(s, t) s * cos(0.5 * pi * t))

mean_fnc <- c()
for(s in 1:S){
    mean_fnc <- cbind(mean_fnc, mu(s,t))
}

basis = create.bspline.basis(norder = 2, breaks = t)
fd_mean <- Data2fd(t, mean_fnc, basisobj = basis)

# mean function end



### TECHNICAL FUNCTIONS ###

get_lambda_subsample <- function(lambda = 1, samp, second_half = F){
    # returns a subsample for specific value of lambda

    # parameters:
    # lambda:       the proportion of the sample that will be returned
    # samp:         the sample of which a subsample will be returned
    # second_half:  if true, returns the latter proportion of the sample
    
    n <- length(samp$coefs[1,])/S
    n_lam = floor(n * lambda)
    
    # if n_lam is zero return empty sample
    if(n_lam == 0){
        if(second_half){
            return(samp)
        }else{
            return(fd(coef = replicate(samp$basis$nbasis, 0), basisobj = samp$basis))   
        }
    }else if(n_lam == n){
        if(second_half){
            return(fd(coef = replicate(samp$basis$nbasis, 0), basisobj = samp$basis))
        }else{
            return(samp)
        }
    }else{
        new_sample <- samp
        
        if(second_half){
            new_sample$coefs <- new_sample$coefs[,((n_lam*S)+1):length(new_sample$coefs[1,])]
            new_sample$fdnames$reps <- new_sample$fdnames$reps[((n_lam*S)+1):length(new_sample$fdnames$reps)]
        }else{
            new_sample$coefs <- new_sample$coefs[,1:(n_lam*S)]
            new_sample$fdnames$reps <- new_sample$fdnames$reps[1:(n_lam*S)]
        }
        return(new_sample)
    }
}



evaluate_sample_in_s <- function(s, lambda = 1, samp){
    # returns a subsample for specific value of lambda

    # parameters:
    # s:            argument in which sample will be evaulated
    # lambda:       the proportion of the sample that will be returned
    # samp:         the sample of which a subsample will be returned
  
  
    if (S == 1){
        return(get_lambda_subsample(lambda = lambda, samp = samp))
    }

    n <- length(samp$coefs[1,])/S
    n_lam = floor(n * lambda)

    # if n_lam is zero return empty sample
    if(n_lam == 0){
        return(fd(coef = replicate(samp$basis$nbasis, 0), basisobj = samp$basis))
    }else{
        # evaluation in specific s corresponds to picking every s-th element
        # from the sample, as they are stored as follows:
        # WHOLE SAMPLE = [FIRST_ELM_EVAL_in_s1, FIRST_ELM_EVAL_in_s2,
        #                 FIRST_ELM_EVAL_in_s3, FIRST_ELM_EVAL_in_s4,
        #                 SECOND_ELM_EVAL_in_s1, SECOND_ELM_EVAL_in_s2,
        #                 SECOND_ELM_EVAL_in_s3, ... ]
        # therefore the needed indices (referred to as selected_cols) will be picked
        # (only up to n_lam, since it might be taken a smaller part of the sample)
        selected_cols = replicate(n_lam, 0)
        
        for(i in 1:n_lam){
            selected_cols[i] = (i-1) * S + s
        }
        
        # take the old sample and only keep the previously chosen columns
        subsample <- samp
        subsample$fdnames$reps <- samp$fdnames$reps[selected_cols]
        
        # special case: when n_lam = 1 must be treated differently
        # because of internal fda.R stuff, where rowlabels are not
        # copied and must therefore be copied manually
        if(n_lam == 1){
            coefs <- matrix(samp$coefs[,selected_cols], ncol = 1)
            rownames(coefs) <- rownames(samp$coefs)
            subsample$coefs <- coefs
        }else{
            subsample$coefs <- samp$coefs[,selected_cols]
        }
        
        return(subsample)
    }
}

get_mean <- function(lambda = 1, samp){
    # returns the mean for an fda-object with spatial arguments

    # parameters:
    # lambda:   if specified returns the mean of a subsample of the sample,
    #           corresponding to the lambda-th proportion of it
    # samp:     sample of which the mean will be returned

    if (S == 1){
        return(mean.fd(evaluate_sample_in_s(1, lambda = lambda, samp = samp)))
    }else{
        # the mean function evaluated at s = 1
        mean_1 = mean.fd(evaluate_sample_in_s(1, lambda = lambda, samp = samp))
        mean_coefs = mean_1$coefs
        fdnames_reps = c(mean_1$fdnames$reps, rep(0, S-1))
        
        # append means for s = 2, ..., S
        for(s in 2:S){
          mean_s = mean.fd(evaluate_sample_in_s(s, lambda = lambda, samp = samp))
          mean_coefs <- cbind(mean_coefs, mean_s$coefs[,1])
          fdnames_reps[s] <- mean_s$fdnames$reps
        }

        # copy original sample to keep correct labels of coefficients and fdnames_reps
        # and replace calculated mean coefficients
        sample_mean <- samp
        sample_mean$coefs <- mean_coefs
        sample_mean$fdnames$reps <- fdnames_reps
        
        return(sample_mean)
    }
    
}




# CHANGE POINT SPECIFIC FUNCTIONS #

get_shifted_sample <- function(samp, cp_at = 0.6, factor = 1){

    if(factor == 0){
        return(samp)
    }
    n <- length(samp$coefs[1,])/S
    
    zerof <- matrix(0, nrow=301, ncol=S)
    func_vals <- c()
    
    for(i in 1:floor(cp_at * n)){func_vals <- cbind(func_vals, zerof)}
    for(i in (floor(cp_at * n) + 1):n){
        func_vals <- cbind(func_vals, factor * mean_fnc)
    }

    return(samp + Data2fd(t, func_vals, basisobj = samp$basis))
}

get_cp_estimate <- function(samp){

    n <- length(samp$coefs[1,])/S

    mean_x <- (function(k) get_mean(lambda = k/n, samp = samp))
    mean_y <- (function(k) get_mean(samp = get_lambda_subsample(lambda = k/n, 
                                                                samp = samp,
                                                                second_half = T)))
    
    # mean_x <- (function(k) mean_spatial(get_lambda_subsample(lambda = k/n, samp = samp)))
    # mean_y <- (function(k) mean_spatial(get_lambda_subsample(lambda = k/n, samp = samp, second_half = T)))
    f_k <- (function(k) k * (n-k) / n^2 * norm_spatial(mean_x(k) - mean_y(k))^2)
    
    eps = 0.05
    f_vals <- c()
    
    for(k in 1:n){
        if(k >= (floor(n * eps) + 1) && k <= n - floor(n * eps)){
            f_vals <- append(f_vals, f_k(k))
        }else{
            f_vals <- append(f_vals, 0)
        }
    }
    
    return(1/n * which.max(f_vals))
}



### COVARIANCE MATRIX AND EIGENSYSTEM CALCULATION ###
### 1. HELPER FUNCTIONS FOR COMPUTING EIGENSYSTEMS IN SPATIAL SETUP ###


get_cov_for_s_pair_cp <- function(s1, s2, lambda=1, samp){
  
  # n <- length(samp$coefs[1,])/S
  # cp_est <- get_cp_estimate(samp)
  
  if (CP_est == -1) {
      print("Estimating CP")
      CP_est <<- get_cp_estimate(samp)
      print(paste("...", CP_est))
  }
  
  samp1 <- get_lambda_subsample(lambda = CP_est, samp = samp)
  samp2 <- get_lambda_subsample(lambda = CP_est, samp = samp, second_half = T)
  
  cf1 <- covariance_fnc(evaluate_sample_in_s(s1, lambda, samp1), 
                        evaluate_sample_in_s(s2, lambda, samp1))
  cf2 <- covariance_fnc(evaluate_sample_in_s(s1, lambda, samp2),
                        evaluate_sample_in_s(s2, lambda, samp2))
  
  jointcf <- cf1
  jointcf$coefs <- CP_est * cf1$coefs + (1 - CP_est) * cf2$coefs
  
  return(jointcf)
}

get_cov_for_s_pair <- function(s1, s2, lambda = 1, samp){
    # returns the covariance function only in t_1, t_2 for specific
    # values of s_1 and s_2

    # parameters:
    # s1:       first spatial argument to be evaluated
    # s2:       second spatial argument to be evaluated 
    # lambda:   if specified, returns covariance function only of subsample
    # samp:     the sample of which the covariance function will be returned

    return(covariance_fnc(evaluate_sample_in_s(s1, lambda, samp), 
                          evaluate_sample_in_s(s2, lambda, samp)))
}

covariance_fnc <- function(func1, func2){
    # returns the covariance function of fd-object that only takes values in
    # the time domain. Must be same dimension and actually responds to the cross-
    # covariance function

    # parameters:
    # func1:    first function to take the cross-covariance of
    # func2:    second function to take the cross-covariance of
    
    # read sample size and center functions
    sample_size = dim(func1$coefs)[2]
    c_func1 <- center.fd(func1)
    c_func2 <- center.fd(func2)
    
    final_coef_mat <- 0
    
    # multiply coefficients and divide by sample size
    for(i in 1:sample_size){
        coef_mat <- c_func1$coefs[,i] %*% t(c_func2$coefs[,i])
        rownames(coef_mat) <- colnames(coef_mat)
        
        final_coef_mat <- final_coef_mat + coef_mat
    }
    final_coef_mat <- 1/sample_size * final_coef_mat
    
    # take coefficients to create fd-object in two arguments
    coefnames <- list(dimnames(c_func1$coefs)[[1]], dimnames(c_func2$coefs)[[1]])
    varbifd   <- bifd(final_coef_mat, c_func1$basis, c_func2$basis, coefnames)
    
    return(varbifd)
}

get_Ai <- function(i, lambda = 1, samp){
    # helper function for get_joint_cov_matrix

    A_i <- get_cov_for_s_pair_cp(1, i, lambda = lambda, samp)$coefs
    
    if (S > 1){
        for(s in 2:S){
            A_i <- rbind(A_i, get_cov_for_s_pair_cp(s, i, lambda = lambda, samp)$coefs)
        }
    }

    # check if fourier basis is used, if so then no gram matrix is needed
    # as it is the identity matrix
    if (samp$basis$names[2] == "sin1"){
        return(A_i)
    }else{
        # calculate Gram matrix
        N <- length(samp$coefs[,1])     # = N+1 (for above N)
        all_basis <- fd(coef = diag(1, nrow = N, ncol = N), basisobj = samp$basis) 
        gram_mat <- inprod(all_basis, all_basis)

        return(A_i %*% gram_mat)
    }
}

get_joint_cov_matrix <- function(lambda = 1, samp){
    # assembles the full covariance matrix of which the eigenvalues
    # will be calculated that correspond to the spatial inner product

    # parameters:
    # lambda:   if specified, returns the "joint" covariance matrix of a subsample
    # samp:     sample of which the "joint" covariance matrix will be returned

    final_mat <- get_Ai(1, lambda = lambda, samp)
    
    if (S > 1){
        for(s in 2:S){
            final_mat <- cbind(final_mat, get_Ai(s, lambda = lambda, samp))
        }
    }
    return(final_mat)
}


### 2. ACTUAL COMPUTATION AND STORING OF EIGENSYSTEMS IN SPATIAL SETUP ###

# this is needed only in cp case in order not to compute new eigensystems all the time
# eigensystems will be stored in this list, in order to reuse them
super_eigs <- list()
eigs <- list()
current_eigs <- -1  # stores the cp that this eigs represents


set_eigensystems_by_cp <- function(cp) {
    
    if(current_eigs == cp) {
        return() # no need to change anything since cp is the same
    }
    if(is.null(super_eigs[[paste(cp, sep = "")]])){
        eigs <<- list()
        current_eigs <<- cp
    }else{
        eigs <<- super_eigs[[paste(cp, sep = "")]]
        current_eigs <<- cp
    }
}




calculate_eigensystem <- function(i, samp){
    # computes the specified eigensystem, stores it in list and returns it.
    # If already done it will not calculate it again

    # parameters:
    # i:        corresponds to the value of n_lam, i.e. if the sample has a size of
    #           150 and lambda is equal to 0.5 then n_lam will be 75 and the eigensys.
    #           has number 75. If lambda is equal to 1, then n_lam will be 150 and the
    #           eigensys. has number 150, etc...
    # samp:     sample of which the eigensystem will be computed

    set_eigensystems_by_cp(CP_est)

    # return nothing, i.e. end function call, if nothing has to be done.
    if(is.null(eigs[[paste("vecs_", i, sep="")]]) == FALSE){
        print("Already computed")
        return()
    }
    
    # otherwise compute the corresponding eigensystem.
    # Here the function eigs from the package RSpectra is used
    # as one can choose to compute only the first d eigenvalues/-vectors
    # which will speed up computation
    n <- length(samp$coefs[1,])/S
    
    print(paste("Computing Eigensystem ", i, sep=""))
    t0 <- Sys.time()
    eigen <- eigs(get_joint_cov_matrix(lambda = i/n, samp), k = max(ds), which="LM")
    # eigen <- eigen(getJointMatrixOfCovariances(lambda=i/n))
    
    # store system in global list
    eigs[[paste("vecs_", i, sep="")]] <<- eigen$vectors
    eigs[[paste("vals_", i, sep="")]] <<- eigen$values
    
    # also store list that contains all systems for specific cp in super list
    super_eigs[[paste(CP_est, sep = "")]] <<- eigs

    # print message to let user know when finished
    print(paste("Finished ", i, "/", n,
                " (", round(Sys.time()-t0, digits = 5), " seconds)", sep=""))
}

get_eigensystem_by_nlam <- function(n_lam, samp){
    # returns eigensystem by nlam, based on given sample.
    # If it is not computed yet, it will be computed beforehand

    # parameters:
    # n_lam:    corresponds to the number of the eigensystem, i.e. if the sample has
    #           a size of 150 and lambda is equal to 0.5 then n_lam will be 75 and
    #           the eigensys. has number 75. If lambda is equal to 1, then n_lam will
    #           be 150 and the eigensys. has number 150, etc...
    # samp:     sample of which the eigensystem will be returned

    # set correct eigs for current estimated CP
    set_eigensystems_by_cp(CP_est)
  
    # check if system is already computed, if not then compute
    if(is.null(eigs[[paste("vals_", n_lam, sep="")]])){
        calculate_eigensystem(n_lam, samp)
    }
  
    return(list(vals=eigs[[paste("vals_", n_lam, sep="")]],
                vecs=eigs[[paste("vecs_", n_lam, sep="")]]))
}

### 3. RETRIEVING EIGENVALUES AND -VECTORS ###

get_eigenvalues <- function(range, eig_system){
    # returns range of eigenvalues of specified eigensystem

    # parameters:
    # range:        range of eigenvalues that will be returned,
    #               e.g. 1:10 returns the first ten eigenvalues
    # eig_system:   the eigensystem of which the eigenvalues will be returned from

    return(eig_system$vals[range])
}

evaluate_eigenvector_in_s <- function(i, s, only_realpart = T, 
                                      as_fdobj = T, eig_system, samp){
    # returns the i-th eigenvector evaluated in s

    # parameters:
    # i:                the number of eigenvector
    # s:                spatial argument to be evaulated
    # only_realpart:    if true, then vector will be real valued
    # as_fdobj:         if true, then vector will be fd-object and not list of coefs
    # eig_system:       eigensystem, from which the eigenvector will be taken
    # samp:             sample corresponding to the eigensystem (this is only needed
    #                   for the number of grid points)

    N = length(samp$coefs[,1])     # = N+1 (for above N)
  
    start_index = (s-1) * N + 1
    end_index = s * N
    
    vec <- eig_system$vecs[,i][start_index:end_index]
    vec <- if(only_realpart) Re(vec) else vec
    
    if(as_fdobj){
        vec <- fd(vec, basisobj = samp$basis)
    }
    return(vec)
}

get_eigenvector <- function(i, only_realpart = T, eig_system, normalize = T, samp){
    # returns the i-th eigenvector

    # parameters:
    # i:                the number of eigenvector
    # only_realpart:    if true, then vector will be real valued
    # eig_system:       the eigensystem of which the vector will be taken
    # normalize:        if true, then the vector will have unit spatial norm
    # samp:             sample of which the eigenvector will be returned


    # start assembling the eigenvector, starting with the coefficients for the
    # eigenvector evaulated in s = 1 and store in eigenvec_data
    N = length(samp$coefs[,1])
    eigenvec_data <- c(evaluate_eigenvector_in_s(i, 1, only_realpart = only_realpart,
                                                 as_fdobj = F, eig_system, samp = samp))
    
    if (S > 1){
        # now for s = 2, ..., S
        for(s in 2:S){
            v_s <- evaluate_eigenvector_in_s(i, s, only_realpart = only_realpart, 
                                             as_fdobj = F, eig_system, samp = samp)
            eigenvec_data <- c(eigenvec_data, v_s)
        }
    }
    
    eigvecs_coef <- matrix(eigenvec_data, nrow = N)
    eigvecs <- fd(eigvecs_coef, basisobj = samp$basis) 

    # if normalize, then divide by norm, unless it is zero
    if(normalize){
        nrm <- norm_spatial(eigvecs)
        if(nrm != 0){
            eigvecs <- 1/norm_spatial(eigvecs) * eigvecs
        }
    }

    return(eigvecs)
}










### MATHEMATICAL TECHNICAL FUNCTIONS ###

inprod_spatial <- function(f, g){
    # returns the inner product of f and g, where f and g are functions
    # in L^2 (SxT) and the inner product is
    # sum_s integral f(s,t) * g(s,t) dt

    # parameters:
    # f:    first function of which the inner product will be taken of
    # g:    second function of which the inner product will be taken of

    
    # inprods will be a matrix of all standard inner products (no spatial component)
    # <f(1, .), g(1, .)>  <f(1, .), g(2, .)>  ...  <f(1, .), g(S, .)>
    # <f(2, .), g(1, .)>  <f(2, .), g(2, .)>  ...  <f(2, .), g(S, .)>
    # ...
    # <f(S, .), g(1, .)>  <f(S, .), g(2, .)>  ...  <f(S, .), g(S, .)>
    inprods = inprod(f, g)
    
    # to compute the inner product on L^2(SxT) we only need tu sum up
    # the diagonal entries, where the values of s coincide for f and g
    inprod_sum <- 0
    for(s in 1:S){
        inprod_sum <- inprod_sum + inprods[s,s]
    }
    return(inprod_sum)
}

norm_spatial <- function(f){
    # returns the norm induced by the spatial inner product

    # parameters:
    # f:    function of which the norm will be taken of

    return(sqrt(inprod_spatial(f, f)))
}




### OTHER TECHNICAL FUNCTIONS ###
get_export_matrix <- function(run_num, samp){
    # returns export matrix to be saved in .csv file
    
    # parameters:
    # run_num:  a number for specifying rowlabels
    # samp:     sample which will be used to make these computations
    
    
    support_lam <- seq(1/num_eigsys, (1-1/num_eigsys), by = 1/num_eigsys)
    
    rnames <- c()
    cnames1 = cnames2 = c()
    
    n <- length(samp$coefs[1,])/S
    
    # this will hold results for the first half of the rows (D_n)
    # this will hold results for the second half of the rows (V_n)
    vns = means_diffs = matrix(0, nrow = length(ds), ncol = length(range_data_mean))
    vns_DEL = means_diffs_DEL = matrix(0, nrow = length(ds), ncol = length(range_data_mean))
    
    d_mean_counter <- 1
    
    # for each row (i.e. d) store the d_mean value when cp was correctly est.
    when_cp_correct <- replicate(length(ds), -1)
    
    # for faster computation, remember if cp already has been estimated correctly
    # must be specific for d value. If d=15 was correct, d=20 need not necessarily
    # be also correct, but if data mean increases and d=15 was correct, then
    # for the higher data mean will also be correct
    correct_cp <- replicate(length(ds), FALSE)
    CP_est <<- -1
    prev_cp <- -1
    
    print("Estimating Changepoints...")
    # since the CPs are most of the time the same and have monotonicic properties
    # one can use that to spare time in between:
    # if d=A and d=B (A < B) have the same CP, then all in between will be the same (for const d_mean)
    # d const. and d_mean = A and d_mean = B (A<B) have same CP, all in between will be same
    estimated_cps <- matrix(replicate(length(ds) * length(range_data_mean), -1), nrow = length(ds))
    # estimated_cps <- matrix(replicate(length(ds) * length(range_data_mean), CP_at), nrow = length(ds))
    
    # up to this value will be compared if cps are the same. If so, then all in between
    # can be filled up
    d_mean_special <- range_data_mean[floor(2/3 * length(range_data_mean))]
    gamm <- sqrt(get_gamma(d_mean_special, min(ds), S = S, type = type))
    cp_special <- get_cp_estimate(get_shifted_sample(samp = samp, cp_at = CP_at,
                                                     factor = gamm))
    special_save <- FALSE
    
    d_mean_last <- range_data_mean[length(range_data_mean)]
    gamm <- sqrt(get_gamma(d_mean_last, min(ds), S = S, type = type))
    cp_last <- get_cp_estimate(get_shifted_sample(samp = samp, cp_at = CP_at,
                                                  factor = gamm))
    last_save <- FALSE
    
    # first estimate all cps to save time
    for (d_mean in range_data_mean) {
        
        if (estimated_cps[1, match(d_mean, range_data_mean)] == -1) {
            gamma_low <- sqrt(get_gamma(d_mean, min(ds), S = S, type = type))
            cp_low <- get_cp_estimate(get_shifted_sample(samp = samp, cp_at = CP_at,
                                                         factor = gamma_low))
        } else {
            cp_low <- estimated_cps[1, match(d_mean, range_data_mean)]
        }
      
        if (estimated_cps[length(ds), match(d_mean, range_data_mean)] == -1) {
            gamma_high <- sqrt(get_gamma(d_mean, max(ds), S = S, type = type))
            cp_high <- get_cp_estimate(get_shifted_sample(samp = samp, cp_at = CP_at,
                                                          factor = gamma_high))
        } else {
            cp_high <- estimated_cps[length(ds), match(d_mean, range_data_mean)]
        }
        
        # print(paste("cp_low =", cp_low, "cp_high =", cp_high))
        
        if (d_mean < d_mean_special && special_save == FALSE) {
            if (cp_low == cp_special) {
                # fill all in between
                # print("saved at special")
                for (i in match(d_mean, range_data_mean): (match(d_mean_special, range_data_mean) - 1)) {
                    estimated_cps[1, i] <- cp_low
                }
                special_save <- TRUE
            }

        }else if(last_save == FALSE) {
            if (cp_low == cp_last) {
                # fill all in between
                # print("saved at last")
                for (i in match(d_mean, range_data_mean): (length(range_data_mean) - 1)) {
                    estimated_cps[1, i] <- cp_low
                }
                last_save <- TRUE
            }
        }
      
        
        # Not necessary
        if (cp_low == CP_at) {
            # fill rest of first row as correct cp
            print("from here correct cp (saved time)")
            for (i in match(d_mean, range_data_mean): length(range_data_mean)) {
                estimated_cps[1, i] <- cp_low
            }
        }
        if (cp_high == CP_at) {
            print("from here correct cp (saved time 2)")
            # only fill last row until end with correct cp
            for (i in match(d_mean, range_data_mean): length(range_data_mean)) {
                estimated_cps[length(ds), i] <- cp_high
            }
        }        
        # Not necessary
      
      
        if (cp_low == cp_high) {
            # print("saved time")
            for (d in ds) {
                estimated_cps[match(d, ds), match(d_mean, range_data_mean)] <- cp_low
            }
        } else {
            # print("has to compute new")
            estimated_cps[1, match(d_mean, range_data_mean)] <- cp_low
            estimated_cps[length(ds), match(d_mean, range_data_mean)] <- cp_high
            
            for (j in 2:(length(ds)-1)) {
                gamm <- sqrt(get_gamma(d_mean, d = ds[j], S = S, type = type))
                shifted_samp <- get_shifted_sample(samp = samp, cp_at = CP_at,
                                                   factor = gamm)
                est_cp <- get_cp_estimate(shifted_samp)
                estimated_cps[j, match(d_mean, range_data_mean)] <- est_cp
                
                # print(paste("NEW", est_cp))
            }
        }
        
    }
    
    
    print(estimated_cps)
    print("Finished estimating Changepoints")
    
    # create list with all different cps
    unique_cps <- unique(unlist(as.list(estimated_cps)))

    # BIS HIER WAR OK

    for (ucp in unique_cps) {

        CP_est <<- ucp
        n_cp_est <- floor(n * CP_est)
        print(paste("CP:", CP_est))


        truefalse_mat <- matrix(grepl(ucp, estimated_cps), nrow = length(ds))

        egs <- get_eigensystem_by_nlam(n, samp)

        for (d_mean in range_data_mean) {
            j = match(d_mean, range_data_mean)

            for (d in ds) {
                i = match(d, ds)

                if (truefalse_mat[i, j]) {

                    gamma <- get_gamma(d_mean, d = d, S = S, type = type)

                    # print("A")
                    shift_samp <- get_shifted_sample(samp = samp, cp_at = CP_at,
                                                     factor = sqrt(gamma))

                    samp1 <- get_lambda_subsample(lambda = CP_est, samp = shift_samp)
                    samp2 <- get_lambda_subsample(lambda = CP_est, samp = shift_samp,
                                                  second_half = T)
                    # print("B")

                    m1 <- get_mean(lambda = 1, samp = samp1)
                    m2 <- get_mean(lambda = 1, samp = samp2)

                    # this is dn
                    a_k <- (function(k) inprod_spatial(m1 - m2, 
                                                       get_eigenvector(k, 
                                                                       eig_system = egs,
                                                                       samp = samp))^2)
                    # for (k in 1:d) {
                    #     ev_k <- get_eigenvector(k, eig_system = egs, samp = samp)
                    #     m_diff <- m_diff + inprod_spatial(m1 - m2, ev_k)^2
                    # }
                    m_diff <- Reduce("+", lapply(c(1:d), a_k))
                    # print("C")

                    # this is vn
                    sum_lam <- 0
                    for (lam in support_lam) {

                        corr1 <- floor(lam * n_cp_est) / n_cp_est
                        corr2 <- floor(lam * (n - n_cp_est)) / (n - n_cp_est)

                        # print("C.1")

                        m1_l <- corr1 * get_mean(lambda = lam, samp = samp1)
                        m2_l <- corr2 * get_mean(lambda = lam, samp = samp2)

                        egs_l <- get_eigensystem_by_nlam(floor(n * lam), samp = samp)
                        # print("C.2")

                        b_k <- (function(k) inprod_spatial(m1_l - m2_l, 
                                                           get_eigenvector(k, eig_system = egs_l,
                                                                           samp = samp))^2)
                        
                        left <- Reduce("+", lapply(c(1:d), b_k))
                        right <- lam^2 * m_diff
                        # for (k in 1:d) {
                        #     ev_kl <- get_eigenvector(k, eig_system = egs_l, samp = samp)
                        #     left <- left + inprod_spatial(m1_l - m2_l, ev_kl)^2
                        # }
                        # print("C.4")
                        sum_lam <- sum_lam + (left - right)^2
                        # print("C.4")
                    }
                    
                    #print("D")
                    means_diffs_DEL[i, j] <- m_diff
                    vns_DEL[i, j] <- sqrt(sum_lam / length(support_lam))
                }

            }
        }
    }


















    # AB HIER WURDE GECUTTET
    
    # BIS HIER WURDE GECUTTET

    # row/col labels
    for (d_mean in range_data_mean) {
        cnames1 <- append(cnames1, paste(d_mean, "_mdiff", sep = ""))
        cnames2 <- append(cnames2, paste(d_mean, "_vn", sep = ""))
    }
    
    for (d in ds) {
        rnames <- append(rnames, paste(run_num, "_", d, sep = ""))
    }
    
    fin_mat <- cbind(means_diffs_DEL, vns_DEL)
    
    # also add the last column, which stores the d_mean when cp was correct
    # t(t(.)) makes an (n times 1)-matrix from a list
    fin_mat <- cbind(fin_mat, t(t(when_cp_correct)))
   
    # set row and column labels
    colnames(fin_mat) <- c(cnames1, cnames2, "correct_cp")
    rownames(fin_mat) <- rnames
    
    return(fin_mat)
}



generate_data <- function(num_runs = 200){
    # generates data based on parameters set at the top

    fin_mat = file_paths = c()
    n_max <- max(ns)
    
    dashes <- "------------------------"
    
    for(run in 1:num_runs){
        t0 <- Sys.time()
        print(paste(dashes, " Run ", run, "/", num_runs, " ", dashes, sep = ""))

        # rename sample_max to max_sample or maximum_sample or biggest_sample (largest)
        # put into switch(...)
        if(type == "bm1"){
            sample_max <- generate_brownian_motion_sample1(n = n_max, S = S)
        }else if(type == "nonsep"){
            sample_max <- gen_non_sep(n = n_max, S = S)
        }else if(type == "bm"){
            sample_max <- generate_brownian_motion_sample(n = n_max, S = S)
        }else if(type == "bb"){
            sample_max <- generate_brownian_bridge_sample(n = n_max, S = S)
        }else if(type == "ma"){
            sample_max <- generate_ma_sample(n = n_max, S = S)
        }
        
        for(n_i in ns){
            print(paste("n=", n_i, sep=""))
            samp <- get_lambda_subsample(n_i/n_max, samp = sample_max)

            exp_mat <- get_export_matrix(run_num = run, samp = samp)

            # create data table and remove colname "rn" in 1st col
            dt <- data.table(exp_mat, keep.rownames = T)
            colnames(dt) <- c("", colnames(dt)[2:length(colnames(dt))])
            
            if(run == 1){
                f_name <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S", ".csv")
                file_paths <- append(file_paths, 
                                     paste(data_path, "/", n_i, "/", f_name,
                                           ".csv", sep=""))
                fwrite(dt, file_paths[match(n_i, ns)], row.names = F, col.names = T)
            }else{
                fwrite(dt, file_paths[match(n_i, ns)], row.names = F, 
                       col.names = F, append = T)
            }
            print(paste("Saved data to ", file_paths[match(n_i, ns)],
                        " (", round(Sys.time()-t0, digits = 5), ")", sep = ""))
        }
        # clear up memory
        super_eigs <<- list()
        gc()
        print("")
    }
}



# generate_data(150)