source("../shared/gen_sample.R")
source("../shared/gamma.R")
library(R.utils)        # ?? 
library(data.table)     # for exporting data faster



# parameters

S = 4                         # number of spatial components
num_eigsys = 20               # number of eigensystems
ns = c(150, 250, 500)   # number of sample sizes

CP_at <- 0.6

test = "cp"             # test for cp
type = "bm"             # type of data (bm, ns, ma)

data_path = paste(getwd(), "/data/fully functional/", test, "/", type, sep="")

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
    
    eps = 0.02
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
    range_data_mean <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06,
                         0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15,
                         0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
                         0.25, 0.26, 0.27, 0.28, 0.29, 0.3)

    # initialize matrix which will later be exported. It has columns
    # <_,v_1>, ..., <_,v_d>         where _ will be chosen by type (means or mu)
    # or lambda_1, ..., lambda_d    if type is "vals"
    # length +1 because 1 is not included in the support
    len_matrix <- length(range_data_mean) * 2
    values <- matrix(replicate(len_matrix, 0), ncol = len_matrix)
    
    rnames <- c(paste(run_num, sep = ""))
    cnames1 = cnames2 = c()
    
    n <- length(samp$coefs[1,])/S
   
    
    # this will hold results for the first half of the rows (D_n)
    means_diffs <- c()

    # this will hold results for the second half of the rows (V_n)
    vns <- c()
    
    # for faster computation, remember if cp already has been estimated correctly
    correct_cp <- FALSE
    
    
    for (d_mean in range_data_mean){
      
        gamma <- sqrt(get_gamma_norm(data_mean = d_mean, S = S, type = type))
        shifted_samp <- get_shifted_sample(samp = samp, cp_at = CP_at, factor = gamma)

        # estimate change point (but only if it was incorrect before)
        # if it was correct the estimation becomes only better
        if (correct_cp == FALSE) { 
            cp_est <- get_cp_estimate(samp = shifted_samp)
            if (cp_est == CP_at) { correct_cp <- TRUE }
        }
    
        # add column names
        cnames1 <- append(cnames1, paste(d_mean, "_mdiff", sep = ""))
        cnames2 <- append(cnames2, paste(d_mean, "_vn", sep = ""))

        # this will be needed for both cases and is denoted as 
        # ||D_n(1, theta)||^2 in TeX doc
        # get mean of first estimated part of the sample
        m1 <- get_mean(lambda = cp_est, samp = shifted_samp)
        # get second estimated part of the sample and take mean
        second_samp <- get_lambda_subsample(lambda = cp_est, samp = shifted_samp,
                                            second_half = T)
        m2 <- get_mean(lambda = 1, samp = second_samp)
        means_diff <- norm_spatial(m1 - m2)^2
        means_diffs <- append(means_diffs, means_diff)
        
        # Computation of Vn
        sum_v <- 0
        for (lam in support_lam) {

            corr_factor <- floor(lam * floor(n * cp_est)) / floor(n * cp_est)

            m1_l <- corr_factor * get_mean(lambda = (floor(n * cp_est) * lam / n),
                                           samp = shifted_samp)
            corr_factor2 <- floor(lam * (n - floor(n * cp_est)))/(n - floor(n * cp_est))
            m2_l <- corr_factor2 * get_mean(lambda = lam, samp = second_samp)

            left_side <- norm_spatial(m1_l - m2_l)^2
            right_side <- lam^2 * means_diff

            sum_v <- sum_v + (left_side - right_side)^2
        }

        # append result for specific mean value to vns list
        vns <- append(vns, sqrt(sum_v / length(support_lam)))
    }
    
    # add all entries to matrix
    values[1, ] <- c(means_diffs, vns)
    
    # set row and column labels
    colnames(values) <- c(cnames1, cnames2)
    rownames(values) <- rnames
    return(values)
}

generate_data <- function(num_runs = 200){
    # generates data based on parameters set at the top

    for (n_i in ns) dir.create(paste(data_path, "/", n_i, sep = ""), recursive = TRUE, showWarnings = FALSE)

    fin_mat = file_paths = c()
    n_max <- max(ns)

    for(run in 1:num_runs){
        t0 <- Sys.time()
        print(paste("Run ", run, "/", num_runs, sep = ""))

        # rename sample_max to max_sample or maximum_sample or biggest_sample (largest)
        # put into switch(...)
        if(type == "bm"){
            sample_max <- generate_brownian_motion_sample(n = n_max, S = S)
        }else if(type == "ns"){
            sample_max <- generate_nonseparable_sample(n = n_max, S = S)
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
                        " (time: ", (Sys.time()-t0), ")", sep = ""))
        }
        # clear up memory
        gc()
        print("")
    }
}
