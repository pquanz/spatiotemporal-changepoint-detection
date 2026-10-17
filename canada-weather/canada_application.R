library(fda)
source("quantiles_W.R")
source("gen_data_proj_cp.R")



ids <- c(2205, 2273, 2925, 4333)

years <- seq(1891, 2007)
years <- years[-c(20, 21, 103, 105)]        # herausnehmen der Jahre 1910, 1911, 1993, 1995


# simulate quantiles once (seeded for reproducibility; quantiles_W.R's simulate_w_quantiles()
# uses parallel::mclapply internally, so L'Ecuyer-CMRG + set.seed() is needed
# for deterministic results regardless of core count)
RNGkind("L'Ecuyer-CMRG")
set.seed(319742)
quants <- simulate_w_quantiles(alpha = c(0.9, 0.95, 0.99, 0.975, 0.025))
names(quants) <- c("90", "95", "99", "975", "025")
print(quants)


get_file <- function(id){
    dir <- "station-data/"
    return(read.csv(paste(dir, id, ".csv", sep = "")))
}




get_num_data_points <- function(file, col_number=10){
    
    
    # 10: max temp
    # 11: max temp flag
    # 12: min temp
    # 13: min temp flag
    # 14: mean temp
    # 15: mean temp flag
    # 20: Total Rain (mm)
    # 22: Total Snow (cm)
    # 24: Total Precip (mm)

    return(length(which(!is.na(file[,col_number]))))
}


get_data_for_year <- function(id, year){
    fl <- get_file(id)
    
    mn = which(fl[, "Date.Time"] == paste(year, "-01-01", sep=""))
    mx = which(fl[, "Date.Time"] == paste(year, "-12-31", sep=""))

    # this part is needed to make fdata of same length (365), i.e. schaltjahr weg
    year_data <- (fl[mn:mx, "Mean.Temp...C."])
    if(length(year_data) == 365){
        return(year_data)
    }else{
        return(year_data[-60])
    }
}

create_fdata_from_ids <- function(ids, years){
    library(zoo)
    
    N = 365
    n = 0
    S = length(ids)
    
    vals <- NULL
    
    for(year in years){
        print(year)
        for(id in ids){
            took_to_many_points = F
            
            dat <- get_data_for_year(id, year)
            
            # if data at beginning or end has NA, copy the one in the year before/after and interpolate
            if(is.na(dat[1])){
                more_dat <- get_data_for_year(id, year-1)
                if(is.na(tail(more_dat, n=1))){
                    print("PROBLEM 1")
                }
                dat <- c(tail(more_dat, n=1), dat)
                dataz <- zoo(dat)
                appr <- na.approx(dataz)
                mtrx <- matrix(appr[-1])
                
            }else if(is.na(tail(dat, n=1))){
                more_dat <- get_data_for_year(id, year+1)
                if(is.na(more_dat[1])){
                    print("PROBLEM 2")
                }

                dat <- append(dat, more_dat[1])
                dataz <- zoo(dat)
                appr <- na.approx(dataz)
                mtrx <- matrix(appr[-length(appr)])
                
            }else{
                dataz <- zoo(dat)
                mtrx <- matrix(na.approx(dataz))
            }
            
            
            
            # print(mtrx)
            vals <- cbind(vals, mtrx)
        }
        n <- n + 1
    }
    basis <- create.fourier.basis(nbasis = 149)
    process = Data2fd(argvals = 1:N/N, y = vals, basisobj = basis, lambda = 0)
    
    print(paste("This sample has parameters: N = ", N, ", n = ", n, ", S = ", S, sep = ""))
    
    # return(vals)
    return(process)
}


# hier die sample machen
create_sample <- function(){
    sm <<- create_fdata_from_ids(ids, years)
}







center_fnc <- function(samp, S, n){

    new_s <- samp
    mn <- get_mean(lambda=1, samp=samp)

    for(s in 1:S){

        for(i in 1:n){
            index <- s + (i-1)*S
            # print(index)
            new_s$coefs[, index] <- samp$coefs[, index] - mn$coefs[, s]
        }
    }

    return(new_s)
}







test <- function(samp, n){

    # samp <- samp - get_mean(lambda=1, samp=samp)
    cp_est <- get_cp_estimate(samp=samp)
    print(cp_est)
    # num_eigsys = 10
    
    # erste stat
    first_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = F)
    m1 <- get_mean(lambda = 1, samp = first_samp)
    
    second_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = T)
    m2 <- get_mean(lambda = 1, samp = second_samp)


    LEFTS <- c()
    RIGHTS <- c()


    means_diff <- norm_spatial(m1 - m2)^2
    
    
    # zweite statistik
    support_lam <- seq(1/num_eigsys, (1-1/num_eigsys), by = 1/num_eigsys)
    sum_v <- 0

    sum_v2 <- 0

    for (lam in support_lam) {

        corr_factor <- floor(lam * floor(n * cp_est)) / floor(n * cp_est)
        m1_l <- corr_factor * get_mean(lambda = lam, samp = first_samp)

        corr_factor2 <- floor(lam * (n - floor(n * cp_est)))/(n - floor(n * cp_est))
        m2_l <- corr_factor2 * get_mean(lambda = lam, samp = second_samp)


        el <- floor(lam * (n - floor(n * cp_est)))/(n - floor(n * cp_est))
        er <- floor(lam * floor(n * cp_est)) / floor(n * cp_est)


        # print(paste("C:", norm_spatial(m1_l)^2, ";", norm_spatial(m2_l)^2))
        # print(paste("D:", " LEFT:", el, "RIGHT:", er, "DIFF:", (el - er), "lam", lam))


        left_side2 <- norm_spatial(m1_l - m2_l + (el - er) * get_mean(1, samp))^2
        left_side <- norm_spatial(m1_l - m2_l)^2

        alpha_nl <- el - er
        beta_nl <- left_side2 - left_side

        right_side <- lam^2 * means_diff


        print(paste("AUSGABE:", round(lam, digits=3), "&", round(alpha_nl, digits=4), "&", round(beta_nl, digits=3), "& \\"))
        # print(paste("LETZTE:", round(lam, digits=3), "&", round(inprod_spatial(m1_l - m2_l, get_mean(1, samp)), digits=3), "\\"))

        LEFTS <- append(LEFTS, left_side)
        RIGHTS <- append(RIGHTS, right_side)

        # LEFTS2 <- append(LEFTS2, left_side)
        # RIGHTS2 <- append(RIGHTS2, right_side)


        sum_v <- sum_v + (left_side - right_side)^2
        sum_v2 <- sum_v2 + (left_side2 - right_side)^2



        # print("ZENTRIERT:::")
        # c_samp <- center_fnc(samp)

        # corr_factor <- floor(lam * floor(n * cp_est)) / floor(n * cp_est)
        
        # m1_l <- corr_factor * get_mean(lambda = (floor(n * cp_est) * lam / n),
        #                                samp = c_samp)
        # corr_factor2 <- floor(lam * (n - floor(n * cp_est)))/(n - floor(n * cp_est))
        # m2_l <- corr_factor2 * get_mean(lambda = lam, samp = second_samp)

        # print("VERGLEICH")

        # print(paste("nichtzentr.", left_side, "vs.", ))


    }

    #png("im1.png")
    plot(m1_l - m2_l, ylim=c(-20,20))
    lines(get_mean(1, samp))
    legend("bottomright", legend=c("s=1", "s=2", "s=3", "s=4"), col=c("black", "red", "green", "blue"), pch=c(15,15,15,15))

    #dev.off()

    print(paste("A:", LEFTS))
    print(paste("B:", RIGHTS))
    
    vn <- sqrt(sum_v / length(support_lam))
    vn2 <- sqrt(sum_v2 / length(support_lam))
    
    
    print(means_diff)
    print(vn)

    diff <- vn - vn2

    print(paste("percent increase", diff/vn2))
    print(paste("DF", vn - vn2))


    q_90 <- quants[1]
    q_95 <- quants[2]
    q_99 <- quants[3]
    

    
    print(paste("MAX DELTA (90):", means_diff - q_90 * vn, "MAX DELTA (95):", means_diff - q_95 * vn, "MAX DELTA (99):", means_diff - q_99 * vn))
        
    
    # quantile
    # print("Simulating quantiles...")
    # quants <- simulate_w_quantiles(alpha = c(0.9, 0.95, 0.99))
    
    print(paste("Quant.", q_90, q_95, q_99))
    print("TEST?")
    
    Deltas <- seq(0.01, 1, 0.01)
    
    res <- matrix(0, nrow=length(Deltas), ncol=4)
    
    for(Delta_index in 1:length(Deltas)){
        
        # print(paste("Delta =", Deltas[Delta_index]))
        res[Delta_index, 1] <- Deltas[Delta_index]
        
        # print("90%")
        right_side <- Deltas[Delta_index] + (q_90 * vn)
        # print(means_diff > right_side)
        res[Delta_index, 2] <- means_diff > right_side
        
        
        # print("95%")
        right_side <- Deltas[Delta_index] + (q_95 * vn)
        # print(means_diff > right_side)
        res[Delta_index, 3] <- means_diff > right_side
        
        
        # print("99%")
        right_side <- Deltas[Delta_index] + (q_99 * vn)
        # print(means_diff > right_side)
        res[Delta_index, 4] <- means_diff > right_side
    }

    # print(res, max=100000)
}


### CLEAN REIMPLEMENTATION ###
# relevant_change_test(): computes D_hat = ||D_n(1,cp_est)||^2, V_hat (the
# self-normalizer, 3.8), the maximal threshold still rejected at each level
# (3.9, inverted for Delta), and one-/two-sided confidence intervals for
# ||delta|| (3.11/3.12), following Dette & Quanz (2023) Theorem 3.1/3.2 and
# Remark 3.4 directly. D_n(lambda, cp_est) uses PARTIAL SUMS divided by the
# FULL (not partial) segment size -- see eq. (3.3) -- which is why each term
# below is get_mean() rescaled by (partial count / full count), not a bare
# call to get_mean().
#
# Caveat: the paper does not state what lambda-grid resolution (num_eigsys
# here) it used for the real-data application (only for the simulation
# study). num_eigsys here defaults to the simulation study's own default
# (20); this reproduces the right order of magnitude and qualitative
# behaviour as the paper's reported threshold (2.974) but does not match it
# exactly -- likely due to this or another unstated hyperparameter (e.g.
# basis dimension in create_fdata_from_ids).
relevant_change_test <- function(samp, n, num_eigsys_grid = num_eigsys){

    cp_est <- get_cp_estimate(samp = samp)
    first_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = FALSE)
    second_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = TRUE)

    seg1 <- floor(n * cp_est)
    seg2 <- n - seg1

    D_n <- function(lam) {
        cf1 <- floor(lam * seg1) / seg1
        cf2 <- floor(lam * seg2) / seg2
        m1 <- cf1 * get_mean(lambda = lam, samp = first_samp)
        m2 <- cf2 * get_mean(lambda = lam, samp = second_samp)
        norm_spatial(m1 - m2)^2
    }

    D_hat <- D_n(1)
    support_lam <- seq(1/num_eigsys_grid, (1 - 1/num_eigsys_grid), by = 1/num_eigsys_grid)
    sq_devs <- sapply(support_lam, function(lam) (D_n(lam) - lam^2 * D_hat)^2)
    V_hat <- sqrt(mean(sq_devs))

    max_delta <- c(
        "90%" = D_hat - quants["90"] * V_hat,
        "95%" = D_hat - quants["95"] * V_hat,
        "99%" = D_hat - quants["99"] * V_hat
    )

    ci_one_sided <- c(0, sqrt(D_hat + quants["95"] * V_hat))
    ci_two_sided <- c(sqrt(max(0, D_hat - quants["975"] * V_hat)),
                       sqrt(D_hat + quants["975"] * V_hat))

    list(cp_est = cp_est, D_hat = D_hat, V_hat = V_hat,
         max_delta = max_delta, ci_one_sided = ci_one_sided, ci_two_sided = ci_two_sided)
}


# relevant_change_test_projection(): same idea as relevant_change_test(),
# for the projection-based test (Section 4, Theorem 4.1-4.3). D_n(lambda) is
# projected onto the first d eigenvectors of the lambda-indexed covariance
# estimator (4.3) instead of using the full function norm. d defaults to 51,
# chosen (as in the paper) so e(d) = cumsum(eigenvalues)/trace ~ 0.85 --
# this repo's own eigensystem gives e(51) = 0.8493, matching the paper's
# reported ~0.849 almost exactly, which is good independent evidence this
# part of the construction is right.
relevant_change_test_projection <- function(samp, n, d = 51, num_eigsys_grid = num_eigsys){

    cp_est <- get_cp_estimate(samp = samp)
    first_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = FALSE)
    second_samp <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = TRUE)

    seg1 <- floor(n * cp_est)
    seg2 <- n - seg1

    proj_norm_sq <- function(diff_fn, eig_system, d){
        sum(sapply(1:d, function(k)
            inprod_spatial(diff_fn, get_eigenvector(k, eig_system = eig_system, samp = samp))^2))
    }

    D_n <- function(lam) {
        cf1 <- floor(lam * seg1) / seg1
        cf2 <- floor(lam * seg2) / seg2
        m1 <- cf1 * get_mean(lambda = lam, samp = first_samp)
        m2 <- cf2 * get_mean(lambda = lam, samp = second_samp)
        egs_lam <- get_eigensystem_by_nlam(floor(n * lam), samp)
        proj_norm_sq(m1 - m2, egs_lam, d)
    }

    egs_full <- get_eigensystem_by_nlam(n, samp)
    e_d <- sum(egs_full$vals[1:d]) / sum(egs_full$vals)

    D_hat <- D_n(1)
    support_lam <- seq(1/num_eigsys_grid, (1 - 1/num_eigsys_grid), by = 1/num_eigsys_grid)
    sq_devs <- sapply(support_lam, function(lam) (D_n(lam) - lam^2 * D_hat)^2)
    V_hat <- sqrt(mean(sq_devs))

    max_delta <- c(
        "90%" = D_hat - quants["90"] * V_hat,
        "95%" = D_hat - quants["95"] * V_hat,
        "99%" = D_hat - quants["99"] * V_hat
    )

    ci_one_sided <- c(0, sqrt(D_hat + quants["95"] * V_hat))
    ci_two_sided <- c(sqrt(max(0, D_hat - quants["975"] * V_hat)),
                       sqrt(D_hat + quants["975"] * V_hat))

    list(cp_est = cp_est, d = d, e_d = e_d, D_hat = D_hat, V_hat = V_hat,
         max_delta = max_delta, ci_one_sided = ci_one_sided, ci_two_sided = ci_two_sided)
}



















# ab hier test für projektionen

dees <- c(1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 80, 90, 100)

test_proj <- function(samp, n){
    # quantile
    # print("Simulating quantiles...")
    # quants <- simulate_w_quantiles(alpha = c(0.9, 0.95, 0.99))

    q_90 <- quants[1]
    q_95 <- quants[2]
    q_99 <- quants[3]
    
    print(paste("Quant.", q_90, q_95, q_99))
    
    
    # num_eigsys <- 20
    support_lam <- seq(1/num_eigsys, (1-1/num_eigsys), by = 1/num_eigsys)

    
    cp_est <- get_cp_estimate(samp)
    
    egs <- get_eigensystem_by_nlam(n, samp)
    
    all_res <- NULL
    
    for(d in dees){
        print(paste("d =", d))
        samp1 <- get_lambda_subsample(lambda = cp_est, samp = samp)
        samp2 <- get_lambda_subsample(lambda = cp_est, samp = samp,
                                      second_half = T)
        # print("B")
        
        m1 <- get_mean(lambda = 1, samp = samp1)
        m2 <- get_mean(lambda = 1, samp = samp2)
        
        # this is dn
        a_k <- (function(k) inprod_spatial(m1 - m2, 
                                           get_eigenvector(k, 
                                                           eig_system = egs,
                                                           samp = samp))^2)
        m_diff <- Reduce("+", lapply(c(1:d), a_k))
        
        LEFTS <- c()
        RIGHTS <- c()

        # this is vn
        n_cp_est <- floor(n * cp_est)
        sum_lam <- 0
        for (lam in support_lam) {

            corr1 <- floor(lam * n_cp_est) / n_cp_est
            corr2 <- floor(lam * (n - n_cp_est)) / (n - n_cp_est)

            m1_l <- corr1 * get_mean(lambda = lam, samp = samp1)
            m2_l <- corr2 * get_mean(lambda = lam, samp = samp2)

            egs_l <- get_eigensystem_by_nlam(floor(n * lam), samp = samp)


            # diese konstruktion zu addieren ist getestet und in ordnung.
            b_k <- (function(k) inprod_spatial(m1_l - m2_l, 
                                               get_eigenvector(k, eig_system = egs_l,
                                                               samp = samp))^2)
            left <- Reduce("+", lapply(c(1:d), b_k))
            # ende diese konstruktion


            # print(paste("UNTERSCHIED:", left, "vs.", norm_spatial(m1_l - m2_l)^2 ))

            right <- lam^2 * m_diff

            LEFTS <- append(LEFTS, left)
            RIGHTS <- append(RIGHTS, right)
            # print(paste("A:", left))
            # print(paste("B:", right))

            # print(paste("A:", (left - right)^2))
            sum_lam <- sum_lam + (left - right)^2
        }

        vn <- sqrt(sum_lam / length(support_lam))

        print(paste("A:", LEFTS))
        print(paste("B:", RIGHTS))

        # vn <- 0.8

        print(m_diff)
        print(vn)
        print(paste("MAX DELTA (90):", m_diff - q_90 * vn, "MAX DELTA (95):", m_diff - q_95 * vn, "MAX DELTA (99):", m_diff - q_99 * vn))
        
        
        print("TEST?")
        
        Deltas <- seq(0.01, 2, 0.001)
        
        res <- matrix(0, nrow=length(Deltas), ncol=4)
        
        for(Delta_index in 1:length(Deltas)){
            
            # print(paste("Delta =", Deltas[Delta_index]))
            res[Delta_index, 1] <- Deltas[Delta_index]
            
            # print("90%")
            right_side <- Deltas[Delta_index] + (q_90 * vn)
            # print(means_diff > right_side)
            res[Delta_index, 2] <- m_diff > right_side
            
            
            # print("95%")
            right_side <- Deltas[Delta_index] + (q_95 * vn)
            # print(means_diff > right_side)
            res[Delta_index, 3] <- m_diff > right_side
            
            
            # print("99%")
            right_side <- Deltas[Delta_index] + (q_99 * vn)
            # print(means_diff > right_side)
            res[Delta_index, 4] <- m_diff > right_side
        }
        
        all_res <- cbind(all_res, res)
        # print(res, max=100000)
    }
    # print(all_res, max=1000000)
    
    write.table(all_res, file="mymatrix.txt", row.names=FALSE, col.names=FALSE)
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
    # set_eigensystems_by_cp(CP_est)
    
    # check if system is already computed, if not then compute
    if(is.null(eigs[[paste("vals_", n_lam, sep="")]])){
        calculate_eigensystem(n_lam, samp)
    }
    
    return(list(vals=eigs[[paste("vals_", n_lam, sep="")]],
                vecs=eigs[[paste("vecs_", n_lam, sep="")]]))
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
    
    # set_eigensystems_by_cp(CP_est)
    
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
    eigen <- eigs(get_joint_cov_matrix(lambda = i/n, samp), k = max(dees), which="LM")
    # eigen <- eigen(getJointMatrixOfCovariances(lambda=i/n))
    
    # store system in global list
    eigs[[paste("vecs_", i, sep="")]] <<- eigen$vectors
    eigs[[paste("vals_", i, sep="")]] <<- eigen$values
    
    # also store list that contains all systems for specific cp in super list
    # super_eigs[[paste(CP_est, sep = "")]] <<- eigs
    
    # print message to let user know when finished
    print(paste("Finished ", i, "/", n,
                " (", round(Sys.time()-t0, digits = 5), " seconds)", sep=""))
}




get_cov_for_s_pair_cp <- function(s1, s2, lambda=1, samp){
    
    # n <- length(samp$coefs[1,])/S
    # cp_est <- get_cp_estimate(samp)
    
    cp_est <- get_cp_estimate(samp)
    
    samp1 <- get_lambda_subsample(lambda = cp_est, samp = samp)
    samp2 <- get_lambda_subsample(lambda = cp_est, samp = samp, second_half = T)
    
    cf1 <- covariance_fnc(evaluate_sample_in_s(s1, lambda, samp1), 
                          evaluate_sample_in_s(s2, lambda, samp1))
    cf2 <- covariance_fnc(evaluate_sample_in_s(s1, lambda, samp2),
                          evaluate_sample_in_s(s2, lambda, samp2))
    
    jointcf <- cf1
    jointcf$coefs <- cp_est * cf1$coefs + (1 - cp_est) * cf2$coefs
    
    return(jointcf)
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
    # print(paste("BLA", "n=",n))
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
        # print(paste("1:", length(samp$fdnames$reps), "2:", selected_cols))
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


# # n=250

# 0.05 & -0.0057 & 0.025 & \\
# 0.1 & -0.0038 & 0.045 & \\
# 0.15 & -0.0019 & 0.039 & \\
# 0.2 & 0 & 0 & \\
# 0.25 & -0.0057 & 0.195 & \\
# 0.3 & -0.0038 & 0.154 & \\
# 0.35 & -0.0019 & 0.091 & \\
# 0.4 & 0 & 0 & \\
# 0.45 & -0.0057 & 0.333 & 0.242 \\
# 0.5 & -0.0038 & 0.249 & \\
# 0.55 & -0.0019 & 0.139 & \\
# 0.6 & 0 & 0 & \\
# 0.65 & -0.0057 & 0.475 & \\
# 0.7 & -0.0038 & 0.347 & \\
# 0.75 & -0.0019 & 0.191 & \\
# 0.8 & 0 & 0 & \\
# 0.85 & -0.0057 & 0.623 & \\
# 0.9 & -0.0038 & 0.445 & \\
# 0.95 & -0.0019 & 0.241 & \\




# # n=500

# 0.05 & -0.0019 & 0.008 & \\
# 0.1 & 0 & 0 & \\
# 0.15 & -0.0019 & 0.032 & \\
# 0.2 & 0 & 0 & \\
# 0.25 & -0.0019 & 0.058 & \\
# 0.3 & 0 & 0 & \\
# 0.35 & -0.0019 & 0.086 & \\
# 0.4 & 0 & 0 & \\
# 0.45 & -0.0019 & 0.11 & 0.084 \\
# 0.5 & 0 & 0 & \\
# 0.55 & -0.0019 & 0.137 & \\
# 0.6 & 0 & 0 & \\
# 0.65 & -0.0019 & 0.163 & \\
# 0.7 & 0 & 0 & \\
# 0.75 & -0.0019 & 0.191 & \\
# 0.8 & 0 & 0 & \\
# 0.85 & -0.0019 & 0.217 & \\
# 0.9 & 0 & 0 & \\
# 0.95 & -0.0019 & 0.243 & \\


# # n = 1000

# 0.05 & 0 & 0 & \\
# 0.1 & 0 & 0 & \\
# 0.15 & 0 & 0 & \\
# 0.2 & 0 & 0 & \\
# 0.25 & 0 & 0 & \\
# 0.3 & 0 & 0 & \\
# 0.35 & 0 & 0 & \\
# 0.4 & 0 & 0 & \\
# 0.45 & 0 & 0 & 0 \\
# 0.5 & 0 & 0 & \\
# 0.55 & 0 & 0 & \\
# 0.6 & 0 & 0 & \\
# 0.65 & 0 & 0 & \\
# 0.7 & 0 & 0 & \\
# 0.75 & 0 & 0 & \\
# 0.8 & 0 & 0 & \\
# 0.85 & 0 & 0 & \\
# 0.9 & 0 & 0 & \\
# 0.95 & 0 & 0 & \\