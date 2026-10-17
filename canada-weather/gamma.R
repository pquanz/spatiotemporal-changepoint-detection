# This file solely purposes the calculation of specific values
# for gamma in the mean function


hyp = "proj"    # the hypothesis that should be tested




get_gamma <- function(data_mean, d, S, type = "bm"){

    # for fast return
    if (data_mean == 0){
        return (0)
    }

    # these are proj. hyp
    switch(type, 
        "bm" = {
            # gives the sequence 1/pi, 1/pi, 1/(3*pi), 1/(3*pi), 1/(5*pi), ...
            # second row takes the sum of a_k^2 up to d
            a_k <- (function(k) 1/pi * 1/(k - (k %% 2 == 0)))
            sum_sq_ak <- Reduce("+", a_k(seq(1:d))^2)

            # square of sum of 1,2,3,4,5,...,n
            sq_sum_s <- Reduce("+", c(1:S))^2

            return(data_mean * S / (2 * sq_sum_s * sum_sq_ak))
        },
        "bm1" = {       # this is the latest
            # gives the sequence 1/pi, 1/pi, 1/(3*pi), 1/(3*pi), 1/(5*pi), ...
            # second row takes the sum of a_k^2 up to d
            a_k <- (function(k) 1/pi * 1/(k - (k %% 2 == 0)))
            sum_sq_ak <- Reduce("+", a_k(seq(1:d))^2)

            # square of sum of 1,2,3,4,5,...,n
            sq_sum_s <- Reduce("+", c(1:S)^2)^2

            return(15 * data_mean / (sq_sum_s * sum_sq_ak))
        },
        "nonsep" = {       # this is the latest
            # gives the sequence 1/pi, 1/pi, 1/(3*pi), 1/(3*pi), 1/(5*pi), ...
            # second row takes the sum of a_k^2 up to d
            s0 <- S
            s1 <- Reduce("+", (seq(1:S)))
            s2 <- Reduce("+", (seq(1:S))^2)
            norm_factor <- (function(k) 1/sqrt(2) * sqrt(s0 + s2))
            a_k <- (function(k) 1/norm_factor(k) * (s1 * 8 * k / (pi * (16 * k^2 - 1)) + s2 * 2/(pi * (1-16*k^2)) ) )
            # a_k <- (function(k) 1/())
            # a_k <- (function(k) 1/pi * 1/(k - (k %% 2 == 0)))
            sum_sq_ak <- Reduce("+", a_k(seq(1:d))^2)

            # square of sum of 1,2,3,4,5,...,n
            # sq_sum_s <- Reduce("+", c(1:S)^2)^2

            return(data_mean / (sum_sq_ak))
        },
        "bb" = {
            # second row takes the sum of a_k^2 up to d
            a_k <- (function(k) 1/pi * (2 * (sin(pi * k) - 2*k) ) / (1 - 4*k^2))
            sum_sq_ak <- Reduce("+", a_k(seq(1:d))^2)

            # square of sum of 1,2,3,4,5,...,n
            sq_sum_s <- Reduce("+", c(1:S))^2

            return(data_mean * S / (2 * sq_sum_s * sum_sq_ak))
        },
        "ma" = {
            a_k <- (function(k) 1/pi * 2*(sin(2*pi*k) - 4*k)/(1 - 16 * k^2))
            sum_sq_ak <- Reduce("+", a_k(seq(1:d))^2)
            
            sum_sq_s <- Reduce("+", c(1:S)^2)

            return(data_mean / (2 * sum_sq_s * sum_sq_ak))
        })
}



get_gamma_norm <- function(data_mean, S, type = "bm"){
    
    if (data_mean == 0){
        return (0)
    }
    
    # these are norm hyp
    switch(type, 
           "bm" = {
               sum_sq_s <- Reduce("+", c(1:S)^2)
               return(2 * data_mean / sum_sq_s)
           },
           "bm1" = {
               sum_sq_s <- Reduce("+", c(1:S)^2)
               return(2 * data_mean / sum_sq_s)
           },
           "nonsep" = {
               sum_sq_s <- Reduce("+", c(1:S)^2)
               return(2 * data_mean / sum_sq_s)
           },
           "bb" = {
               sum_sq_s <- Reduce("+", c(1:S)^2)
               return(2 * data_mean / sum_sq_s)
           },
           "ma" = {
               sum_sq_s <- Reduce("+", c(1:S)^2)
               return(2 * data_mean / sum_sq_s)
           })
}