source("file_handling_proj_cp.R")
source("../shared/quantiles_W.R")                 # simulates the limiting distribution
library(latex2exp)                      # for labels in plot


type = "bb"
testing = "cp"


# quantiles (seeded for reproducibility; simulate_w_quantiles() uses parallel::mclapply
# internally, so L'Ecuyer-CMRG + set.seed() is needed for determinism
# regardless of core count)
print("Simulating quantiles...")
RNGkind("L'Ecuyer-CMRG")
set.seed(582017)
quants <- simulate_w_quantiles(alpha = c(0.9, 0.95, 0.99))

q_90 <- quants[1]
q_95 <- quants[2]
q_99 <- quants[3]


test <- function(Delta = 0.02, num_run, data_mean, quantile, d) {

    left_side <- get_mdiff(num_run, data_mean, d)
    right_side <- Delta + (quantile * get_vn(num_run, data_mean, d))

    return(if(left_side > right_side) 1 else 0)
}



analyze_by_sample_size <- function(sample_size = 150, quantile = 95, S = 4, d = 15,
                                   Delta = 0.02, plot = TRUE, add_to_plot = TRUE,
                                   lty = 1){
    
    directory = paste(getwd(), "/data/projection/", testing, "/", type, "/",
                      sample_size, sep = "")
    set_directory(directory)

    read_data()

    q <- 0

    switch(as.character(quantile),
           "90" = {q <- q_90},
           "95" = {q <- q_95},
           "99" = {q <- q_99})

    print(paste("Delta=", Delta, "; n=", sample_size, "; S=", S,
                "; quantile=", quantile, " (", q, ")", sep=""))
    print(paste("Simulation runs:", get_num_of_runs()))

    results = percents = c()
    
    for(actual_mean in x_axis){
        
        res_mean_category <- c()

        for(run in 1:get_num_of_runs()){
            test_result = test(Delta, run, actual_mean, q, d)
            res_mean_category <- append(res_mean_category, test_result)
        }
        results[[paste(actual_mean)]] <- res_mean_category
        
        percent <- Reduce("+", results[[paste(actual_mean)]]) / length(results[[paste(actual_mean)]])
        percents <- append(percents, percent)
        print(paste("ERPr. for norm(mu)=", actual_mean, " is: ", percent, sep=""))
    }


    if(plot){
        x <- x_axis
        y <- percents
        
        if(add_to_plot == FALSE){
            par(mar = c(7.5, 7, 4, 2) + 0.1)
            plot(x = x_axis, y = percents, "l",
                 xlab = "", ylab = "Empirical Rejection Probabilities", lty = lty,
                 ylim = c(0,1), cex.lab = 2.2, cex.axis = 1.8)
            mtext(TeX("$\\sum_{k=1}^d < \\delta, v_k >^2$"), side = 1, line = 7, cex = 2.2)
            # plot(x=x, y=y, type="l")

            abline(v = Delta, lwd = 0.5)
            lines(x = x, y = replicate(length(x), (1-quantile/100)), lwd = 0.5)
        }else{
            lines(x = x, y = percents, lty = lty)
        }
    }
}


analyze_all <- function(save_as_png = F, img_name = "plot.png", Delta = 0.02, d = 15){

    if (save_as_png) { png(img_name) }
    analyze_by_sample_size(sample_size = 150, add_to_plot = F, lty = 1, Delta = Delta, d = d)
    analyze_by_sample_size(sample_size = 250, add_to_plot = T, lty = 2, Delta = Delta, d = d)
    analyze_by_sample_size(sample_size = 500, add_to_plot = T, lty = 3, Delta = Delta, d = d)
    # analyze_by_sample_size(sample_size = 1000, add_to_plot = T, lty = 4)
    legend("topleft", c(expression("n = 150"), expression("n = 250"),
                         expression("n = 500")), lty=1:3, cex = 1.7)
    if (save_as_png) { dev.off() }
}


# New function: like analyze_by_sample_size(), but takes an explicit data
# directory instead of constructing one from the global type/testing
# variables. Used for one-off comparisons (e.g. the nu-sensitivity figure)
# where data lives outside the usual data/projection/cp/<type>/<n>/ layout.
analyze_directory <- function(directory, quantile = 95, d = 15, Delta = 0.15,
                              plot = TRUE, add_to_plot = TRUE, lty = 1){

    set_directory(directory)
    read_data()

    q <- switch(as.character(quantile), "90" = q_90, "95" = q_95, "99" = q_99)

    percents <- sapply(x_axis, function(actual_mean)
        mean(sapply(1:get_num_of_runs(), function(run) test(Delta, run, actual_mean, q, d))))

    if (plot){
        if (!add_to_plot){
            par(mar = c(7.5, 7, 4, 2) + 0.1)
            plot(x = x_axis, y = percents, type = "l", lty = lty, ylim = c(0, 1),
                 xlab = "", ylab = "Empirical Rejection Probabilities",
                 cex.lab = 2.2, cex.axis = 1.8)
            mtext(TeX("$\\sum_{k=1}^d < \\delta, v_k >^2$"), side = 1, line = 7, cex = 2.2)
            abline(v = Delta, lwd = 0.5)
            lines(x = x_axis, y = rep(1 - quantile / 100, length(x_axis)), lwd = 0.5)
        } else {
            lines(x = x_axis, y = percents, lty = lty)
        }
    }

    percents
}


# New function: analyze_all()'s logic, applied to a directory layout of the
# form <path_prefix>/<n>/ for n in 150/250/500, instead of the usual
# data/projection/cp/<type>/<n>/ layout.
analyze_all_directory <- function(path_prefix, d = 15, Delta = 0.15){
    analyze_directory(file.path(path_prefix, "150"), add_to_plot = FALSE, lty = 1, d = d, Delta = Delta)
    analyze_directory(file.path(path_prefix, "250"), add_to_plot = TRUE, lty = 2, d = d, Delta = Delta)
    analyze_directory(file.path(path_prefix, "500"), add_to_plot = TRUE, lty = 3, d = d, Delta = Delta)
    legend("topleft", c(expression("n = 150"), expression("n = 250"),
                        expression("n = 500")), lty = 1:3, cex = 1.7)
}



plot_hist <- function(n = 150, mean_val = 0.05){

    directory = paste(getwd(), "/data/projection/", testing, "/", type, "/",
                      n, sep = "")
    set_directory(directory)

    read_data()

    S = 4

    print(paste("simulation runs:", get_num_of_runs()))

    results <- c()

    for(run in 1:get_num_of_runs()){

        x1 <- get_mdiff(run, mean_val)
        v1 <- get_vn(run, mean_val)

        results <- append(results, (x1 - mean_val)/v1)
        # results <- append(results, v1)
        # results <- append(results, (x1 - mean_val))
        # print(x1-mean_val)
    }
    hist(results, breaks=200) 

}