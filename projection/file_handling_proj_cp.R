library(data.table)








x_axis <- c()
m_diffs <- c()
vns <- c()
ds <- c()



read_data <- function(){
    
    read_x_axis <- FALSE        # will be true after first loop run

    m_diffs_temp <- c()
    vns_temp <- c()

    for (path in files) {
        data <- fread(path)

        half <- ((length(data[1,]) - 1) / 2)        # marks the end of the D_n's

        # read various data means (i.e. 0, 0.005, 0.01, 0.015, 0.02, ..., 0.15)
        if (read_x_axis == FALSE) {
            temp_x_axis <- colnames(data)[2:(half + 1)]
            x_axis <<- as.numeric(gsub("_mdiff", "", temp_x_axis))

            # also read ds
            TEMP_d <- which(grepl(paste("^1_", sep=""), data[["V1"]], fixed = FALSE))
            ds <<- as.integer(gsub("._", "", data[["V1"]][1:length(TEMP_d)]))

            read_x_axis <- TRUE     # only read x-axis once, so set to true
        }

        for(x_val in x_axis) {
            m_diffs_temp[[paste(x_val)]] <- append(m_diffs_temp[[paste(x_val)]], 
                                                   data[[paste(x_val, "_mdiff",
                                                               sep = "")]])
            vns_temp[[paste(x_val)]] <- append(vns_temp[[paste(x_val)]], 
                                               data[[paste(x_val, "_vn", sep = "")]])
        }
    }

    m_diffs <<- m_diffs_temp
    vns <<- vns_temp
}




files <- c()

set_directory <- function(path){
    files <<- c()
    files <<- list.files(path, pattern="*.csv", full.names=TRUE)
    print(paste("Found", length(files), "files"))
}

get_num_of_runs <- function() {
    return(length(vns[[paste(x_axis[1], sep = "")]])/length(ds))
}

get_vn <- function(num_run, data_mean, d) {
    if (d %in% ds){
        index <- length(ds) * (num_run - 1) + match(d, ds)
        return(vns[[paste(data_mean, sep = "")]][index])
    }
    print("no possible d value")
}

get_mdiff <- function(num_run, data_mean, d) {
    if (d %in% ds){
        index <- length(ds) * (num_run - 1) + match(d, ds)
        return(m_diffs[[paste(data_mean, sep = "")]][index])
    }
    print("no possible d value")
    
}