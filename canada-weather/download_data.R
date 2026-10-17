# Downloads daily climate data from Environment and Climate Change Canada
# (climate.weather.gc.ca) for the station IDs used in canada_application.R,
# and concatenates each station's years into a single CSV under
# data/<id>.csv, matching the format get_file() expects.
#
# Usage: Rscript download_data.R

ids <- c(2205, 2273, 2925, 4333)
years <- 1890:2008   # small buffer around the 1891-2007 range used in canada_application.R

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE)

download_station <- function(id){
    out_file <- file.path(out_dir, paste0(id, ".csv"))
    if (file.exists(out_file)) {
        cat("skip", id, "(already downloaded)\n")
        return(invisible())
    }

    all_years <- list()
    for (year in years) {
        url <- paste0("https://climate.weather.gc.ca/climate_data/bulk_data_e.html",
                       "?format=csv&stationID=", id,
                       "&Year=", year, "&Month=1&Day=1&timeframe=2&submit=Download+Data")
        dat <- tryCatch(read.csv(url, check.names = FALSE), error = function(e) NULL)
        if (!is.null(dat) && nrow(dat) > 0) {
            all_years[[length(all_years) + 1]] <- dat
        }
        Sys.sleep(0.5)   # be polite to the server
    }

    combined <- do.call(rbind, all_years)
    write.csv(combined, out_file, row.names = FALSE)
    cat("wrote", out_file, "(", nrow(combined), "rows )\n")
}

for (id in ids) download_station(id)
