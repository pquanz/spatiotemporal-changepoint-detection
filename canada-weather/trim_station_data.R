# Trims the full ECCC station downloads (data/, gitignored, from
# download_data.R) down to just the columns this repo's code actually uses
# (Longitude, Latitude, Date/Time, Mean Temp), and writes them to
# station-data/, which IS committed to git (~5.5MB total for 4 stations,
# vs. ~22MB for the full 31-column downloads).
#
# Usage: Rscript trim_station_data.R

ids <- c(2205, 2273, 2925, 4333)
cols <- c("Longitude..x.", "Latitude..y.", "Date.Time", "Mean.Temp...C.")

dir.create("station-data", showWarnings = FALSE)

for (id in ids) {
    fl <- read.csv(file.path("data", paste0(id, ".csv")))
    write.csv(fl[, cols], file.path("station-data", paste0(id, ".csv")), row.names = FALSE)
    cat("wrote station-data/", id, ".csv\n", sep = "")
}
