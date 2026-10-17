# Reproduces results/station_map.png: a map of Canada showing the 4
# weather stations used in canada_application.R.
#
# Usage (from repo root): Rscript canada-weather/scripts/station_map.R
#
# Reads coordinates from the committed station-data/ (not the gitignored
# full data/ download), so this runs with no download step.

setwd("canada-weather")

library(maps)

stations <- list(
    "2205" = "Calgary Int'l A",
    "2273" = "Medicine Hat A",
    "2925" = "Indian Head CDA",
    "4333" = "Ottawa CDA"
)

get_coords <- function(id){
    fl <- read.csv(file.path("station-data", paste0(id, ".csv")))
    c(fl[1, "Longitude..x."], fl[1, "Latitude..y."])
}

dir.create("results", showWarnings = FALSE)
png("results/station_map.png", width = 800, height = 600)

maps::map("world", "Canada", xlim = c(-140, -55), ylim = c(42, 70), fill = TRUE, col = "grey95")

for (id in names(stations)){
    coords <- get_coords(id)
    points(coords[1], coords[2], col = "red", pch = 19, cex = 1.4)
    text(coords[1], coords[2], labels = stations[[id]], pos = 4, cex = 1.2)
}

title("Weather stations used in the real-data application", cex.main = 1.8)
dev.off()
cat("wrote results/station_map.png\n")
