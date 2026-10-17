# Reproduces results/yearly_temperature.png: yearly mean daily temperature
# at the 4 stations, with the estimated change point marked.
#
# Usage (from repo root): Rscript canada-weather/scripts/yearly_temperature.R

setwd("canada-weather")
source("canada_application.R")

station_names <- c("Calgary Intl A", "Medicine Hat A", "Indian Head CDA", "Ottawa CDA")
yearly_means <- matrix(NA, nrow = length(years), ncol = length(ids))
for (i in seq_along(ids)) {
    for (j in seq_along(years)) {
        yearly_means[j, i] <- mean(get_data_for_year(ids[i], years[j]), na.rm = TRUE)
    }
}

create_sample()
cp_est <- get_cp_estimate(samp = sm)
cp_year <- years[floor(length(years) * cp_est)]

# Overall average across the 4 stations, before/after the change point (the
# paper's model (2.1) convention: cp_year is the LAST pre-change year).
overall_mean <- rowMeans(yearly_means, na.rm = TRUE)
pre_mean <- mean(overall_mean[years <= cp_year])
post_mean <- mean(overall_mean[years > cp_year])

dir.create("results", showWarnings = FALSE)
png("results/yearly_temperature.png", width = 1000, height = 600)
par(mar = c(5.5, 7, 4, 2) + 0.1)
matplot(years, yearly_means, type = "l", lty = 1, col = 1:4, lwd = 1.5,
        xlab = "Year", ylab = "Mean daily temperature (°C)",
        main = "Yearly mean temperature at the 4 Canadian stations",
        cex.lab = 2.2, cex.axis = 1.8, cex.main = 1.8)
abline(v = cp_year, col = "darkgray", lwd = 2, lty = 2)
segments(min(years), pre_mean, cp_year, pre_mean, lwd = 4, col = "purple")
segments(cp_year, post_mean, max(years), post_mean, lwd = 4, col = "purple")
text(cp_year, max(yearly_means, na.rm = TRUE), paste("estimated change point (", cp_year, ")", sep = ""),
     pos = 4, col = "darkgray", cex = 1.2)
legend("bottomright", legend = c(station_names, "Overall average (pre/post CP)"),
       col = c(1:4, "purple"), lty = 1, lwd = c(1.5, 1.5, 1.5, 1.5, 4), bty = "n", cex = 1.7)
dev.off()
cat("wrote results/yearly_temperature.png\n")
