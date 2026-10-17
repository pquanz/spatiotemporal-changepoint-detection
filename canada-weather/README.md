# canada-weather

Real-data application from the paper, testing for a relevant change in the
spatiotemporal mean function of daily mean temperature across 4 Canadian
weather stations.

![Stations used](results/station_map.png)

| Station ID | Name             |
|-----------:|------------------|
| 2205       | Calgary Int'l A  |
| 2273       | Medicine Hat A   |
| 2925       | Indian Head CDA  |
| 4333       | Ottawa CDA       |

Data: daily climate records from Environment and Climate Change Canada
(climate.weather.gc.ca), years 1891-2007 (minus 1910, 1911, 1993, 1995, which
`canada_application.R` excludes).

**Attribution**: based on Environment and Climate Change Canada data, used
under the [ECCC data licence](https://climate.weather.gc.ca/prods_servs/attachment1_e.html).
A trimmed copy (Longitude, Latitude, Date/Time, Mean Temp only, the columns
this code actually uses) is committed in `station-data/`, ~5.2MB total. If
you redistribute this data further, the licence requires you to pass along
the same redistribution terms.

## Usage

Each script in `scripts/` reproduces one result; run from the repo root:

```r
Rscript canada-weather/scripts/changepoint_estimate.R   # theta_hat_n
Rscript canada-weather/scripts/fully_functional_test.R  # Table 1
Rscript canada-weather/scripts/projection_test.R        # Table 2
Rscript canada-weather/scripts/eigenvalue_curve.R        # Figure 7
Rscript canada-weather/scripts/station_map.R
Rscript canada-weather/scripts/yearly_temperature.R
```

All of them work directly against the committed `station-data/`, no
download needed.

To refresh `station-data/` from the live ECCC source instead:

```r
Rscript download_data.R        # full columns, writes to data/ (gitignored)
Rscript trim_station_data.R    # trims to station-data/ (the committed subset)
```

## Files

- `station-data/`: the committed, trimmed station data (see Attribution).
- `download_data.R`: fetches the full station CSVs from ECCC.
- `trim_station_data.R`: trims a full download down to `station-data/`.
- `scripts/`: one script per reproduced result (see Usage above).
- `canada_application.R`: the application itself.
- `gen_data_proj_cp.R`, `gen_sample.R`, `gamma.R`, `quantiles_W.R`: this
  folder's own self-contained copies of the simulation helpers
  `canada_application.R` depends on (kept separate from `../shared/` and
  `../projection/` since this snapshot may have drifted slightly from
  those).
- `results/`: station map and output plots.
