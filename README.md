# spatiotemporal-changepoint-detection

Simulation and application code for:

Dette, H. and Quanz, P. (2023). Detecting relevant changes in the spatiotemporal
mean function. *Journal of Time Series Analysis*, 44(5-6), 505-532.

Tests whether the mean function of a spatiotemporal process has changed by more
than a given threshold (a "relevant" change), using projections onto a
finite-dimensional basis.

## Quick start

```r
setwd("projection")
source("gen_data_proj_cp.R")      # sets parameters, defines generate_data()
generate_data(num_runs = 50)      # writes CSVs to data/<n>/
source("data_analysis_proj_cp.R") # defines analyze_by_sample_size()
analyze_by_sample_size(sample_size = 500, d = 15, plot = TRUE)
```

`fully-functional/` follows the same pattern (`gen_data_cp.R` +
`data_analysis_ff_cp.R`, no `d` argument). For many runs in parallel, see
`cluster/run_parallel.R`. For the real-data application, see
`canada-weather/README.md`.

## The test

Let $X_1(t,s),\dots,X_n(t,s)$, $t\in[0,1]$, $s\in S$, be the observed
curves, so that $X_1,\dots,X_{k_0}$ have mean $\mu_1(t,s)$ and
$X_{k_0+1},\dots,X_n$ have mean $\mu_2(t,s)$ for some unknown change point
$k_0$. Write $\delta(t,s) = \mu_2(t,s) - \mu_1(t,s)$ for the (unknown)
change and fix a relevance threshold $\Delta>0$. We test

$$H_0: \|\delta\|^2 \le \Delta \quad \text{against} \quad H_1: \|\delta\|^2 > \Delta$$

With a change-point estimate $\hat k_n$ and $\lambda\in(0,1)$, define

$$D_n(\lambda)(t,s) = \frac{1}{\hat k_n}\sum_{i=1}^{\lfloor \lambda \hat k_n \rfloor} X_i(t,s) - \frac{1}{n-\hat k_n}\sum_{j=\hat k_n+1}^{\hat k_n+\lfloor \lambda(n-\hat k_n) \rfloor} X_j(t,s)$$

We estimate $\|\delta\|^2$ by $\hat{\mathbb D}_n = \|D_n(1)\|^2$. We use
self-normalization via

$$\hat{\mathbb V}_n = \left(\int_0^1 \left(\|D_n(\lambda)\|^2 - \lambda^2 \hat{\mathbb D}_n\right)^2 \, d\nu(\lambda)\right)^{1/2}$$

for a probability measure $\nu$ on $(0,1)$ (our code discretizes this
integral on a grid; the paper's own definition has no grid, just the
integral); $\hat{\mathbb V}_n$ needs no separate estimate of the error
process's variance. We show that the decision rule, reject $H_0$ if

$$\hat{\mathbb D}_n > \Delta + q_{1-\alpha}(\mathbb W)\,\hat{\mathbb V}_n$$

has asymptotic level $\alpha$, where $q_{1-\alpha}(\mathbb W)$ is the
quantile of

$$\mathbb W = \frac{\mathbb B(1)}{\left(\int_0^1 \lambda^2(\mathbb B(\lambda) - \lambda \mathbb B(1))^2 \, d\nu(\lambda)\right)^{1/2}}$$

with $\mathbb B$ a standard Brownian motion on $[0,1]$.

The projection-based test applies the same idea after replacing $\delta$ by
its projection onto the first $d$ eigenfunctions of the error process's
covariance, motivated by the Karhunen-Loève decomposition: truncating to
$d$ terms still captures most of the process's variance (we pick $d$ so the
explained-variance ratio $e(d) \approx 0.85$, as in the paper) while
keeping the test's reference distribution finite-dimensional.

## Results

Full writeup, all three error processes, and the data example in
[`results/README.md`](results/README.md). Highlights:

![Projection power curve](results/projection/power-curve-3process.png)
![Fully-functional power curve](results/fully-functional/power-curve-3process.png)

*Empirical rejection probability vs. the true squared norm of the change,
one panel per error process (Brownian motion / non-separable / fMA(1), as
in the paper), for the projection test (top, d=15/11/11) and the
fully-functional test (bottom). Both cross close to the nominal 5% level at
the threshold $\Delta=0.15$.*

The fully-functional test uses the whole function (the "usual" $L^2$ norm);
the projection test reduces to its first $d$ components instead. $n=500$,
Brownian-motion process:

| $\|\delta\|^2$ | 0 | 0.1 | **0.15 ($\Delta$)** | 0.2 | 0.3 |
|---|---:|---:|---:|---:|---:|
| Projection (d=15) | 0.000 | 0.005 | **0.045** | 0.170 | 0.700 |
| Fully-functional (full $L^2$) | 0.000 | 0.000 | **0.075** | 0.275 | 0.850 |

**Canada weather data example**: the change-point estimator finds a relevant
change in mean temperature at **$\hat\vartheta_n=0.7079646$** (matching
the paper's reported 0.708 almost exactly), around 1971/1972; see
`results/README.md` for how this differs from the paper's own stated
"~1969".

![Stations used](canada-weather/results/station_map.png)
![Yearly mean temperature](canada-weather/results/yearly_temperature.png)

See `results/README.md` for the full threshold/confidence-interval results
and the station details.
