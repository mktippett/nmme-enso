# Reproducing NOAA CPC's Relative ONI (RONI) from ERSSTv6

This note investigates how to reproduce NOAA CPC's operational **Relative
Oceanic Niño Index (RONI)** — both its monthly and seasonal published
series — from the published ERSSTv6 SST dataset. CPC's RONI is a distinct
product from the relative Niño-3.4 index this project computes for NMME
forecast models (`docs/relative_nino34.md`); it is examined here purely to
reproduce an external, independently verifiable series. CPC's own
description of RONI fixes the anomaly recipe (SST source, box definitions,
base period, 3-month running mean) but not the numeric scaling factor that
"adjusts \[the difference\] so the variance equals the original Niño 3.4
index." The factor is not published, and CPC's description does not state
what period, or whether detrending, entered its computation. The
motivating question is therefore: given only ERSSTv6 and CPC's stated
recipe, can the published RONI values — monthly and seasonal — be
reproduced?

## 1. Data and the published recipe

CPC publishes RONI at two temporal resolutions:

- **Monthly** relative Niño-3.4 anomaly, `Rnino34.ascii.txt`
  (https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt) — one value
  per calendar month, no smoothing.
- **Seasonal** RONI, `RONI.ascii.txt`
  (https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt) — a 3-month
  running mean, labeled with the standard overlapping season code (DJF,
  JFM, ..., NDJ) and the calendar year of the season's center month.

Both are stated to derive from ERSSTv6
(https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc):
the Niño-3.4 box (5°S–5°N, 170°W–120°W) minus the tropical-mean box
(20°S–20°N, all longitudes), each anomaly relative to a 1991–2020 monthly
climatology. Cosine-latitude weighting versus a plain grid-cell mean makes
a negligible difference to the results below (under 0.0001 in every fitted
factor), so the choice between them is not otherwise a source of
uncertainty here.

CPC's description amounts to the following formula, at monthly resolution:

$$
\text{RONI}(\text{month}) \;=\; f(\text{month}) \times \Big(\underbrace{\text{N34}_{\text{anom}}}_{\text{known}} \;-\; \underbrace{\text{Trop}_{\text{anom}}}_{\text{known}}\Big)
$$

$\text{N34}_{\text{anom}}$ and $\text{Trop}_{\text{anom}}$ — the two box
anomalies above — are both computable directly from ERSSTv6. The
per-calendar-month factor $f(\text{month})$ is the one unknown quantity:
CPC states only that it is chosen so that "the variance equals the
original Niño 3.4 index," not its numeric value, the period over which
that variance is computed, or whether the underlying series is
detrended first. Backing out $f(\text{month})$ from the known left- and
right-hand sides is the task of Sections 3–4.

## 2. The published monthly and seasonal series are mutually consistent

Before fitting anything, the two published series were checked against
each other directly: the 3-month centered running mean of the published
*monthly* series was compared with the published *seasonal* series, with
no involvement of ERSSTv6 or any fitted factor. All 919 published seasonal
values (1950–2026) matched a corresponding monthly-average value. The
maximum discrepancy is 0.0067°C and the mean discrepancy is 0.00007°C,
indistinguishable from zero. A discrepancy of exactly 0.0067°C (two-thirds
of the 0.01°C rounding step) is precisely the bound expected from two
series that are each independently rounded to two decimal places before
publication, which means that no unmodeled processing separates the two
products.

This confirms that CPC's scaling factor is applied at the **monthly**
level, with the published seasonal values simply the 3-month average of
already-scaled monthly values — not a separately-fit seasonal quantity. An
important consequence is that a season's effective factor is *not* the
simple average of its three constituent monthly factors. Within a season,
the three months' unscaled differences are serially correlated and
generally unequal in magnitude (ENSO SST anomalies evolve on a time scale
of several months), so a factor fit directly to 3-month-averaged data
recovers a variance-weighted combination of the three monthly factors,
not their arithmetic mean. Sections 3 and 4 quantify how much this
distinction matters in practice.

## 3. Backing out the monthly factors

Each calendar month's factor was fit independently by ordinary least
squares, forced through the origin:

$$\hat f(\text{month}) = \frac{\sum_i x_i y_i}{\sum_i x_i^2}, \qquad x_i = \text{ERSSTv6 monthly diff (Niño-3.4 anom} - \text{tropical-mean anom)}, \quad y_i = \text{published monthly RONI}$$

using the 20 most recent years of each calendar month's own record
(2006/2007–2025/2026, depending on data availability at the time of
writing). The 12 factors, with their standard errors, are:

| Month | Factor | SE | Month | Factor | SE |
|---|---|---|---|---|---|
| Jan | 1.248 | 0.0013 | Jul | 1.187 | 0.0025 |
| Feb | 1.274 | 0.0016 | Aug | 1.167 | 0.0017 |
| Mar | 1.339 | 0.0025 | Sep | 1.202 | 0.0017 |
| Apr | 1.390 | 0.0027 | Oct | 1.215 | 0.0013 |
| May | 1.356 | 0.0026 | Nov | 1.221 | 0.0012 |
| Jun | 1.248 | 0.0034 | Dec | 1.232 | 0.0013 |

The factor follows a pronounced annual cycle, peaking in April (1.390) and
reaching a minimum in August (1.167) — a range of 0.22, far larger than
any of the standard errors above (at most 0.0034). It is therefore a
genuinely month-dependent quantity, not a constant obscured by sampling
noise.

Applying these factors reproduces both published series with a small,
fairly uniform residual. Pooled across all 240 fitted (month, year) points,
the RMS error against the published monthly series is 0.0054°C, with a
modest systematic bias of −0.0038°C (the fit averages about 0.004°C above
the published values). Applying the same monthly factors and then taking a
3-month running mean — the same operation CPC applies, per Section 2 —
reproduces the published seasonal series with an overall RMS of 0.0049°C,
essentially unchanged from the monthly figure. Aggregating to seasons
neither helps nor hurts the fit, which is expected once Section 2
established that seasonal values are nothing more than the monthly
average.

The summary statistics above describe the fit's accuracy in aggregate;
Tables 2 and 3 show what that accuracy looks like against the actual
published values, for the same 2020–2026 window as CPC's own RONI product
page (https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/),
which displays these values rounded to one decimal place — the two-decimal
values below are from the underlying ASCII files, not visible on that page
directly, and are what make the small remaining errors visible at all.

*Table 2. Published seasonal RONI, 2020–2026 (`RONI.ascii.txt`).*

| Year | DJF | JFM | FMA | MAM | AMJ | MJJ | JJA | JAS | ASO | SON | OND | NDJ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2020 | 0.14 | 0.10 | 0.02 | -0.25 | -0.55 | -0.72 | -0.76 | -0.85 | -1.16 | -1.37 | -1.46 | -1.39 |
| 2021 | -1.32 | -1.18 | -1.08 | -0.89 | -0.68 | -0.57 | -0.60 | -0.76 | -0.95 | -1.08 | -1.24 | -1.22 |
| 2022 | -1.18 | -1.13 | -1.19 | -1.26 | -1.17 | -0.95 | -0.86 | -0.94 | -1.07 | -1.10 | -1.05 | -0.99 |
| 2023 | -0.86 | -0.68 | -0.55 | -0.31 | -0.04 | 0.28 | 0.54 | 0.79 | 1.04 | 1.30 | 1.42 | 1.40 |
| 2024 | 1.13 | 0.78 | 0.42 | 0.04 | -0.28 | -0.45 | -0.49 | -0.56 | -0.67 | -0.76 | -0.88 | -1.05 |
| 2025 | -1.10 | -0.85 | -0.66 | -0.50 | -0.48 | -0.36 | -0.41 | -0.59 | -0.78 | -0.93 | -0.98 | -1.04 |
| 2026 | -0.91 | -0.76 | -0.44 | -0.04 | 0.49 | 0.97 | 1.36 | — | — | — | — | — |

*Table 3. This note's reconstruction, from ERSSTv6 and the Section 3
factor table alone.*

| Year | DJF | JFM | FMA | MAM | AMJ | MJJ | JJA | JAS | ASO | SON | OND | NDJ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2020 | 0.14 | 0.10 | 0.02 | -0.25 | -0.55 | -0.72 | -0.75 | -0.84 | -1.16 | -1.36 | -1.46 | -1.39 |
| 2021 | -1.32 | -1.18 | -1.08 | -0.89 | -0.68 | -0.57 | -0.60 | -0.76 | -0.94 | -1.08 | -1.24 | -1.22 |
| 2022 | -1.17 | -1.13 | -1.19 | -1.26 | -1.17 | -0.95 | -0.86 | -0.94 | -1.06 | -1.10 | -1.05 | -0.99 |
| 2023 | -0.86 | -0.68 | -0.55 | -0.31 | -0.04 | 0.28 | 0.55 | 0.79 | 1.05 | 1.30 | 1.43 | 1.41 |
| 2024 | 1.14 | 0.78 | 0.43 | 0.04 | -0.27 | -0.45 | -0.49 | -0.56 | -0.67 | -0.76 | -0.88 | -1.04 |
| 2025 | -1.09 | -0.84 | -0.65 | -0.50 | -0.47 | -0.36 | -0.40 | -0.58 | -0.78 | -0.93 | -0.97 | -1.03 |
| 2026 | -0.90 | -0.75 | -0.44 | -0.04 | 0.48 | 0.97 | 1.36 | — | — | — | — | — |

Over these 79 values the largest single discrepancy is 0.012°C (AMJ 2026:
published 0.49 vs. 0.478 reconstructed) and the RMS is 0.0047°C, consistent
with the pooled statistics above. At the one-decimal precision CPC's own
page displays, 75 of 79 values match exactly; the 4 exceptions are cases
where the published value sits almost exactly on a rounding boundary
(e.g. MAM 2020: published −0.25, reconstructed −0.253, which round to
different neighbors depending on the rounding rule applied), not genuine
divergence.

## 4. A seasonal-only fit, for comparison

CPC's RONI is not built this way (Section 2), but fitting a single factor
per season directly — ignoring the monthly structure and regressing the
3-month-averaged ERSSTv6 difference against the published seasonal
values — is a natural simpler alternative, and is instructive for
quantifying what the monthly structure buys. The seasonal RMS from this
simpler fit ranges from 0.0041°C (SON) to 0.0128°C (MJJ), averaging
0.0065°C across the 12 seasons — about a third higher than the
monthly-factor model's 0.0049°C, and more than double in the two
worst-fitting seasons, AMJ (0.0124°C) and MJJ (0.0128°C). These are
exactly the seasons straddling the sharpest part of the monthly factor's
annual cycle (April's peak of 1.390 falling to June's 1.248), consistent
with the serial-correlation mechanism described in Section 2: where the
monthly factor changes quickly within a season, treating it as constant
costs measurable accuracy.

## 5. Recipe summary

To reproduce RONI from ERSSTv6 alone:

1. **Box averages**: Niño-3.4 (5°S–5°N, 170°W–120°W) and tropical-mean
   (20°S–20°N, all longitudes) SST, cosine-latitude weighted.
2. **Monthly anomalies**: subtract a 1991–2020 monthly climatology from
   each box average.
3. **Unscaled monthly difference**: `diff = n34_anom − trop_anom`, at
   monthly resolution (no smoothing).
4. **Scale monthly**, using the fixed per-calendar-month factor from the
   table in Section 3: `relative_monthly = factor(month) * diff`.
5. **Seasonal RONI**: a 3-month centered running mean of
   `relative_monthly`.

Implementation: `scripts/roni_factor_backout.py`. Behavioral detail and
saved output tables: `specs/roni_factor_backout.md`,
`plots/roni_factor_backout/`.

## 6. Limitations

Even under the correct (monthly-factor) model, a residual of about
0.005°C RMS remains — larger than the roughly 0.003°C RMS expected from
rounding alone if the two series were rounded independently and matched
exactly otherwise — together with the small systematic negative bias noted
in Section 3. Two explanations are plausible but not distinguished by the
analysis here: a difference between the rounding convention used for the
published values and simple round-to-nearest, or a genuine mismatch
between the period (and possibly detrending) used in the calculation
above and whatever CPC actually uses for its own variance-matching step,
which is not documented on the RONI product page. Since the residual is
small relative to the factor's own annual cycle and well within the
0.01°C publication precision, it does not affect the practical
reproduction recipe in Section 5.

## References

- CPC Relative Oceanic Niño Index (RONI) product page:
  https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/
- L'Heureux, M. L., et al., 2024: A Relative Sea Surface Temperature Index
  for Classifying ENSO Events in a Changing Climate. *J. Climate*, **37**,
  1197–1211,
  [doi:10.1175/JCLI-D-23-0406.1](https://doi.org/10.1175/JCLI-D-23-0406.1).
- ERSSTv6 monthly SST:
  https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc
