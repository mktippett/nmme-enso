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
almost no difference to how well the recipe reproduces the published
series: across the 12 calendar months the two choices differ by at most
0.0002 in RMS error. They do differ in the fitted factor itself, by 0.0023
to 0.0034 — larger than that factor's own standard error, and in the same
direction for every month. The factors quoted below are the
cosine-weighted ones; each would fall by about 0.003 if computed from a
plain grid-cell mean, which is how the reference implementation listed
under References computes its box averages. That shift is well inside the
0.01 precision at which RONI is published, so it is not a source of
uncertainty in the reproduction, but it is not negligible in the factor.

CPC's description amounts to the following formula, at monthly resolution:

$$
\text{RONI}(\text{month}) = f(\text{month}) \times \Big(\underbrace{\text{N34}_{\text{anom}}}_{\text{known}} - \underbrace{\text{Trop}_{\text{anom}}}_{\text{known}}\Big)
$$

The two box anomalies above, $\text{N34}_{\text{anom}}$ and
$\text{Trop}_{\text{anom}}$, are both computable directly from ERSSTv6. The
per-calendar-month factor $f(\text{month})$ is the one unknown quantity:
CPC states only that it is chosen so that "the variance equals the
original Niño 3.4 index," not its numeric value, the period over which
that variance is computed, or whether the underlying series is
detrended first. Backing out $f(\text{month})$ from the known left- and
right-hand sides is the task of Sections 3–5.

## 2. The published monthly and seasonal series are mutually consistent

Before fitting anything, the two published series were checked against
each other directly: the 3-month centered running mean of the published
*monthly* series was compared with the published *seasonal* series, with
no involvement of ERSSTv6 or any fitted factor. All 919 published seasonal
values (1950–2026) matched a corresponding monthly-average value. The
maximum discrepancy is 0.0067°C and the mean discrepancy is 0.00007°C,
indistinguishable from zero.

That maximum is exactly the bound expected from two series that are each
independently rounded to two decimal places before publication, and the
argument takes two steps. A published seasonal value is a multiple of
0.01, and the mean of three published monthly values is a multiple of
0.01/3, so the discrepancy between them is necessarily a multiple of
0.0033. It is also strictly smaller than one full 0.01 step, being the
difference between the seasonal value's own rounding error and the average
of the three monthly ones, both drawn from the same 0.01-wide window. The
largest attainable discrepancy is therefore two-thirds of a step, 0.0067,
which is what is observed. The observed RMS of 0.0032 likewise matches the
0.0033 predicted for two independently rounded series. No unmodeled
processing separates the two products.

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
the published values). That bias is systematic rather than sampling noise,
and Section 5 identifies where it comes from. Applying the same monthly
factors and then taking a
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

The disagreements are also one-sided, which is the clue Section 5 follows.
Of the 79 values, 24 are reconstructed more than 0.005°C *above* the
published value and only 1 more than 0.005°C below.

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

## 5. The published values are truncated, not rounded

The residual left by Section 3 is not scatter. It is a bias of one
consistent sign and roughly constant size, and rounding to nearest cannot
produce one: rounding is symmetric, so it adds noise with zero mean.
Truncating downward is not symmetric. It subtracts, on average, exactly
half of the 0.01 publication step — 0.005°C, which is the size of the
observed bias.

Two tests confirm this. The first refits the same model with one intercept
shared by all 12 calendar months, `RONI = f(month) × diff + c`, so that any
constant offset is separated from the factors instead of being absorbed
into them. The fitted intercept is −0.0042°C on the same 20-year sample
Section 3 uses, and −0.0051°C over the full 1950–2026 monthly record. The
residual scatter around it falls to 0.0037°C and 0.0033°C respectively,
against a quantization floor of 0.01/√12 = 0.0029°C for values published
to two decimals. The factors themselves barely move:

| Month | Through origin (Section 3) | Refit with intercept | Change |
|---|---|---|---|
| Jan | 1.2476 | 1.2461 | −0.0015 |
| Feb | 1.2737 | 1.2717 | −0.0019 |
| Mar | 1.3391 | 1.3365 | −0.0026 |
| Apr | 1.3900 | 1.3859 | −0.0041 |
| May | 1.3563 | 1.3516 | −0.0047 |
| Jun | 1.2479 | 1.2456 | −0.0023 |
| Jul | 1.1870 | 1.1858 | −0.0013 |
| Aug | 1.1671 | 1.1658 | −0.0014 |
| Sep | 1.2022 | 1.2005 | −0.0018 |
| Oct | 1.2149 | 1.2136 | −0.0014 |
| Nov | 1.2205 | 1.2196 | −0.0009 |
| Dec | 1.2322 | 1.2313 | −0.0009 |

The second test applies each candidate convention to the unrounded
prediction and counts how often it reproduces the published two-decimal
value *exactly*. Over all 920 published monthly values, using the refit
factors:

| Convention | Exact reproductions | RMS |
|---|---|---|
| Floor (round down) | 88.8% | 0.0035 |
| Round to nearest | 48.7% | 0.0072 |
| Truncate toward zero | 45.3% | 0.0079 |

The published *seasonal* series behaves the same way — flooring the
3-month mean of the unrounded scaled monthly values reproduces 89.3% of
the 919 published values exactly, against 50.1% for rounding to nearest.
Note that the seasonal value is floored from the unrounded monthly values,
not assembled from the already-published two-decimal ones.

The remaining 11% is what our own inputs cost. Our Niño-3.4 minus
tropical-mean difference is not identical to CPC's; the ~0.0015°C by which
the residual exceeds the quantization floor is enough to push a prediction
that lands near a 0.01 boundary onto the wrong side of it, which accounts
for roughly one value in nine.

One caveat. Flooring with no offset, and rounding to nearest applied to a
series that genuinely sits 0.005°C below ours, are observationally
identical — both put the residual in the same place. What favours flooring
is that the fitted offset lands at *exactly* half a publication step
rather than at some arbitrary value, which would be a coincidence under
the alternative. This also settles the other hypothesis raised below: a
mismatch in the period used for CPC's variance matching would rescale the
factor, and cannot introduce a constant offset.

## 6. Recipe summary

To reproduce RONI from ERSSTv6 alone:

1. **Box averages**: Niño-3.4 (5°S–5°N, 170°W–120°W) and tropical-mean
   (20°S–20°N, all longitudes) SST, cosine-latitude weighted.
2. **Monthly anomalies**: subtract a 1991–2020 monthly climatology from
   each box average.
3. **Unscaled monthly difference**: `diff = n34_anom − trop_anom`, at
   monthly resolution (no smoothing).
4. **Scale monthly**, using the fixed per-calendar-month factor from the
   refit column of the table in Section 5:
   `relative_monthly = factor(month) * diff`. The Section 3 factors absorb
   part of the offset that step 6 handles, and reproduce 83.5% of the
   published monthly values exactly rather than 88.8%.
5. **Seasonal RONI**: a 3-month centered running mean of
   `relative_monthly`.
6. **Publication convention**: to match the published two-decimal values,
   round *down* to two decimals, at both monthly and seasonal resolution
   (Section 5). Skipping this step leaves every value about 0.005°C high.

Implementation: `scripts/roni_factor_backout.py`. Behavioral detail and
saved output tables: `specs/roni_factor_backout.md`,
`plots/roni_factor_backout/`.

## 7. Limitations

With the publication convention accounted for (Section 5), the residual
falls to 0.0033°C RMS, against the 0.0029°C floor imposed by two-decimal
publication. What is left is close enough to that floor that it carries
little information, but three things remain unresolved.

The flooring interpretation is not proven, only strongly favoured, for the
reason given at the end of Section 5. The mean residual also drifts slowly
across the record, from −0.0065°C in the 1950s to −0.0041°C in the 2020s,
where a pure publication convention would hold it at −0.005°C throughout;
a slowly varying difference of order 0.002°C between our inputs and CPC's,
such as a dataset revision, would produce that. And the period, and any
detrending, behind CPC's variance-matching step remains undocumented —
Section 5 rules it out as the source of the *bias*, but not as a reason
the factors themselves might be revised in future.

None of this affects the reproduction recipe in Section 6. Of the 920
published monthly values it reproduces 817 exactly and all but one of the
rest to within one unit in the last published digit; the seasonal series
behaves the same way, 821 of 919 exact. The single exception is June 2026,
published as 1.07 against 1.047 reconstructed, which also drags the AMJ
2026 season out by two units. Both are at the very end of the record,
where ERSSTv6 is still subject to revision after CPC computes its
published value.

## References

- CPC Relative Oceanic Niño Index (RONI) product page:
  https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/
- L'Heureux, M. L., et al., 2024: A Relative Sea Surface Temperature Index
  for Classifying ENSO Events in a Changing Climate. *J. Climate*, **37**,
  1197–1211,
  [doi:10.1175/JCLI-D-23-0406.1](https://doi.org/10.1175/JCLI-D-23-0406.1).
- Reference implementation of the relative Niño-3.4 recipe, M. L'Heureux:
  https://github.com/michellelheureux/Relative-SST
- ERSSTv6 monthly SST:
  https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc
