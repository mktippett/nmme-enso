# roni_factor_backout.py — Behavioral Specification

> Last reviewed against code: 2026-09-10 (initial spec, monthly-factor fit + saved output added same session)

## Purpose

CPC's Relative Oceanic Niño Index (RONI) page
(https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/)
describes RONI's anomaly recipe (ERSSTv6, Niño-3.4 minus tropical-mean SST
anomaly, 1991-2020 base period, 3-month running mean) but not the numeric
scaling factor that adjusts the difference "so the variance equals the
original Niño 3.4 index." This script backs that factor out empirically
from ERSSTv6 SST plus CPC's own published `RONI.ascii.txt` table, first at
seasonal granularity (12 season-block fits) and then — per the user's
2026-09-10 correction that CPC scales *monthly* relative-index values and
only then averages 3 months into a season — via one joint linear fit over
12 shared calendar-month factors. Exploratory/diagnostic only: not part of
the production pipeline, not referenced by any other script.

## Inputs

| File | Relevant columns | Filters applied |
|------|-----------------|-----------------|
| `config.ERSSTV6_NC` (`observations/ERSSTv6.sst.mnmean.nc`) | `sst(time, lat, lon)` | Niño-3.4 box (`config.N34_LAT`/`N34_LON`, 5S-5N, 170W-120W) and tropical-mean box (`config.TROPICS_LAT`, 20S-20N, all lon), cosine-latitude weighted; anomalies vs. a 1991-2020 (`config.CLIM_START_YEAR`/`CLIM_END_YEAR`) monthly climatology |
| `config.RONI_TXT` (`observations/RONI.ascii.txt`, downloaded from https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt) | `SEAS, YR, ANOM` | none — full published record (DJF 1950 .. latest) |

## Outputs

All written to `config.PLOTS_DIR_RONI_FACTOR_BACKOUT` (`plots/roni_factor_backout/`):

| File | Contents | Format |
|------|----------|--------|
| `report.txt` | Full console report (all tables below, concatenated) | text |
| `seasonal_factors_singlerow.csv` | `backout_factors()` output: per-season factor backed out from the single largest-\|diff\| year, validated against every other year | CSV |
| `seasonal_factors_recent.csv` | `backout_factors_recent()` output: per-season OLS-through-origin factor + SE/95% CI from the 20 most recent years | CSV |
| `monthly_factors.csv` | `backout_monthly_factors()` output: 12 calendar-month factors (joint fit) + SE/95% CI | CSV |
| `seasonal_vs_monthly_rms.csv` | Per-season RMS residual under the seasonal-block fit vs. the joint monthly fit | CSV |
| `factor_by_year.csv` | `factor_table_by_year()`: per-(year, season) `RONI / diff_3mo` ratio, pivoted year x season (CPC page layout) | CSV |

## Algorithm

### 1. Load ERSSTv6 anomalies

`_ersstv6_monthly_anomalies()`: cosine-latitude-weighted box means (`_box_mean`, `weighted=True` by default) for the Niño-3.4 and tropical-mean boxes, monthly, anomalies vs. a 1991-2020 monthly climatology. `weighted=False` (plain grid-cell mean, no cos-lat weighting) is available as a parameter for testing only — kept as a fallback for reproducing CPC's number if it turns out they skip area weighting; **not the default** (confirmed 2026-09-10 not to matter — see Edge Cases).

### 2. Two anomaly representations

- `load_ersstv6_3mo_anomalies()`: 3-month centered rolling mean of the monthly anomalies, labeled with CPC's 3-letter season code (`_season_label`, center-month convention) and year (`YR` = center month's own calendar year). Produces `diff_3mo = n34_3mo - trop_3mo`, the unscaled 3-month-averaged relative index.
- `monthly_diff_lookup()`: the same anomalies with **no** running mean, as a `{(year, month): diff}` dict — the per-month building block for the joint monthly fit.

### 3. Match to published RONI

`load_merged()` inner-joins `load_ersstv6_3mo_anomalies()` to `load_roni_table()` (parsed `RONI.ascii.txt`) on `(SEAS, YR)`.

### 4. Seasonal-block factor (two estimators, superseded by §5 for the primary answer)

- `backout_factors()`: per season, the reference row is the year with the largest `|diff_3mo|` (minimizes the relative effect of RONI's 2-decimal-place rounding); `factor_single_row = RONI_ref / diff_3mo_ref`, then validated by applying it to every *other* year of that season and reporting `max_abs_err`/`rms_err`. Also reports the naive mean/std of `RONI/diff_3mo` across all years of the season — included only to show this estimator is unstable (blows up when `diff_3mo` is near zero for a neutral year).
- `backout_factors_recent(n_recent=20)`: per season, OLS-through-origin (`factor = Σxy/Σx²`, x=`diff_3mo`, y=`RONI`) over that season's 20 most recent published years, with SE and 95% CI from `sigma² = Σresid²/(n-1)`, `SE = sqrt(sigma²/Σx²)`.

Both assume one constant factor per *season* — i.e. that the factor is uniform across a season's 3 months. Superseded as the best estimate by §5, which found this assumption measurably wrong for AMJ/MJJ (see Constants).

### 5. Joint monthly factor fit (primary result)

`backout_monthly_factors(n_recent=20)`: CPC scales the *monthly* relative-index value first, `factor(month) * diff_monthly(month)`, and only then averages 3 consecutive months into a season. So each published seasonal RONI value is a linear combination of 3 (of 12 total) monthly factors:

```
RONI(season, year) = (1/3) * [ f(m-1)*diff(m-1) + f(m)*diff(m) + f(m+1)*diff(m+1) ]
```

where `m` is the season's center calendar month and `diff(k)` is `monthly_diff_lookup()`'s value at the appropriate (year, month) (with year rollover at Dec/Jan via `_neighbor_month`). One combined design matrix `X` (n_obs x 12, 3 nonzero entries per row) is built over all season-year rows in the most recent `n_recent` years; solved via `numpy.linalg.lstsq`. Standard errors: `sigma2 = SSR/(n_obs-12)`, `Cov(f) = sigma2 * (X'X)^-1`.

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `N34_LAT`/`N34_LON`, `TROPICS_LAT` | Reused from `config.py` | Guarantees this script's boxes match the rest of the project's relative-Niño-3.4 recipe; also matches CPC's own stated RONI box definitions |
| `BASE_START`/`BASE_END` | 1991-2020 (`config.CLIM_START_YEAR`/`CLIM_END_YEAR`) | Matches CPC's stated RONI base period |
| `n_recent=20` (both `backout_factors_recent` and `backout_monthly_factors`) | 20 most recent published years per season | Balances enough samples for a stable regression against staying representative of the current record; single-row backout (§4) confirmed no material drift vs. older reference years (1982-1997), so 20 years is not a sensitive choice |
| Reference row selection in `backout_factors()` | Largest \|diff_3mo\| | 2-decimal rounding in published RONI is a fixed absolute error; dividing by the largest-amplitude diff minimizes its relative effect. The naive alternative (mean of all years' ratios) is unstable — e.g. JAS: mean 2.72, std 13.1 — because near-neutral years have `diff_3mo` close to zero |
| Seasonal-block factor is only an approximation | — | Joint monthly fit (§5) reduced overall RMS residual from a per-season range of 0.0041-0.0128 to a uniform ~0.0043-0.0056; the two worst-fitting seasons under the block model, AMJ (0.0124→0.0056) and MJJ (0.0128→0.0044), straddle the sharpest part of the factor's annual cycle (April peak 1.3868 → June 1.2503), confirming the block model's single-factor-per-season assumption was measurably wrong there, not just noisier data |
| `weighted=False` tested, not default | Cos-lat area weighting vs. plain grid-cell mean | Tested 2026-09-10 per user request ("scientifically wrong but we are trying to replicate [CPC]"): unweighted gave *no* improvement (joint-fit RMS 0.004884 vs. 0.004819 weighted) — ruled out as the source of the residual. Left in as a `weighted` parameter (default `True`) rather than removed, in case it's revisited |

## Edge Cases & Error Handling

- **`YR` label for Dec-centered seasons**: CPC labels each season by its center month's own calendar year (e.g. NDJ 1997 = Nov97,Dec97,Jan98, center month Dec97 → "1997"; DJF 1998 = Dec97,Jan98,Feb98, center month Jan98 → "1998") — verified against `RONI.ascii.txt`'s 1997-98 El Niño rows. `load_ersstv6_3mo_anomalies()` sets `YR = year` directly from the center-month timestamp, which already has this convention baked in — no adjustment needed.
- **Near-zero `diff_3mo` years blow up simple ratios**: `backout_factors()`'s `factor_mean_all_years` diagnostic and `factor_table_by_year()`'s per-year grid both show this directly (e.g. 1967 JAS ratio = 116.51, 1986 MJJ = -4.44) — expected and documented, not a bug. The single-row (largest-\|diff\|) and OLS-through-origin estimators are unaffected since both naturally downweight near-zero-`diff` observations.
- **Unresolved residual structure**: even the joint monthly fit (§5) does not reach zero error (overall RMS 0.0048, dof=223) — above the pure-rounding floor (`0.01/sqrt(12)` ≈ 0.0029). Residuals show a small systematic negative mean in every season (-0.0026 to -0.0042) and a negative correlation with the published anomaly itself (corr ≈ -0.32 with ANOM), i.e. the model tends to overshoot the magnitude of large excursions in both directions. Neither area-weighting (tested, ruled out — see Constants) nor the monthly-vs-seasonal factor order (tested, substantially improved but did not eliminate the residual) fully explains this. Live hypotheses, not yet tested: CPC truncates published RONI toward zero rather than rounding to nearest (would produce exactly this sign pattern given the recent sample's skew toward stronger positive ENSO states); a genuinely amplitude-dependent (nonlinear) factor; or a mismatch between this project's ERSSTv6 file/processing and CPC's internal one. See Open Items.
- **2026 partial coverage**: `RONI.ascii.txt` runs through JJA 2026 at time of writing, so DJF-JJA use published years 2007-2026 in the 20-year windows while JAS-NDJ use 2006-2025 (one season lagged) — not a script bug, just the data's cutoff (`backout_factors_recent`/`backout_monthly_factors`'s `.tail(n_recent)`/year-cutoff logic handles this per season automatically).

## Open Items

- **Residual source unresolved** (see Edge Cases): truncation-vs-rounding convention and amplitude-dependent nonlinearity are both live, untested hypotheses for the ~0.0048 RMS residual remaining after the joint monthly fit. Next diagnostic step discussed but not run: check whether residual sign flips with published-value sign the way a truncate-toward-zero convention would predict.
- **`observations/RONI.ascii.txt` is a live CPC product**: unlike ERSSTv5/v6 (refreshed via the README Quickstart's `curl -z`), there is no refresh command for this file yet. Re-download manually from https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt if this script needs to be rerun against a newer record.

## Verification Snippet

```python
# Run after changes to confirm key invariants
import sys
sys.path.insert(0, "scripts")
from roni_factor_backout import load_merged, backout_factors, backout_factors_recent, backout_monthly_factors
import config

merged = load_merged()
assert len(merged) > 900, f"expected ~919 matched rows, got {len(merged)}"

single = backout_factors(merged=merged)
assert len(single) == 12
assert single.rms_err_other_years.max() < 0.05, "single-row backout should validate tightly"

recent = backout_factors_recent(n_recent=20, merged=merged)
assert len(recent) == 12
assert (recent.se > 0).all()

monthly, rms, n_obs = backout_monthly_factors(n_recent=20, merged=merged)
assert len(monthly) == 12
assert 1.0 < monthly.factor.min() and monthly.factor.max() < 1.5, "factors should be O(1), not blown up"
assert rms < 0.01

for f in ["report.txt", "seasonal_factors_singlerow.csv", "seasonal_factors_recent.csv",
          "monthly_factors.csv", "seasonal_vs_monthly_rms.csv", "factor_by_year.csv"]:
    assert (config.PLOTS_DIR_RONI_FACTOR_BACKOUT / f).exists(), f"missing output {f}"

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-09-10 | Initial `scripts/roni_factor_backout.py` written: seasonal-block factor backout (single-row + OLS-through-origin variants), then superseded/supplemented by a joint 12-calendar-month linear fit (per user correction that CPC scales monthly values before seasonal averaging) — reduced worst-season RMS (AMJ/MJJ) from 0.0124-0.0128 to 0.0044-0.0056. Tested and ruled out area-weighting as the source of the small remaining residual. Added `config.ERSSTV6_NC`, `config.RONI_TXT`, `config.PLOTS_DIR_RONI_FACTOR_BACKOUT`. Script now saves its report + all tables to `plots/roni_factor_backout/`. | ✓ |
