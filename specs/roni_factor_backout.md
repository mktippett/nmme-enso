# roni_factor_backout.py — Behavioral Specification

> Last reviewed against code: 2026-09-10 (direct monthly fit against Rnino34.ascii.txt + monthly/seasonal self-consistency check added same session)

## Purpose

CPC's Relative Oceanic Niño Index (RONI) page
(https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/)
describes RONI's anomaly recipe (ERSSTv6, Niño-3.4 minus tropical-mean SST
anomaly, 1991-2020 base period, 3-month running mean) but not the numeric
scaling factor that adjusts the difference "so the variance equals the
original Niño 3.4 index," nor the period or detrending treatment behind
it. This script backs that factor out empirically from ERSSTv6 SST plus
CPC's own published series — the seasonal `RONI.ascii.txt` and the monthly
`Rnino34.ascii.txt` — in three stages: (1) a self-consistency check
between CPC's two published products, confirming the factor is applied
monthly, before any temporal averaging; (2) a per-calendar-month factor
fit directly against the published monthly series (the primary result);
(3) a simpler per-season factor fit, kept for comparison since it is *not*
CPC's actual method (established by stage 1). Exploratory/diagnostic
only: not part of the production pipeline, not referenced by any other
script. Full write-up for an external reader: `docs/CPC_RONI_scaling.md`.

## Inputs

| File | Relevant columns | Filters applied |
|------|-----------------|-----------------|
| `config.ERSSTV6_NC` (`observations/ERSSTv6.sst.mnmean.nc`) | `sst(time, lat, lon)` | Niño-3.4 box (`config.N34_LAT`/`N34_LON`, 5S-5N, 170W-120W) and tropical-mean box (`config.TROPICS_LAT`, 20S-20N, all lon), cosine-latitude weighted; anomalies vs. a 1991-2020 (`config.CLIM_START_YEAR`/`CLIM_END_YEAR`) monthly climatology |
| `config.RONI_TXT` (`observations/RONI.ascii.txt`, downloaded from https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt) | `SEAS, YR, ANOM` | none — full published record (DJF 1950 .. latest) |
| `config.RNINO34_TXT` (`observations/Rnino34.ascii.txt`, downloaded from https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt) | `YR, MTH, ANOM` (renamed `year, month, ANOM` on load) | none — full published record (Dec 1949 .. latest), continuous with no gaps |

## Outputs

All written to `config.PLOTS_DIR_RONI_FACTOR_BACKOUT` (`plots/roni_factor_backout/`):

| File | Contents | Format |
|------|----------|--------|
| `report.txt` | Full console report (all tables below, concatenated) | text |
| `monthly_seasonal_consistency.csv` | `check_monthly_seasonal_consistency()` output: per-(season, year), published seasonal RONI vs. 3-month rolling mean of published monthly RONI, and their difference | CSV |
| `monthly_factors_direct.csv` | `backout_monthly_factors_direct()` output: 12 calendar-month factors fit directly against `Rnino34.ascii.txt` + SE/95% CI/RMS — **the primary factor estimate** | CSV |
| `direct_monthly_vs_seasonal_rms.csv` | `validate_direct_monthly_against_seasonal()` output: per-season RMS when the direct monthly factors are applied monthly, 3-month averaged, and compared to the published *seasonal* RONI | CSV |
| `seasonal_factors_singlerow.csv` | `backout_factors()` output: per-season factor backed out from the single largest-\|diff\| year, validated against every other year (comparison only — not CPC's method) | CSV |
| `seasonal_factors_recent.csv` | `backout_factors_recent()` output: per-season OLS-through-origin factor + SE/95% CI from the 20 most recent years (comparison only — not CPC's method) | CSV |
| `monthly_factors.csv` | `backout_monthly_factors()` output: 12 calendar-month factors from the joint linear fit through seasonal aggregates only (superseded by `monthly_factors_direct.csv` now that the monthly target is available directly; kept for the pre/post comparison) + SE/95% CI | CSV |
| `seasonal_vs_monthly_rms.csv` | Per-season RMS residual under the seasonal-block fit vs. the (seasonal-aggregate-only) joint monthly fit | CSV |
| `factor_by_year.csv` | `factor_table_by_year()`: per-(year, season) `RONI / diff_3mo` ratio, pivoted year x season (CPC page layout) | CSV |

## Algorithm

### 1. Self-consistency check (logically first — no ERSSTv6 involved)

`check_monthly_seasonal_consistency()`: does a 3-month centered rolling mean of the published *monthly* series (`Rnino34.ascii.txt`) reproduce the published *seasonal* series (`RONI.ascii.txt`)? Labels the monthly series with CPC's season code (`_season_label`) and merges on `(SEAS, YR)`. This uses only CPC's two published products — no ERSSTv6 recomputation — and establishes whether CPC's factor is applied monthly (then averaged) or independently per season. See docs/CPC_RONI_scaling.md §2 for the result and its interpretation.

### 2. Load ERSSTv6 anomalies

`_ersstv6_monthly_anomalies()`: cosine-latitude-weighted box means (`_box_mean`, `weighted=True` by default) for the Niño-3.4 and tropical-mean boxes, monthly, anomalies vs. a 1991-2020 monthly climatology. `weighted=False` (plain grid-cell mean, no cos-lat weighting) is available as a parameter for testing only — kept as a fallback for reproducing CPC's number if it turns out they skip area weighting; **not the default** (confirmed 2026-09-10 not to matter — see Edge Cases; also matches the reference implementation's own choice, see Constants).

### 3. Two anomaly representations

- `load_ersstv6_3mo_anomalies()`: 3-month centered rolling mean of the monthly anomalies, labeled with CPC's 3-letter season code (`_season_label`, center-month convention) and year (`YR` = center month's own calendar year). Produces `diff_3mo = n34_3mo - trop_3mo`, the unscaled 3-month-averaged relative index.
- `monthly_diff_lookup()`: the same anomalies with **no** running mean, as a `{(year, month): diff}` dict — the per-month building block for the monthly fits (§4, §6).

### 4. Direct monthly factor fit (primary result)

`backout_monthly_factors_direct(n_recent=20)`: each calendar month's factor is fit independently by OLS-through-origin (`factor = Σxy/Σx²`) directly against `load_rnino34_table()`'s published monthly value, `x = diff_monthly` (from `monthly_diff_lookup()`), `y = published monthly RONI`, over that month's 20 most recent published years. No coupling across months is needed — every observation maps to exactly one calendar month's factor — since the monthly target is now available directly (§1 established that seasonal aggregation is just a 3-month average of this, so nothing is lost by fitting it this way rather than through seasonal aggregates as in §6). Standard errors as in §6's per-season formula, now per month.

`validate_direct_monthly_against_seasonal()`: applies the §4 factors monthly, takes a 3-month centered running mean of the *scaled* values, and compares against the published *seasonal* series (`load_merged()`'s `ANOM`) — the full round-trip check, independent confirmation that the monthly fit is also consistent with the seasonal product.

### 5. Match ERSSTv6 to published seasonal RONI

`load_merged()` inner-joins `load_ersstv6_3mo_anomalies()` to `load_roni_table()` (parsed `RONI.ascii.txt`) on `(SEAS, YR)`. Used by §6 and by the seasonal-block comparison fits (§7).

### 6. Joint monthly factor fit through seasonal aggregates only (superseded by §4)

`backout_monthly_factors(n_recent=20)`: solves for the same 12 monthly factors as §4, but *without* the monthly target — using only `load_merged()`'s seasonal aggregates, before `Rnino34.ascii.txt` was available to this analysis. Each published seasonal RONI value is treated as a linear combination of 3 (of 12 total) monthly factors:

```
RONI(season, year) = (1/3) * [ f(m-1)*diff(m-1) + f(m)*diff(m) + f(m+1)*diff(m+1) ]
```

where `m` is the season's center calendar month and `diff(k)` is `monthly_diff_lookup()`'s value at the appropriate (year, month) (with year rollover at Dec/Jan via `_neighbor_month`). One combined design matrix `X` (n_obs x 12, 3 nonzero entries per row) is built over all season-year rows in the most recent `n_recent` years; solved via `numpy.linalg.lstsq`. Standard errors: `sigma2 = SSR/(n_obs-12)`, `Cov(f) = sigma2 * (X'X)^-1`. Kept for the before/after comparison in `seasonal_vs_monthly_rms.csv`; its fitted factors and RMS are close to but not identical to §4's (direct fit against the true monthly target), confirming the seasonal-aggregate-only approach was already a good approximation once the correct monthly-then-average model was used — see Constants.

### 7. Seasonal-block factor (two estimators, for comparison only — not CPC's method)

- `backout_factors()`: per season, the reference row is the year with the largest `|diff_3mo|` (minimizes the relative effect of RONI's 2-decimal-place rounding); `factor_single_row = RONI_ref / diff_3mo_ref`, then validated by applying it to every *other* year of that season and reporting `max_abs_err`/`rms_err`. Also reports the naive mean/std of `RONI/diff_3mo` across all years of the season — included only to show this estimator is unstable (blows up when `diff_3mo` is near zero for a neutral year).
- `backout_factors_recent(n_recent=20)`: per season, OLS-through-origin (`factor = Σxy/Σx²`, x=`diff_3mo`, y=`RONI`) over that season's 20 most recent published years, with SE and 95% CI from `sigma² = Σresid²/(n-1)`, `SE = sqrt(sigma²/Σx²)`.

Both assume one constant factor per *season* — i.e. that the factor is uniform across a season's 3 months. §1 established this is not how CPC actually computes RONI; kept only to quantify the cost of that simplification (see Constants).

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `N34_LAT`/`N34_LON`, `TROPICS_LAT` | Reused from `config.py` | Guarantees this script's boxes match the rest of the project's relative-Niño-3.4 recipe; also matches CPC's own stated RONI box definitions |
| `BASE_START`/`BASE_END` | 1991-2020 (`config.CLIM_START_YEAR`/`CLIM_END_YEAR`) | Matches CPC's stated RONI base period |
| `n_recent=20` (§4, §6, §7) | 20 most recent published years per calendar month/season | Balances enough samples for a stable regression against staying representative of the current record; single-row backout (§7) confirmed no material drift vs. older reference years (1982-1997), so 20 years is not a sensitive choice |
| Reference row selection in `backout_factors()` (§7) | Largest \|diff_3mo\| | 2-decimal rounding in published RONI is a fixed absolute error; dividing by the largest-amplitude diff minimizes its relative effect. The naive alternative (mean of all years' ratios) is unstable — e.g. JAS: mean 2.72, std 13.1 — because near-neutral years have `diff_3mo` close to zero |
| Order of operations: monthly, then average (§1) | Confirmed, not assumed | `check_monthly_seasonal_consistency()`: all 919 published seasonal rows match a 3-month rolling mean of the published monthly series to within 0.0067 (exactly the two-independently-rounded-series bound) and a mean bias of 0.00007 — indistinguishable from zero, using only CPC's own two products (no ERSSTv6 involved). This is the basis for treating §4 (direct monthly fit) as authoritative over §7 (seasonal-block fit) |
| Seasonal-block factor (§7) is only an approximation | — | Confirmed wrong, not just noisier: the seasonal-block model's worst two seasons, AMJ (RMS 0.0124) and MJJ (0.0128), straddle the sharpest part of the monthly factor's annual cycle (April peak 1.390 falling to June 1.248) — treating the factor as constant across those 3 months measurably costs accuracy, consistent with §1's finding that a season's effective factor is a variance-weighted combination of 3 serially-correlated monthly factors, not their simple average |
| §4 (direct monthly fit) vs. §6 (joint fit through seasonal aggregates only) | Close but not identical | E.g. Apr: 1.390 (§4) vs. 1.387 (§6); Aug: 1.167 (§4) vs. 1.175 (§6) — within ~0.01 everywhere. Confirms §6's seasonal-aggregate-only approach (used before `Rnino34.ascii.txt` was pulled into this analysis) was already a good approximation, once the correct monthly-then-average model was used instead of §7's per-season blocks |
| `weighted=False` tested, not default | Cos-lat area weighting vs. plain grid-cell mean | Tested 2026-09-10: unweighted gave *no* improvement (joint-fit RMS 0.004884 vs. 0.004819 weighted) — ruled out as the source of the residual discussed below. The reference implementation this recipe is modeled on (M. L'Heureux's public notebook, github.com/michellelheureux/Relative-SST) also uses a plain, unweighted grid-cell mean for both boxes, consistent with this being immaterial rather than a discrepancy. Left in as a `weighted` parameter (default `True`) rather than removed |

## Edge Cases & Error Handling

- **`YR` label for Dec-centered seasons**: CPC labels each season by its center month's own calendar year (e.g. NDJ 1997 = Nov97,Dec97,Jan98, center month Dec97 → "1997"; DJF 1998 = Dec97,Jan98,Feb98, center month Jan98 → "1998") — verified against `RONI.ascii.txt`'s 1997-98 El Niño rows. `load_ersstv6_3mo_anomalies()` sets `YR = year` directly from the center-month timestamp, which already has this convention baked in — no adjustment needed.
- **Near-zero `diff_3mo` years blow up simple ratios**: `backout_factors()`'s `factor_mean_all_years` diagnostic and `factor_table_by_year()`'s per-year grid both show this directly (e.g. 1967 JAS ratio = 116.51, 1986 MJJ = -4.44) — expected and documented, not a bug. The single-row (largest-\|diff\|) and OLS-through-origin estimators are unaffected since both naturally downweight near-zero-`diff` observations.
- **Unresolved residual structure, ruled out to not be a seasonal-aggregation artifact**: even the direct monthly fit (§4, the most direct possible test — fit and target are both monthly, no seasonal aggregation involved at all) does not reach zero error: pooled RMS 0.0054 across 240 fitted points, with a small systematic negative bias (mean -0.0038) and a negative correlation with the published anomaly itself (corr ≈ -0.28 to -0.32 with ANOM, i.e. the model tends to overshoot the magnitude of large excursions in both directions). Since §4 matches §6's seasonal-aggregate-only fit closely (see Constants) and §1 confirms the two published CPC products are mutually consistent to near the rounding floor, this residual is not an artifact of season-vs-month aggregation, and not explained by area-weighting (tested, ruled out — see Constants) either. See Open Items for remaining hypotheses.
- **2026 partial coverage**: `RONI.ascii.txt`/`Rnino34.ascii.txt` run through JJA/Aug 2026 at time of writing, so DJF-JJA (or Jan-Aug) use published years 2007-2026 in the 20-year windows while JAS-NDJ (or Sep-Dec) use 2006-2025 (one year lagged) — not a script bug, just the data's cutoff (`.tail(n_recent)`/year-cutoff logic handles this per month/season automatically).

## Open Items

- **Residual source unresolved** (see Edge Cases): truncation-vs-rounding convention and amplitude-dependent nonlinearity are both live, untested hypotheses for the ~0.005 RMS residual remaining after the direct monthly fit (§4) — now confirmed not attributable to seasonal aggregation or area weighting. Next diagnostic step discussed but not run: check whether residual sign flips with published-value sign the way a truncate-toward-zero convention would predict.
- **`observations/RONI.ascii.txt` and `Rnino34.ascii.txt` are live CPC products**: unlike ERSSTv5/v6 (refreshed via the README Quickstart's `curl -z`), there is no refresh command for either file yet. Re-download manually from https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt / https://www.cpc.ncep.noaa.gov/data/indices/Rnino34.ascii.txt if this script needs to be rerun against a newer record.

## Verification Snippet

```python
# Run after changes to confirm key invariants
import sys
sys.path.insert(0, "scripts")
from roni_factor_backout import (
    check_monthly_seasonal_consistency, load_merged, backout_factors,
    backout_factors_recent, backout_monthly_factors,
    backout_monthly_factors_direct, validate_direct_monthly_against_seasonal,
)
import config

consistency = check_monthly_seasonal_consistency()
assert len(consistency) > 900, f"expected ~919 matched rows, got {len(consistency)}"
assert consistency["diff"].abs().max() < 0.01, "monthly-vs-seasonal mismatch exceeds one rounding step"

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

monthly_direct, rnino34_matched = backout_monthly_factors_direct(n_recent=20)
assert len(monthly_direct) == 12
assert 1.0 < monthly_direct.factor.min() and monthly_direct.factor.max() < 1.5
seasonal_check, seasonal_rms = validate_direct_monthly_against_seasonal(monthly_direct, merged, n_recent=20)
assert seasonal_rms < 0.01

for f in ["report.txt", "monthly_seasonal_consistency.csv", "monthly_factors_direct.csv",
          "direct_monthly_vs_seasonal_rms.csv", "seasonal_factors_singlerow.csv",
          "seasonal_factors_recent.csv", "monthly_factors.csv", "seasonal_vs_monthly_rms.csv",
          "factor_by_year.csv"]:
    assert (config.PLOTS_DIR_RONI_FACTOR_BACKOUT / f).exists(), f"missing output {f}"

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-09-10 | Initial `scripts/roni_factor_backout.py` written: seasonal-block factor backout (single-row + OLS-through-origin variants), then superseded/supplemented by a joint 12-calendar-month linear fit (per user correction that CPC scales monthly values before seasonal averaging) — reduced worst-season RMS (AMJ/MJJ) from 0.0124-0.0128 to 0.0044-0.0056. Tested and ruled out area-weighting as the source of the small remaining residual. Added `config.ERSSTV6_NC`, `config.RONI_TXT`, `config.PLOTS_DIR_RONI_FACTOR_BACKOUT`. Script now saves its report + all tables to `plots/roni_factor_backout/`. | ✓ |
| 2026-09-10 | Same-day follow-up, three additions: (1) `check_monthly_seasonal_consistency()` — confirmed CPC's published monthly and seasonal products are mutually consistent (max diff 0.0067, mean 0.00007) using no ERSSTv6 data at all, establishing the monthly-then-average order of operations directly rather than assuming it. (2) `backout_monthly_factors_direct()` + `validate_direct_monthly_against_seasonal()` — fit each month's factor directly against the newly-added `Rnino34.ascii.txt` (published monthly series), now the primary factor estimate, superseding the seasonal-aggregate-only joint fit; confirmed close agreement with it (within ~0.01 per month) and confirmed the residual is not a seasonal-aggregation artifact. (3) Cross-checked the recipe (box means, climatology, order of operations) against M. L'Heureux's public reference notebook (github.com/michellelheureux/Relative-SST/blob/main/compute_relativeNino34_ONI_JClimate.ipynb, harvested for reference only, not incorporated into the codebase) — confirmed unweighted box means and the monthly-then-average order, no changes to the recipe resulted. Added `config.RNINO34_TXT`. Wrote `docs/CPC_RONI_scaling.md` (external-reader write-up of all three checks). | ✓ |
