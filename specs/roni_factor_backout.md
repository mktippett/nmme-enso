# roni_factor_backout.py — Behavioral Specification

> Last reviewed against code: 2026-09-10 (correctness review of `docs/CPC_RONI_scaling.md` and this script; publication-convention diagnostics added, residual resolved)

## Purpose

CPC's Relative Oceanic Niño Index (RONI) page
(https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/)
describes RONI's anomaly recipe (ERSSTv6, Niño-3.4 minus tropical-mean SST
anomaly, 1991-2020 base period, 3-month running mean) but not the numeric
scaling factor that adjusts the difference "so the variance equals the
original Niño 3.4 index," nor the period or detrending treatment behind
it. This script backs that factor out empirically from ERSSTv6 SST plus
CPC's own published series — the seasonal `RONI.ascii.txt` and the monthly
`Rnino34.ascii.txt` — in four stages: (1) a self-consistency check
between CPC's two published products, confirming the factor is applied
monthly, before any temporal averaging; (2) a per-calendar-month factor
fit directly against the published monthly series (the primary result);
(3) a simpler per-season factor fit, kept for comparison since it is *not*
CPC's actual method (established by stage 1); (4) publication-convention
diagnostics that identify the small residual left by stage 2 as CPC
flooring its published values to two decimals, plus the forward
reconstruction that closes the loop from ERSSTv6 back to the published
tables. Exploratory/diagnostic only: not part of the production pipeline,
not referenced by any other script. Full write-up for an external reader:
`docs/CPC_RONI_scaling.md`.

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
| `rounding_convention.csv` | `rounding_convention_check()` output, monthly series: per (sample, factor model) raw residual statistics, the fitted shared intercept, and the exact-reproduction rate + RMS + max error of each candidate 2-decimal convention (floor / round-to-nearest / truncate-toward-zero) | CSV |
| `rounding_convention_seasonal.csv` | `seasonal_convention_check()` output: the same convention comparison against the published *seasonal* series | CSV |
| `monthly_factors_offset_refit.csv` | The 12 factors refit with a free shared intercept, beside the §4 through-origin values and their difference | CSV |
| `residual_by_decade.csv` | `residual_by_decade()`: mean and RMS residual per decade under the refit factors — the drift a pure publication convention does not explain | CSV |
| `weighting_sensitivity.csv` | `weighting_sensitivity()`: §4 factors under cos-lat weighting vs. a plain grid-cell mean, with each fit's RMS | CSV |
| `seasonal_published_recent.csv`, `seasonal_reconstructed_recent.csv`, `seasonal_reconstruction_stats.csv` | `reconstruction_tables()`: published and reconstructed seasonal RONI in CPC's year x season layout for 2020-present, plus n / max error / RMS / 1-decimal agreement count (Tables 2-3 of `docs/CPC_RONI_scaling.md`) | CSV |

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

### 8. Publication convention (identifies the §4 residual)

`rounding_convention_check(n_recent=20)`: the §4 fit leaves a systematic *negative* bias, which round-to-nearest cannot produce (it is symmetric) but flooring can (it subtracts half a step on average). Two diagnostics, each run on the §4 20-year sample and on the full 1950-present record, both built from `monthly_fit_frame()` (published monthly value joined to `monthly_diff_lookup()`'s diff, the shared `(x, y)` table for every monthly-resolution fit):

- `_fit_slopes_with_offset()`: refit `y = f(month)*diff + c` with one intercept `c` shared by all 12 months, so a constant offset is separated from the factors instead of absorbed into them. Design matrix is 13 columns (12 per-month slopes + intercept), solved by `numpy.linalg.lstsq`.
- `_quantize_match()`: apply floor, round-to-nearest and truncate-toward-zero to the unrounded prediction and count exact reproductions of the published 2-decimal value, with each convention's RMS and max error.

`seasonal_convention_check(factors)`: the same convention comparison at seasonal resolution, via `reconstruct_from_factors()`. Confirms the published seasonal value is floored from the *unrounded* scaled monthly values, not assembled from the already-published 2-decimal monthly ones.

`residual_by_decade(factors)`: mean/RMS residual per decade. A pure publication convention would hold the mean at half a step in every decade; it drifts instead (see Edge Cases).

### 9. Forward reconstruction (the round trip)

`reconstruct_from_factors(factors)`: scale each ERSSTv6 monthly diff by its calendar month's factor, then take a 3-month centered running mean of the *scaled* values, labeled with CPC's season code — the recipe of `docs/CPC_RONI_scaling.md` §6 steps 4-5, in one function.

`reconstruction_tables(factors, start_year=2020)`: published vs. reconstructed seasonal RONI as two pivot tables in CPC's own year x season layout, with n / max absolute error / RMS / count agreeing at the 1-decimal precision CPC's product page displays / counts of reconstructions more than 0.005 above and below published, plus the rows that disagree at 1 decimal. Reproduces Tables 2-3 of the write-up.

### 10. Weighting sensitivity

`weighting_sensitivity(n_recent=20)`: runs §4 twice, `weighted=True` and `weighted=False`, and reports both factor sets, their difference, and each fit's RMS. Quantifies the claim in `docs/CPC_RONI_scaling.md` §1 (see Constants).

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
| `weighted=False` tested, not default | Cos-lat area weighting vs. plain grid-cell mean | Immaterial for *fit quality*, not for the *factor*. Measured 2026-09-10 by `weighting_sensitivity()` (§10): the two choices differ by at most 0.0002 in any month's RMS, but by 0.0023-0.0034 in the fitted factor — larger than that factor's own SE, and negative for all 12 months. The reference implementation this recipe is modeled on (M. L'Heureux's public notebook, github.com/michellelheureux/Relative-SST) uses a plain, unweighted grid-cell mean, so the published factors likely run ~0.003 above CPC's own. Still well inside RONI's 0.01 publication precision. Default stays `True` for consistency with the rest of the project's box means; the `weighted` parameter is kept, not removed. An earlier version of this row claimed the difference was under 0.0001 — that figure was the *RMS* difference, not the factor difference |
| Publication convention: floor, not round (§8) | Confirmed 2026-09-10 | Refitting with a free shared intercept gives −0.0042 on the 20-year sample and −0.0051 over the full record — half of the 0.01 publication step — and drops residual scatter to 0.0037/0.0033 against a 0.0029 quantization floor, while moving the factors by at most 0.005. Applying each convention to the unrounded prediction, flooring reproduces 88.8% of the 920 published monthly values exactly vs. 48.7% for round-to-nearest and 45.3% for truncate-toward-zero; seasonally, 89.3% vs. 50.1%. Not proof — flooring with no offset and round-to-nearest with a genuine −0.005 offset are observationally identical — but the offset landing at *exactly* half a step is a coincidence under the alternative. Also rules out a variance-matching-period mismatch as the source of the bias: that would rescale the factor, not add a constant |
| Two factor tables (§4 through-origin vs. §8 refit) | Both kept | The §4 factors absorb the flooring offset into the slope; the §8 refit separates them and is the better estimate of CPC's own factor. With flooring applied, §8's factors reproduce 88.8% of published monthly values exactly against §4's 83.5%. `docs/CPC_RONI_scaling.md` Tables 2-3 are built from the §4 factors and are labeled as such |
| Factor stability across the record | Not sensitive to `n_recent` | Measured 2026-09-10: full-record (1950-2026) per-month fits differ from the 20-year fits by at most 0.007, and a 1950-2005-only fit gives nearly identical factors. Supports the `n_recent=20` choice above, and means the recipe is not a recent-decades artifact |

## Edge Cases & Error Handling

- **`YR` label for Dec-centered seasons**: CPC labels each season by its center month's own calendar year (e.g. NDJ 1997 = Nov97,Dec97,Jan98, center month Dec97 → "1997"; DJF 1998 = Dec97,Jan98,Feb98, center month Jan98 → "1998") — verified against `RONI.ascii.txt`'s 1997-98 El Niño rows. `load_ersstv6_3mo_anomalies()` sets `YR = year` directly from the center-month timestamp, which already has this convention baked in — no adjustment needed.
- **Near-zero `diff_3mo` years blow up simple ratios**: `backout_factors()`'s `factor_mean_all_years` diagnostic and `factor_table_by_year()`'s per-year grid both show this directly (e.g. 1967 JAS ratio = 116.51, 1986 MJJ = -4.44) — expected and documented, not a bug. The single-row (largest-\|diff\|) and OLS-through-origin estimators are unaffected since both naturally downweight near-zero-`diff` observations.
- **Residual structure, resolved 2026-09-10**: the direct monthly fit (§4) leaves a pooled RMS of 0.0054 across 240 fitted points with a systematic negative bias (mean -0.0038). §8 identifies this as CPC's publication convention — the published values are floored to two decimals, not rounded to nearest — which accounts for the whole bias and drops the residual to the 0.0029 quantization floor (see Constants). Earlier hypotheses now ruled out: seasonal-vs-monthly aggregation (§4 fits and targets monthly throughout), area weighting (§10), and a variance-matching-period mismatch (which would rescale the factor, not offset it). What the earlier spec described as a negative correlation with the published anomaly is the same effect seen through a fit forced through the origin.
- **Residual drift across decades**: under the §8 refit factors the mean residual is not flat, running from -0.0065 in the 1950s to -0.0041 in the 2020s where a pure publication convention would hold it at -0.005 throughout (`residual_by_decade.csv`). A slowly varying ~0.002 difference between this script's inputs and CPC's, such as an ERSSTv6 revision, would produce it. Left as an Open Item.
- **End-of-record values can miss by more than one unit in the last digit**: with flooring applied, 817 of 920 monthly values reproduce exactly and all but one of the rest are off by exactly 0.01. The exception is June 2026 (published 1.07 vs. 1.047 reconstructed), which also drags AMJ 2026 out by two units — both at the very end of the record, where ERSSTv6 is still subject to revision after CPC computes its published value. Expected, not a bug.
- **2026 partial coverage**: `RONI.ascii.txt`/`Rnino34.ascii.txt` run through JJA/Aug 2026 at time of writing, so DJF-JJA (or Jan-Aug) use published years 2007-2026 in the 20-year windows while JAS-NDJ (or Sep-Dec) use 2006-2025 (one year lagged) — not a script bug, just the data's cutoff (`.tail(n_recent)`/year-cutoff logic handles this per month/season automatically).

## Open Items

- **Flooring vs. a genuine offset is not decided by the data** (see Constants): flooring with no offset, and round-to-nearest applied to a series that genuinely sits 0.005 below this script's, put the residual in exactly the same place. The argument for flooring rests on the fitted offset landing at precisely half a publication step. Deciding it would need something outside these two products — CPC's own code, or a statement of the convention.
- **Decadal drift in the residual** (see Edge Cases): -0.0065 in the 1950s to -0.0041 in the 2020s, ~0.002 of unexplained slow variation. Candidate: an ERSSTv6 revision between CPC's computation and the locally cached file. Decision to make: whether it is worth chasing at all, given it sits below the 0.01 publication precision.
- **Two cosmetic items in `docs/CPC_RONI_scaling.md`, found 2026-09-10, not fixed**: (1) §4 compares a mean of 12 per-season RMS values (0.0065) against a pooled RMS (0.0049); like-for-like the seasonal-block fit is 46% worse pooled (0.00714 vs. 0.00490) or 34% worse by mean of seasons, so the note understates its own case with "about a third higher". (2) The note never states that the factors are stable over the full record (see Constants), so the recipe reads as validated only on 2006-2026.
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
    rounding_convention_check, seasonal_convention_check,
    reconstruction_tables, residual_by_decade, weighting_sensitivity,
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

# Publication convention: flooring must beat round-to-nearest by a wide margin,
# and the fitted intercept must sit near half a 0.01 publication step.
conv, refit = rounding_convention_check(n_recent=20)
full = conv[conv["sample"] == "1950-present"]
assert (full.floor_exact > 0.8).all(), "flooring should reproduce most published values exactly"
assert (full.floor_exact > full.round_exact + 0.25).all(), "flooring should beat round-to-nearest"
assert (full.offset_fit.between(-0.006, -0.004)).all(), "offset should sit near -0.005"
assert refit.delta.abs().max() < 0.01, "refit should barely move the factors"

factors_refit = dict(zip(range(1, 13), refit.factor_offset_refit.values))
seas = seasonal_convention_check(factors_refit)
assert (seas.floor_exact > seas.round_exact).all()

_, recon, stats, _ = reconstruction_tables(factors_refit, start_year=2020)
assert stats.max_abs_err.iloc[0] < 0.02, "reconstruction should track published within one step"
assert residual_by_decade(factors_refit).mean_resid.max() < 0, "residual bias is one-signed"

w = weighting_sensitivity(n_recent=20)
assert (w.delta < 0).all(), "unweighted factors run below cos-lat-weighted ones"

for f in ["report.txt", "monthly_seasonal_consistency.csv", "monthly_factors_direct.csv",
          "direct_monthly_vs_seasonal_rms.csv", "seasonal_factors_singlerow.csv",
          "seasonal_factors_recent.csv", "monthly_factors.csv", "seasonal_vs_monthly_rms.csv",
          "factor_by_year.csv", "rounding_convention.csv", "rounding_convention_seasonal.csv",
          "monthly_factors_offset_refit.csv", "residual_by_decade.csv",
          "weighting_sensitivity.csv", "seasonal_published_recent.csv",
          "seasonal_reconstructed_recent.csv", "seasonal_reconstruction_stats.csv"]:
    assert (config.PLOTS_DIR_RONI_FACTOR_BACKOUT / f).exists(), f"missing output {f}"

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-09-10 | Initial `scripts/roni_factor_backout.py` written: seasonal-block factor backout (single-row + OLS-through-origin variants), then superseded/supplemented by a joint 12-calendar-month linear fit (per user correction that CPC scales monthly values before seasonal averaging) — reduced worst-season RMS (AMJ/MJJ) from 0.0124-0.0128 to 0.0044-0.0056. Tested and ruled out area-weighting as the source of the small remaining residual. Added `config.ERSSTV6_NC`, `config.RONI_TXT`, `config.PLOTS_DIR_RONI_FACTOR_BACKOUT`. Script now saves its report + all tables to `plots/roni_factor_backout/`. | ✓ |
| 2026-09-10 | Correctness review of `docs/CPC_RONI_scaling.md` and this script (every quoted number re-derived independently; all reproduced). Three outcomes. (1) **Residual resolved**: added `rounding_convention_check()`, `seasonal_convention_check()`, `residual_by_decade()` and the `_fit_slopes_with_offset()`/`_quantize_match()` helpers (§8) — CPC floors its published values to two decimals rather than rounding to nearest, which accounts for the entire -0.004 bias; flooring reproduces 88.8% of published monthly values exactly vs. 48.7% for round-to-nearest. Write-up gained a new §5 and a `round down` step in the recipe. (2) **Weighting claim corrected**: `weighting_sensitivity()` (§10) shows cos-lat vs. plain grid-cell means differ by 0.0023-0.0034 in the factor, not the <0.0001 the note claimed (that figure was the RMS difference). (3) **Note made regenerable**: added `monthly_fit_frame()`, `_recent_years()`, `reconstruct_from_factors()`, `reconstruction_tables()` (§9) plus pooled fit statistics, so every number the write-up quotes is now emitted to `report.txt` and CSVs. Also fixed `validate_direct_monthly_against_seasonal()` silently ignoring `weighted=False` (now takes the parameter), and rewrote the stale module docstring, which still described the superseded per-season estimator as the method. Section 2's rounding-bound claim was checked and confirmed correct; its argument is now spelled out in the note rather than changed. | ✓ |
| 2026-09-10 | Same-day follow-up, three additions: (1) `check_monthly_seasonal_consistency()` — confirmed CPC's published monthly and seasonal products are mutually consistent (max diff 0.0067, mean 0.00007) using no ERSSTv6 data at all, establishing the monthly-then-average order of operations directly rather than assuming it. (2) `backout_monthly_factors_direct()` + `validate_direct_monthly_against_seasonal()` — fit each month's factor directly against the newly-added `Rnino34.ascii.txt` (published monthly series), now the primary factor estimate, superseding the seasonal-aggregate-only joint fit; confirmed close agreement with it (within ~0.01 per month) and confirmed the residual is not a seasonal-aggregation artifact. (3) Cross-checked the recipe (box means, climatology, order of operations) against M. L'Heureux's public reference notebook (github.com/michellelheureux/Relative-SST/blob/main/compute_relativeNino34_ONI_JClimate.ipynb, harvested for reference only, not incorporated into the codebase) — confirmed unweighted box means and the monthly-then-average order, no changes to the recipe resulted. Added `config.RNINO34_TXT`. Wrote `docs/CPC_RONI_scaling.md` (external-reader write-up of all three checks). | ✓ |
