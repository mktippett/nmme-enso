# skill.py — Behavioral Specification

> Last reviewed against code: 2026-07-06 (initial version)

## Purpose

Quantify NMME forecast skill at predicting the Niño-3.4 index itself:
anomaly correlation (AC) and MSESS of the ensemble-mean forecast against
ERSSTv5, as a function of calendar month and lead, for each model and for
the multi-model mean (MMM). These are diagnostic/exploratory figures (no
manuscript caption) establishing the baseline skill structure — event-based
verification and reliability diagnostics are separate, later scripts.

## Inputs

| Source | Relevant variables | Filters applied |
|--------|--------------------|------------------|
| `config.load_nino34_verification()` | `ssta` (model, S, M, L), `obsa` (S, L) | Restricted to starts `S` in `[CLIM_START_YEAR-01-01, CLIM_END_YEAR-12-01]` (1991-2020) for both AC and MSESS — a single common sample, not the notebook's original AC-only-1991-2020 / MSESS-all-starts split |

`load_nino34_verification()` itself wraps `config.load_nino34_ssta()` (NMME
ensemble forecast anomalies) with ERSSTv5 observational anomalies aligned
onto the same `(S, L)` grid — see Algorithm §1 and that function's
docstring for the full alignment rationale.

## Outputs

| File | Contents | Format |
|------|----------|--------|
| `plots/skill/n34_ac_start.png` | Anomaly correlation, start month × lead, 8 panels (7 models + MMM), `vmin=0, vmax=1` | PNG, dpi=150 |
| `plots/skill/n34_ac_target.png` | Same, relabeled to target month (= start + lead) | PNG, dpi=150 |
| `plots/skill/n34_msess_start.png` | MSESS, start month × lead, 8 panels, `vmax=1` (no `vmin` — negative/no-skill shows) | PNG, dpi=150 |
| `plots/skill/n34_msess_target.png` | Same, relabeled to target month | PNG, dpi=150 |

Each figure: one panel per model (short label via `config.short_label`) plus
an 8th "MMM" panel, faceted `col_wrap=4` (a clean 2×4 grid — no blank
slot). X-axis = month (Jan-Dec abbreviations), y-axis = lead (months, 0-11
displayed as 0-12 tick positions per `L`'s half-integer values), colorbar =
"correlation" / "MSESS". Suptitle carries the framing (start vs. target
month) and period, since these are exploratory figures without a caption.

## Algorithm

### 1. Observation alignment (in `config.py`, not this script)

`config.load_nino34_verification()` aligns ERSSTv5 (Gregorian, monthly, day
stamped on the 1st) onto the model `(S, L)` grid, where `target = S + L` is
a 360-day cftime valid time stamped on the 16th (`L=0.5` = the init month
itself). Both grids are collapsed to an integer `(year, month)` index
— `ym = (year - 1960) * 12 + (month - 1)` — and the gather is an exact-
integer `.sel(ym=...)`, rather than a float/timestamp join on the two
grids' differing day-of-month and calendar conventions (which would need a
fabricated 360-day mid-month obs axis and floating-point `.sel` equality —
fragile to any half-month or 1st-vs-16th convention drift). See that
function's docstring for the full rationale. Observational anomalies are
computed over ERSSTv5's full record (climatology still fixed to
1991-2020), so `obsa` is defined for every target the model's `S`/`L` grid
reaches; only genuine future forecasts (targets beyond the end of the
ERSSTv5 record) are NaN.

### 2. Anomaly correlation — `ac_by_start(x, y)`

For forecast `x` and observation `y` (already anomalies), masked to where
both are non-missing (so models with different lead counts don't bias the
mean/std), stratified by `S.month`:
```
ok = ~isnan(x) & ~isnan(y)
xa = (x.where(ok) - mean_S(x.where(ok)))   # per start-month group
ya = (y.where(ok) - mean_S(y.where(ok)))   # per start-month group
AC = mean_S(xa * ya) / std_S(xa) / std_S(ya)   # per start-month group
```
`x = ds.ssta.mean("M")` (ensemble-mean forecast, dims model/S/L); `y =
ds.obsa` (dims S/L, broadcasts over model).

### 3. MSESS — `msess_by_start(f, o)`

Mean-square-error skill score vs. a zero-anomaly reference forecast,
stratified by `S.month`:
```
MSESS = 1 - mean_S((f - o)^2) / mean_S(o^2)
```
`.mean("S")` (not `.sum()`) so a (month, lead) cell is NaN only if *all*
its samples are missing (e.g. a model with a shorter lead count stops
contributing at higher leads), never silently treated as zero.

### 4. Multi-model mean (MMM) — `add_mmm`

`mmm_forecast = forecast_mean.mean("model")` — equal-weight mean across
models of the already-per-model ensemble-mean forecast. Its own skill is
computed with the same `ac_by_start`/`msess_by_start` functions, tagged
`model="MMM"`, and `xr.concat`-ed onto the per-model skill array along
`model` — giving 8 panels instead of 7, and letting the MMM be compared
directly against individual models on the same color scale.

### 5. Start → target month relabel — `to_target_month(sk)`

For each lead `l` (half-integer months, `L=0.5` = the init month itself),
target month = start month + `round(l - 0.5)`, computed via
`sk.sel(L=l).roll(month=round(l-0.5), roll_coords=False)` and
concatenated back over `L`. A functional per-lead build-and-concat, not an
in-place `.loc[...] = ...roll(...)` mutation loop — same result, no mutable
intermediate copy.

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| Skill period | `CLIM_START_YEAR`-`CLIM_END_YEAR` (1991-2020) for both AC and MSESS | One common sample for every metric in this script (deliberate choice — the source notebook used 1991-2020 for AC but the full record for MSESS; this script uses the same window for both for consistency) |
| MSESS reference forecast | Zero anomaly (climatology) | Standard skill-score baseline; `msess_by_start` assumes both `f` and `o` are already anomalies against their own climatologies |
| MMM weighting | Equal-weight mean across models (after each model's own member mean) | Simplest multi-model combination; no skill-based or variance-based model weighting attempted here |
| Facet layout | `col_wrap=4` → 2×4 | 7 models + MMM = 8 panels exactly fills a 4-column, 2-row grid with no blank panel |

## Edge Cases & Error Handling

- **Models with shorter lead counts** (`NASA-GEOSS2S`, `NCEP-CFSv2`):
  NaN-padded at longer leads in the store; `.mean("S")` in both skill
  functions returns NaN only for (month, lead) cells with zero valid
  samples, so their panels simply stop (rendered as a blank region by
  `xarray`'s auto color mapping) rather than showing a spurious value.
- **Forecast targets beyond the ERSSTv5 record** (`obsa` NaN via
  `load_nino34_verification`'s `reindex`): excluded from both skill
  functions by the same `~isnan` masking (AC) / all-NaN-collapses-to-NaN
  behavior (MSESS) — not specific to this script.
- **MMM panel look**: expected to show *higher* skill than most individual
  models (ensemble-of-ensembles averaging reduces noise) — this is a
  scientifically expected result, not a bug, and was confirmed visually on
  the first production run (see Synchronization Log).

## Verification Snippet

```python
# Run after changes to confirm key invariants
import sys
sys.path.insert(0, "scripts")
import config
import numpy as np
from skill import ac_by_start, msess_by_start, add_mmm, to_target_month

ds = config.load_nino34_verification()
ds_p = ds.sel(S=slice(f"{config.CLIM_START_YEAR}-01-01", f"{config.CLIM_END_YEAR}-12-01"))
forecast_mean = ds_p.ssta.mean("M")
obs = ds_p.obsa

ac = ac_by_start(forecast_mean, obs)
assert set(ac.dims) == {"model", "month", "L"}
assert (ac.max() <= 1.0 + 1e-6) and (ac.min() >= -1.0 - 1e-6), "correlation out of [-1,1]"

ac_mmm = add_mmm(ac, forecast_mean, obs, ac_by_start)
assert "MMM" in ac_mmm.model.values
assert ac_mmm.sizes["model"] == ac.sizes["model"] + 1

msess = msess_by_start(forecast_mean, obs)
assert set(msess.dims) == {"model", "month", "L"}

ac_target = to_target_month(ac_mmm)
assert ac_target.month.attrs["long_name"] == "target month"
assert ac_target.sizes == ac_mmm.sizes

for f in ["n34_ac_start.png", "n34_ac_target.png", "n34_msess_start.png", "n34_msess_target.png"]:
    assert (config.PLOTS_DIR_SKILL / f).exists(), f"missing output {f}"

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|---------------|
| 2026-07-06 | Initial `scripts/skill.py` + `config.load_nino34_verification()` (and supporting `config._obs_index_anomalies()` refactor, shared with `load_obs_scaling`) written | ✓ |
| 2026-07-07 | Output moved from `plots/` to `plots/skill/` (per-script subdirectory, `config.PLOTS_DIR_SKILL`) | ✓ |
| 2026-07-07 | `config.ERSSTV5_NC` default changed from `/Users/tippett/notebooks/data/ERSSTv5.sst.mnmean.nc` to the repo-local `OBS_DIR / "ERSSTv5.sst.mnmean.nc"` (bootstrapped/refreshed by the README Quickstart `curl -z` command; env var still overrides). 1991-2020 skill statistics verified unchanged; the fresher PSL file extends obs coverage by one month (June 2026). | ✓ |
