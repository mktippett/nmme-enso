# latest_forecast.py — Behavioral Specification

> Last reviewed against code: 2026-07-22 (added synthetic error-covariance spread family)

## Purpose

Monitor the most recent NMME Niño-3.4 forecast: show what each model
currently predicts, how spread out the ensemble is, and how the latest
initialization compares with the previous one. Produces the routine
forecast-watch figures — not manuscript figures. Two indices are plotted,
each as its own figure set:

- **Niño-3.4** (`n34_*`): the standard box anomaly.
- **relative Niño-3.4** (`n34r_*`): Niño-3.4 anomaly minus the tropical-mean
  anomaly, rescaled by a `(model, start-month, lead)`-dependent factor so
  the *model's own* relative-index variance matches the *observed*
  (ERSSTv5) Niño-3.4 variance over 1991-2020 — see L'Heureux, Tippett et
  al. (2024, *J. Climate*, `papers/`) for the index's origin. The specific
  model-relative, (month, lead)-stratified factor is an extension beyond
  the paper's own obs-only ratio (per the first author, 2026-07-06; see
  Constants & Scientific Rationale).

For each index, the four plume figures (compare, spread, spread-synthetic,
mean) are each produced in a **monthly** and a **seasonal** (3-month running
mean) variant, plus one grid figure (monthly only) — 9 outputs per index, 18
total. Spread-synthetic (§4b) replaces the real (underdispersed, model-bias-
dominated) ensemble members with a calibrated Gaussian plume drawn from the
historical MMM forecast-error covariance across leads (Barnston, Tippett,
van den Dool & Unger 2015, *J. Appl. Meteor. Climatol.*, **54**, 1579–1595,
https://doi.org/10.1175/JAMC-D-14-0188.1, Fig. 9 lower panels).

## CLI Arguments

| Flag | Default | Effect |
|------|---------|--------|
| `--init-date YYYY-MM[-DD]` | none (latest init) | Plot this initialization instead of the latest, matched by calendar year/month (day ignored). Raises `ValueError` if no init in the store matches that month. All 18 output filenames get a `_<YYYY-MM>` suffix (see Outputs), so an explicit-init run never overwrites the default latest-init run's files. |

## Inputs

| File | Relevant columns | Filters applied |
|------|-------------------|-----------------|
| `$NMME_STORE_DIR/nmme_sst.zarr` (`config.STORE_SST`) | `sst(S, M, L, Y, X)`, all 7 model groups | Reduced to the Niño-3.4 box average and the tropical-mean box average via `config.load_nino34_ssta()`; no further row filtering — all available `S` (init times) are loaded, only the latest one or two are plotted |
| `config.ERSSTV5_NC` (ERSSTv5 monthly SST) | `sst(time, lat, lon)` | Reduced to a Niño-3.4 box average via `config._obs_index_anomalies()`, aligned onto the `(S, L)` grid as `ds.obsa` by `config.load_nino34_verification()`; used as the numerator of `config.rel_scaling_factor()`'s scaling factor, restricted to `[CLIM_START_YEAR, CLIM_END_YEAR]` (1991-2020) |

### Caching

`config.load_nino34_ssta()` (wrapped by `config.load_nino34_verification()`,
which this script calls) caches to disk (default `use_cache=True`) since the
~76GB SST store's box-average read dominates runtime (~1 min cold) while the
underlying data barely changes between runs — `update_archive.py` (in
`~/claude/NMME-zarr`) only appends new forecast starts and rechecks its
trailing `RECHECK_TAIL=2` starts.

- `load_nino34_ssta()` → `cache/nino34_ssta.nc`, keyed on
  `config._store_fingerprint(store)`: per model group, `S` array size, last
  `S` value, and the `last_updated` zarr attribute — cheap metadata reads,
  never the sst data itself. A mismatch (new/rechecked starts) triggers a
  full recompute and cache rewrite.
- `config.rel_scaling_factor()` has **no disk cache of its own** — it's a
  `groupby('S.month').std('S')` reduction over the already-cached (small)
  verification Dataset, cheap enough to recompute every run.
- `cache/` is gitignored (regenerable, instance-specific). Pass
  `use_cache=False` to `load_nino34_ssta()`/`load_nino34_verification()` to
  force a fresh read.

## Outputs

Filenames below are for the default (no `--init-date`) run. When
`--init-date` is passed, every filename gets a `_<YYYY-MM>` suffix before
`.png` (e.g. `n34_monthly_grid_2026-06.png`) — the selected init's
year-month, not the value typed on the command line.

| File | Contents | Format |
|------|----------|--------|
| `plots/latest_forecast/n34_monthly_grid.png` | Facet grid (one panel per model) of the selected-init Niño-3.4 anomaly, lead x member | PNG, dpi=150 |
| `plots/latest_forecast/n34_monthly_compare.png` | Ensemble-mean plume: selected init (solid) vs. previous init (dashed), monthly | PNG, dpi=200 |
| `plots/latest_forecast/n34_monthly_spread.png` | All members (thin) + ensemble mean (thick), selected init only, monthly | PNG, dpi=200 |
| `plots/latest_forecast/n34_monthly_spread_synthetic.png` | 100 synthetic members + MMM + 10th/90th percentile, drawn from the historical MMM error covariance across leads (see §4b), selected init only, monthly | PNG, dpi=200 |
| `plots/latest_forecast/n34_monthly_mean.png` | Ensemble-mean-only plume, selected init only, monthly | PNG |
| `plots/latest_forecast/n34_seasonal_compare.png` | Same as `n34_monthly_compare.png`, but each series is a centered 3-month running mean over lead, x-axis labeled by target season | PNG, dpi=200 |
| `plots/latest_forecast/n34_seasonal_spread.png` | Same as `n34_monthly_spread.png`, seasonal (members and ensemble mean both smoothed) | PNG, dpi=200 |
| `plots/latest_forecast/n34_seasonal_spread_synthetic.png` | Same as `n34_monthly_spread_synthetic.png`, seasonal (error covariance computed in the smoothed space, see §4b) | PNG, dpi=200 |
| `plots/latest_forecast/n34_seasonal_mean.png` | Same as `n34_monthly_mean.png`, seasonal | PNG |
| `plots/latest_forecast/n34r_monthly_grid.png` | Same as `n34_monthly_grid.png`, for the scaled relative Niño-3.4 anomaly | PNG, dpi=150 |
| `plots/latest_forecast/n34r_monthly_compare.png` | Same as `n34_monthly_compare.png`, relative index | PNG, dpi=200 |
| `plots/latest_forecast/n34r_monthly_spread.png` | Same as `n34_monthly_spread.png`, relative index | PNG, dpi=200 |
| `plots/latest_forecast/n34r_monthly_spread_synthetic.png` | Same as `n34_monthly_spread_synthetic.png`, relative index — error covariance computed in the scaled `ssta_rel` space and verified against **observed absolute** Niño-3.4 (`ds.obsa`), not `ds.obsa_rel` (see §4b) | PNG, dpi=200 |
| `plots/latest_forecast/n34r_monthly_mean.png` | Same as `n34_monthly_mean.png`, relative index | PNG |
| `plots/latest_forecast/n34r_seasonal_compare.png` | Same as `n34_seasonal_compare.png`, relative index (own seasonal scaling factor) | PNG, dpi=200 |
| `plots/latest_forecast/n34r_seasonal_spread.png` | Same as `n34_seasonal_spread.png`, relative index | PNG, dpi=200 |
| `plots/latest_forecast/n34r_seasonal_spread_synthetic.png` | Same as `n34_seasonal_spread_synthetic.png`, relative index (own seasonal scaling factor, verified against `ds.obsa`) | PNG, dpi=200 |
| `plots/latest_forecast/n34r_seasonal_mean.png` | Same as `n34_seasonal_mean.png`, relative index | PNG |
| `plots/latest_forecast/latest_forecast_summary.md` | Monthly and seasonal MMM anomaly tables for both indices (n34, n34r), each value's rank (1 = highest) among all MMM forecasts issued in the same calendar start month, `ANALYSIS_START_YEAR`-present — see §5 | Markdown |

## Algorithm

### 1. Load

Call `config.load_nino34_ssta()`. This is the shared loader (used by any
future script needing the Niño-3.4 index/anomaly) and does the following,
**in this order** — the order is load-bearing, not stylistic:

1. For each model group in `config.nmme_groups()`, open with
   `decode_times=False` (raw float `S`), reduce `sst` to the Niño-3.4 box
   average (cosine-latitude-weighted mean over `config.N34_LON` x
   `config.N34_LAT`) *and* to the tropical-mean box average
   (cosine-latitude-weighted mean over `config.TROPICS_LAT`, all
   longitudes; land points are NaN and skipped by the weighted mean), tag a
   `model` coordinate, and merge all groups.
2. Compute `target = L + S` **while `S` is still a raw float** — this must
   happen before any cftime decoding, since decoding converts `S` to a
   cftime object and `L + S` would no longer be a plain float offset.
3. Decode once, via `xr.decode_cf`, after patching the calendar attribute
   from the store's `'360'` to CF-compliant `'360_day'`. This single call
   decodes both `S` and `target` (both carry matching `units`/`calendar`
   attrs) to cftime `Datetime360Day`.
4. Compute forecast anomalies for both box averages (`ssta`, `ssta_trop`)
   via the shared helper `config._forecast_anomaly()`, which applies **two
   distinct climatology schemes** — see Constants table — identically to
   both. This asymmetry (different baseline windows *and* different time
   coordinates) reflects a real NMME hindcast discontinuity in
   `config.TWO_CLIM_GROUPS` and must not be collapsed into one scheme.
5. `ssta_rel = ssta - ssta_trop` — the **unscaled** relative Niño-3.4
   anomaly. Because both terms are anomalies against the same monthly
   climatology, the (differing) seasonal cycles of the two regions cancel
   rather than leaking into `ssta_rel`.

### 1a. Model-relative scaling factor

`main()` loads with `config.load_nino34_verification()` (not
`load_nino34_ssta()`) so `ds` carries `ds.obsa` — ERSSTv5 Niño-3.4 anomalies
aligned onto the same `(S, L)` grid as the forecasts (see
`specs/skill.md`/`config.load_nino34_verification` docstring for the
year-month alignment mechanism) — then calls `config.rel_scaling_factor(ds)`
once (not inside the loader, since it is a diagnostic quantity derived from
both obs and model data, not a property of the raw forecast itself):

1. Restrict to forecast starts `S` in `[CLIM_START_YEAR, CLIM_END_YEAR]`
   (1991-2020) via `ds.where(..., drop=True)` on `S.dt.year`.
2. Numerator: `std_obs(n34) = ref.obsa.groupby('S.month').std('S')` — dims
   `(month, L)`, no model dependence (obs has no model dim).
3. Denominator: `std_model(n34r) = sqrt(ref.ssta_rel.groupby('S.month').var(['S', 'M']))`
   — dims `(model, month, L)`. This pools ensemble members via the
   **grand mean — the pooled standard deviation over the flattened
   (start, member) sample** — deliberately *not* `std('S')` of the ensemble
   mean (`ssta_rel.mean('M')`), which would shrink wherever skill is low
   and mix predictability into what should be a pure variance-matching
   factor.

   There are **two** valid member-based pooling estimators here, equal in
   expectation but differing in finite samples: this grand-mean form
   (`sqrt(var(['S', 'M']))`, used since 2026-07-07), and a per-member form
   (`sqrt(var('S').mean('M'))`, the per-member temporal variance averaged
   across members — the original choice) — the grand-mean form additionally
   carries the between-member spread of the per-member time-means (law of
   total variance), which the per-member form averages out. Both differ
   from the ensemble-mean estimator above, which is not just
   finite-sample-different but systematically flawed (skill-contaminated).
   See Constants & Scientific Rationale for why the grand-mean estimator is
   now used, and `scripts/rel_scaling_compare.py` for the supporting
   evidence comparing all three pairwise (per-member vs. ensemble-mean,
   per-member vs. grand-mean, ensemble-mean vs. grand-mean).
4. `factor_monthly(model, month, L) = std_obs(n34) / std_model(n34r)`.
5. `factor_seasonal(model, month, L)`: identical ratio, but both `obsa` and
   `ssta_rel` are passed through a centered 3-month running mean over `L`
   first (matching the seasonal plotting variant, since running-meaning
   changes the variance ratio). Endpoint leads (first/last `L`) are NaN
   from the rolling window, same as the seasonal plots' own NaN endpoints.

Returns `(factor_monthly, factor_seasonal)`, each a `DataArray` dims
`(model, month, L)`. The relative-index plot functions multiply the
(optionally seasonally-smoothed) `ssta_rel` by
`factor.sel(model=model, month=start_month)`, where `start_month` is the
calendar month of the forecast's *init* time (`ds.S.dt.month`), **not**
target/valid month — scaling is a plot-time transform, `ssta_rel` in the
loaded Dataset stays unscaled.

### 2. Init-time axis and init selection

`start = ds.indexes['S'].to_datetimeindex(time_unit='ms')` converts the
360-day cftime `S` index to a Gregorian pandas `DatetimeIndex` (with an
xarray-flagged "unsafe" calendar conversion, since 360-day and Gregorian
calendars don't align exactly — inits land on the 16th of each month).

`main()` parses an optional `--init-date YYYY-MM[-DD]` CLI argument and
resolves it to a positional `S` index via `_resolve_init_idx(start,
init_date)`:

- `init_date=None` (default, no flag passed) → `now_idx = -1`, the most
  recent init.
- `init_date` given → matched by **calendar year/month only** (the day is
  ignored, since inits land on the 16th, not necessarily the day passed);
  `now_idx` is the first (only) matching position. No match raises
  `ValueError` listing the available `start[0]:start[-1]` range.

`prev_idx = now_idx - 1` throughout (works identically whether `now_idx` is
`-1` or a resolved positive position). Every function that used to hardcode
`S=-1`/`S=-2` or `start[-1]`/`start[-2]` now takes `now_idx`/`prev_idx` (or
just `now_idx`, for the three latest-init-only plots) as an explicit
parameter — there is no implicit "latest" fallback inside the plot
functions themselves.

`date_suffix = f"_{start[now_idx]:%Y-%m}"` when `--init-date` was passed,
else `""` — appended to every output filename for that run (see Outputs).
This keeps the default (no-flag) run's filenames unchanged, so the routine
monitoring figures stay at fixed paths, while an explicit past-init run
never overwrites them.

### 3. Index specs and available models

`main()` loops over two index specs, `{"var": "ssta", "prefix": "n34",
"name": "Nino 3.4", "scale": None}` and `{"var": "ssta_rel", "prefix":
"n34r", "name": "relative Nino 3.4", "scale": (factor_monthly,
factor_seasonal)}`, running the full grid/compare/spread/mean sequence once
per spec.

`avail = np.where(~np.isnan(ds[spec['var']].isel(S=now_idx, L=0).mean('M')))[0]`
— positional model indices with a live (non-NaN) zero-lead forecast at
`now_idx` for that spec's variable. All figures for a given spec loop only
over its `avail`, so a model missing that init's forecast (e.g. delayed
release, or a model not yet running at an older init) is silently skipped
rather than plotted as a gap.

### 4. Figures

- **Grid** — `ds[spec['var']].isel(S=now_idx)` (scaled by
  `factor_monthly.sel(month=start_month)` first, for the relative index,
  where `start_month` is the selected init's calendar month),
  `.plot(col='model', col_wrap=4)`;
  suptitle shows the init date (`ds.S.isel(S=now_idx)`, `str(...)[:10]`).
  Monthly only — there is no seasonal grid variant. Before plotting, the
  array is renamed to `spec['prefix']` and given fresh `attrs = {"units":
  "degC"}` — replacing (not patching) any inherited `long_name`/
  `standard_name` from the raw store field or the ERSSTv5 scale factor,
  which otherwise silently wins over the array name in xarray's auto
  colorbar label.
- **Compare** — per model in `avail`, two plumes over 12 leads each:
  previous init (`S=prev_idx`, dashed, circle marker) and selected init
  (`S=now_idx`, solid, square marker). X-axis is `pd.date_range(start[i],
  periods=12, freq='MS')` per init. Colors from the plume-specific override
  (see Constants). Legend labels via `config.short_label()`. Title uses
  `spec['name']`.
- **Spread** — selected init (`now_idx`) only; all ensemble members thin
  (alpha 0.35) + ensemble mean thick (alpha 0.85). Same color override array
  as Compare (indexed by the same positional model index), not the plain
  default cycle.
- **Spread-synthetic** — selected init (`now_idx`) only; the current MMM
  (built identically to Spread's, via the same per-model `_index_transform`
  loop) plus `N_SYNTHETIC_MEMBERS=100` Gaussian scenarios drawn from the
  historical MMM forecast-error covariance across leads and added to that
  MMM, plus a dashed 10th/90th percentile envelope. See §4b for the full
  covariance/draw procedure. Uses its own colors
  (`SYNTHETIC_MEMBER_COLOR`/`SYNTHETIC_MMM_COLOR`, not the Compare/Spread/Mean
  override array or `MMM_COLOR`) since `MMM_COLOR = "0.75"` would be
  invisible against the synthetic member cloud.
- **Mean** — selected init (`now_idx`) only, ensemble-mean lines. Same color
  override array as Compare/Spread.
- **Multi-model mean (MMM)** — Compare, Spread, and Mean each also plot an
  equal-weight mean across `avail` models (`MMM_COLOR = "0.75"`, light gray),
  computed from the same per-model ensemble-mean arrays already plotted for
  that figure (i.e. `xr.concat([...], dim="model").mean("model")` over the
  per-model `_index_transform` output collected during the loop — same MMM
  definition as `skill.py`'s `add_mmm`: equal-weight mean of each model's own
  member mean, not a member-pooled grand mean). It is plotted **last**, after
  the per-model loop, with no explicit `zorder`, so it renders on top of the
  colored per-model lines (matplotlib's default increasing z-order by call
  order) rather than being buried under them. In Compare this means two MMM
  lines (previous init dashed, selected init solid), matching the per-model
  convention; legend label `"MMM"`. **Not** added to the Grid figure — tried
  and deliberately reverted (see Synchronization Log 2026-07-09) since Grid's
  panels are member x lead heatmaps and MMM has no member dimension of its
  own.

### 4a. X-axis ticks and limits

`_set_xaxis(fig, ax, ticks, seasonal)` is the single place tick positions,
labels, and xlim are set — every plot function computes its own `ticks`
from its own plotted data and passes them in; none reuse another figure's
limits.

- **Tick array per figure**: Compare uses 13 ticks spanning
  `pd.date_range(start[-2], periods=13, freq='MS')` (covering both the
  previous init's 12 leads and the latest init's 12 leads, which overlap by
  11 months). Spread and Mean each use 12 ticks matching their own single
  `leads = pd.date_range(start[-1], periods=12, freq='MS')` — deliberately
  *not* Compare's wider 13-tick range, since Spread/Mean only plot the
  latest init.
- **Seasonal trimming**: in the seasonal variant, `_set_xaxis` drops the
  first and last tick (`ticks[1:-1]`) before setting ticks/labels/xlim,
  because `_seasonal`'s centered rolling mean leaves those two leads NaN
  (see Edge Cases) — without trimming, the outermost tick shows a season
  label with no line reaching it.
- **`_tight_xlim(ticks)`**: xlim is `(ticks[0] - step/2, ticks[-1] + step/2)`
  where `step = ticks[1] - ticks[0]` (one month) — half a time step of
  padding on each side, computed fresh per figure after any seasonal
  trimming, rather than matplotlib's automatic ~5% margin or a shared/reused
  xlim from another figure.
- **Seasonal variant** (Compare, Spread, Mean): each plotted array (ensemble
  means, and in Spread the individual members) is passed through
  `_index_transform(da, spec, seasonal, model, start_month)`, which applies
  `_seasonal(da) = da.rolling(L=3, center=True).mean()` before plotting — a
  centered 3-month running mean over lead — and, for a spec with a `scale`,
  multiplies by `factor.sel(model=model, month=start_month)` (using
  `factor_seasonal` instead of `factor_monthly` when `seasonal=True`, since
  running-meaning changes the variance ratio). `start_month` is the
  **init's** calendar month (`ds.S.isel(S=...).dt.month`), computed once per
  init (`now_idx`/`prev_idx`) outside the per-model loop; the factor
  additionally selects on `model` since the denominator is model-dependent
  (unlike the prior obs-only factor, which was the same for every model).
  The factor multiply commutes with the ensemble `mean('M')` (it depends
  only on model and lead, not member), so it can be applied before or after
  that reduction without changing the result. The x-axis keeps the same
  monthly datetime tick positions, but tick *labels* are replaced with
  `_season_label(month)`, a 3-letter overlapping-season code centered on
  that tick's month (e.g. January -> `DJF`, December -> `NDJ`), using
  `SEASON_INITIALS = "JFMAMJJASOND"`. The start/end lead marker moves from
  lead 0 to lead 1 (`l0 = 1`), since a centered 3-month window leaves lead 0
  and the last lead undefined (see Edge Cases). Titles/output filenames
  distinguish the two kinds and the two indices;
  `config.load_nino34_verification()` itself is unaffected — smoothing and
  scaling are both plotting-time transforms only.

### 4b. Synthetic error-covariance plume (`_synthetic_plume`)

Implements Barnston, Tippett, van den Dool & Unger (2015, *J. Appl. Meteor.
Climatol.*, **54**, 1579–1595, https://doi.org/10.1175/JAMC-D-14-0188.1,
Fig. 9 lower panels): "generate a plume of equally likely scenarios ...
using a Gaussian random number generator, by employing the MME mean
forecast in combination with the historical covariance of the errors over
the hindcast period." Computed
fresh for each figure (index x monthly/seasonal), **in that figure's own
final transformed/scaled space**, so the seasonal covariance is not derived
from a shared monthly draw:

1. **Historical MMM forecast**, over the `avail` models used by the current
   plume, restricted to hindcast starts `year(S) in [CLIM_START_YEAR,
   CLIM_END_YEAR]` (1991-2020; the same predicate as `rel_scaling_factor`).
   This exact period was chosen **because every NMME model has complete
   forecast coverage there** (verified 360/360 starts per model at zero
   lead — `GFDL-SPEAR`, the newest model, starts exactly in 1991), not
   merely because it's the standard climatology window. `mean('model')`
   later in this step is skipna, so if a future model with a shorter
   hindcast joined `avail_models`, it would silently drop out of part of
   the historical sample instead of raising — `_synthetic_plume` guards
   this explicitly: before building `fc`, it checks every `avail_models`
   member has a non-NaN zero-lead forecast at every hindcast start, and
   raises `ValueError` naming any model that doesn't. Then:
   `fc = ds[spec['var']].sel(model=avail_models).mean('M').where(hindcast, drop=True)`
   → `(model, S, L)`. If `seasonal`, roll `fc` (`L=3, center=True`) before
   any scaling — matching `_index_transform`'s order. For the relative
   index, multiply by `factor.sel(model=avail_models, month=fc.S.dt.month)`
   (`factor_seasonal` if seasonal, else `factor_monthly`) — same factor
   already applied to the plotted plume. `mmm_hist = fc.mean('model')` →
   `(S, L)`.
2. **Observed reference**: `obs = ds.obsa.where(hindcast, drop=True)`
   (rolled the same way if seasonal). **Both n34 and n34r verify against
   observed absolute Niño-3.4** (`ds.obsa`), *not* `ds.obsa_rel` — consistent
   with `rel_scaling_factor`, which calibrates the model relative index to
   observed absolute Niño-3.4 variance (see Constants).
3. **Error**: `err = mmm_hist - obs` → `(S, L)`, used **demeaned** —
   `np.cov` subtracts the per-lead sample mean, so any residual MMM
   conditional bias is deliberately *not* injected into the synthetic
   members; the plume stays centered on the *current* MMM, not a
   bias-corrected one.
4. **Stratify by the current start month**: `err_month = err.where(err.S.dt.month
   == start_month_now, drop=True)` → `(S~30, L)` (≈30 hindcast years sharing
   the same init calendar month as the current forecast).
5. **Drop leads with any missing sample**: `err_valid =
   err_month.dropna('L', how='any')` — labeled, not positional; keeps the
   surviving `L` coordinate (`valid_L`). Monthly: expect 0 dropped (12/12
   leads, since `ds.obsa` is fully populated for 1991-2020 starts and
   `mean('model')` is skipna). Seasonal: expect 2 dropped (the rolling
   mean's NaN first/last lead, same two leads Spread's seasonal variant
   already drops). A print warning fires if the dropped count differs from
   this expectation.
6. **Covariance** (numpy boundary): `cov = np.cov(err_valid.transpose('S',
   'L').values, rowvar=False)` → `(nL, nL)`. Diagonal = per-lead error
   variance (≈ SEE²); off-diagonals = lead-to-lead error correlation — this
   is what gives the synthetic members realistic lead-to-lead coherence
   instead of independent per-lead noise.
7. **Draw**: `rng = np.random.default_rng(SYNTHETIC_SEED)`, `draws =
   rng.multivariate_normal(np.zeros(nL), cov, size=N_SYNTHETIC_MEMBERS)` →
   `(100, nL)`. Fixed seed (`SYNTHETIC_SEED=0`) so a rerun with unchanged
   inputs reproduces byte-identical draws.
8. **Add to the current MMM**: wrap `draws` as a DataArray on `(member,
   L=valid_L)`, add `mmm.sel(L=valid_L)` (labeled align — `mmm` is the same
   `(L,)` array Spread-synthetic already built), then `.reindex(L=mmm['L'])`
   back onto the full 12-lead axis (NaN at any dropped lead, e.g. the
   seasonal endpoints) → synthetic `(member, L)`.

The plotted 10th/90th percentile envelope is `synthetic.quantile([0.1, 0.9],
dim='member')` (labeled, skipna) — at the seasonal endpoint leads (all-NaN
across members) this raises a benign `RuntimeWarning: All-NaN slice
encountered` and correctly returns NaN (see Edge Cases).

**Deliberate deviations from raw Spread**: (a) the plume is *not* meant to
reproduce the total multi-model spread, which mixes model-specific climate
drift/bias with genuine forecast uncertainty (compare `n34_monthly_spread.png`,
where inter-model separation — e.g. GEOSS2S trending to ~4.7°C vs. CanESM5 to
~1°C by Feb — dwarfs any single model's own member spread); it reproduces the
*MMM's own* out-of-sample verification uncertainty instead, which is
typically narrower at short leads. (b) No bias correction is applied — the
plume centers on the current (as-is) MMM, not a debiased one.

### 5. Summary tables (`write_summary_tables`)

After the grid/compare/spread/mean figures are written for both index specs,
`main()` calls `write_summary_tables(ds, start, now_idx, index_specs,
avail_by_prefix, date_suffix)`, where `avail_by_prefix` collects each spec's
`avail` (the same fixed model set used for that spec's figures) during the
existing per-spec loop. Produces one markdown file with two tables —
monthly and seasonal — each with 4 rows (`n34 anom`, `n34 rank`, `n34r anom`,
`n34r rank`) and columns = target period (see Outputs).

1. **Historical MMM pool** (`_historical_mmm(ds, avail, spec, seasonal)`):
   `fc = ds[spec['var']].sel(model=avail).mean('M')` → `(model, S, L)`, over
   **all** available `S`, not just the hindcast period — unlike
   `_synthetic_plume`'s 1991-2020-only `fc`. If `seasonal`, roll first
   (`L=3, center=True`); then, for a spec with a `scale`, multiply by
   `factor.sel(model=avail, month=fc.S.dt.month)` (`factor_seasonal` if
   seasonal else `factor_monthly`) — same order as `_index_transform`/
   `_synthetic_plume`. `mmm = fc.mean('model')` → `(S, L)`, then restricted
   to `S.dt.year >= config.ANALYSIS_START_YEAR` (1991). Using the **same
   fixed `avail` model set for every historical year** (rather than each
   year's actual, varying NMME roster) keeps this MMM numerically identical
   to the "MMM" line already plotted in Compare/Spread/Mean for the current
   forecast — a deliberate choice (confirmed with the user) over
   reproducing each year's true historical model set.
2. **Rank** (`_rank_at_lead(mmm_hist, start_month, current_S)`): pool =
   `mmm_hist` restricted to `S.dt.month == start_month` (the current
   forecast's init month) → `(S', L)`, `S'` spanning
   `ANALYSIS_START_YEAR`-present (pool size `n` printed/shown per row).
   `current = pool.sel(S=current_S)` → `(L,)`. `rank = (pool >
   current).sum('S') + 1` — 1 when nothing in the pool exceeds the current
   value (i.e. the current forecast is the warmest on record for that start
   month/lead); NaN comparisons contribute 0 (don't inflate the count).
   `rank` is then masked to NaN wherever `current` itself is NaN (seasonal
   endpoint leads) via `rank.where(current.notnull())` — otherwise a NaN
   current value would still get a defined (meaningless) rank since
   `pool > NaN` is elementwise `False` everywhere.
3. **Column labels**: monthly = `f"{d:%b %Y}"` for each of the `n_leads`
   `pd.date_range(start_date, periods=n_leads, freq='MS')` entries (year
   included since a mid-year start's 12 leads can cross a calendar-year
   boundary); seasonal = `_season_label(d.month)` for `leads_monthly[1:-1]`
   (the two rolling-mean endpoint leads dropped, same slicing as
   `_set_xaxis`'s seasonal tick trimming — see §4a).
4. **Row values**: anomaly formatted `f"{v:.2f}"` (`"–"` if NaN); rank
   formatted as a plain integer (`"–"` if NaN, i.e. the seasonal endpoints).
   Rank row label includes the pool size, e.g. `"n34 rank (n=36)"` — shown
   per index since `avail` (and therefore pool completeness) can differ
   between `n34` and `n34r` if a model has valid `ssta` but not
   `ssta_rel` data (or vice versa) at the current init.

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `config.N34_LAT`, `config.N34_LON` | `slice(-5,5)`, `slice(190,240)` | Standard Niño-3.4 box (5°S-5°N, 170°W-120°W) |
| `config.TWO_CLIM_GROUPS` | `{COLA-RSMAS-CCSM4, COLA-RSMAS-CESM1, NCEP-CFSv2}` | These models' hindcasts have a configuration change; a single climatology across it would bias early/late forecasts differently |
| `config.TWO_CLIM_SPLIT` | 1999 | Year the configuration change occurs; splits `TWO_CLIM_GROUPS` init times into pre/post segments |
| Two-clim, segment 1 (inits `< TWO_CLIM_SPLIT`) | climatology = per-`S.month` ensemble mean over that segment's own inits | No future-data contamination |
| Two-clim, segment 2 (inits `>= TWO_CLIM_SPLIT`) | climatology = per-`S.month` ensemble mean over inits in `[TWO_CLIM_SPLIT, CLIM_END_YEAR]` only | Fixed baseline applied to all later inits, including forecasts beyond `CLIM_END_YEAR` |
| All other models | climatology = per-`S.month` ensemble mean over **target (valid) time** in `[CLIM_START_YEAR, CLIM_END_YEAR]` | Standard fixed climatology period (1991-2020), defined on valid time rather than init time |
| `config.CLIM_START_YEAR`, `config.CLIM_END_YEAR` | 1991, 2020 | Standard WMO 30-year normal period (truncated at 2020, the last full decade at time of writing) |
| Plume color overrides | `#005030` (model 0), `#F17221` (model 1), `#1f77b4` (model 2) | Notebook convention distinguishing the first three models from the default matplotlib cycle |
| `MMM_COLOR` | `"0.75"` (light gray) | Visually distinct from all per-model colors so the multi-model mean reads as a summary line, not another model |
| `N_SYNTHETIC_MEMBERS` | 100 | Matches the paper's Fig. 9 lower-panel member count |
| `SYNTHETIC_SEED` | 0 | Fixed seed (`np.random.default_rng`) so reruns with unchanged inputs reproduce byte-identical synthetic draws |
| `SYNTHETIC_MEMBER_COLOR` | `"#7a86c8"` (light steel-blue) | Distinct from the real-member Spread palette and from `MMM_COLOR`, echoing the paper's Fig. 9 lower-panel member color |
| `SYNTHETIC_MMM_COLOR` | `"#1a1a1a"` (near-black) | `MMM_COLOR = "0.75"` is invisible against the synthetic member cloud (also light); Spread-synthetic uses a dark line for its MMM/percentile lines instead |
| Legend `ncol` (Compare/Spread/Mean) | `2` (all three) | Standardized 2026-07-09 — was inconsistent (1/2/3) before the MMM entry brought every legend to the same 8-item (7 models + MMM) count |
| Seasonal window | `rolling(L=3, center=True)` | Centered (not trailing) so the season label (e.g. DJF) matches the NOAA ONI overlapping-season convention, which is also centered |
| `SEASON_INITIALS` | `"JFMAMJJASOND"` | Single-letter month initials (index 0 = January) used to build 3-letter season labels |
| `config.TROPICS_LAT` | `slice(-20, 20)`, all longitudes | Tropical-mean region for the relative Niño-3.4 index, per L'Heureux, Tippett et al. (2024); van Oldenborgh et al. (2021) originated the index, the paper tested 15-30° alternatives and confirmed 20°S-20°N |
| Relative index definition | `ssta_rel = ssta - ssta_trop` (unscaled), then `x factor_monthly` (or `factor_seasonal`) | Niño-3.4 anomaly minus the tropical-mean anomaly better tracks local atmospheric instability / deep convection and is less sensitive to reclassification as the 30-yr climatology drifts under tropical-mean warming; subtracting the tropical mean loses variance, which the scale factor restores so fixed ±0.5°C ENSO thresholds stay meaningful |
| `config.rel_scaling_factor` formula | `factor(model, month, L) = std_obs(n34) / std_model(n34r)`, both over forecast starts `S` in 1991-2020, via `groupby('S.month')` on the `(S, L)` grid | Per the first author of L'Heureux, Tippett et al. (2024) (2026-07-06 correspondence): "scale the model relative Niño-3.4 variance to match the observed 1991-2020 Niño-3.4 variance." This is a deliberate **extension beyond the paper**, which defines only an obs-only ratio (`std_obs(n34)/std_obs(n34-trop)`, one number per calendar month, model-independent, formerly `config.load_obs_scaling`); the paper also uses 1950-2020 for its std ratio vs. 1991-2020 for the anomaly climatology — this project uses 1991-2020 for both, unchanged from the prior implementation |
| Denominator is model-dependent | `std_model(n34r)`, not `std_obs(n34-trop)` | The prior implementation's scale was a pure observational ratio (same factor applied to every model); per the first author, the factor should instead restore *each model's own* relative-index variance to the observed Niño-3.4 variance, so models with more/less variance in their `ssta_rel` get correspondingly different factors |
| Factor stratified by `(start month, L)`, not target month | `groupby('S.month')` on the `(S, L)` grid, applied via `factor.sel(month=ds.S.dt.month)` | The prior implementation keyed the factor on *target* month (`ds.target.dt.month`); the revision keys on *start* month and lets `L` vary the factor directly, since NMME model climatology/variance structure is known to vary by both init month and lead, not just by the calendar month being verified |
| Ensemble-member pooling | `sqrt(ssta_rel.groupby('S.month').var(['S', 'M']))` — grand-mean pooled std over the flattened start x member sample | Deliberately **not** `std('S')` of the ensemble mean (`ssta_rel.mean('M')`): the ensemble mean shrinks variance wherever skill is low, which would let predictability (a model-skill property) contaminate a factor meant to capture only variance (a model-climatology property). Switched 2026-07-07 from the equivalent-in-expectation per-member form (`sqrt(ssta_rel.groupby('S.month').var('S').mean('M'))`, the original choice) — see `scripts/rel_scaling_compare.py`, which found all three candidate denominators give the same anomaly correlation (as expected — correlation is scale-invariant); MSESS for the grand-mean-pooled denominator is *higher* than the ensemble-mean alternative (+0.067 averaged over model/month/lead, 1991-2020, a systematic effect of that estimator's skill contamination) and *slightly higher* than the per-member alternative too (-0.015 mean MSESS(per-member) − MSESS(grand-mean), i.e. grand-mean is marginally better in this sample — the residual of the two member-based estimators' finite-sample difference, small relative to the ensemble-mean gap, consistent with the two being equal in expectation). Grand-mean was adopted because it is not worse on this evidence and is the simpler estimator to describe ("pool all members together" vs. "average each member's own variance across members") |

## Edge Cases & Error Handling

- **Model missing its selected-init forecast** (e.g. delayed release, or a
  model not yet contributing at an older `--init-date`): excluded from all
  figures via the `avail` filter rather than plotted with NaN gaps.
- **`--init-date` matches no init in the store**: `_resolve_init_idx` raises
  `ValueError` naming the requested date and the store's full
  `start[0]:start[-1]` range, rather than silently falling back to the
  latest init.
- **`--init-date` matching month has multiple `S` entries** (shouldn't occur
  in practice — inits are one per month): `_resolve_init_idx` takes the
  first match; not otherwise guarded against.
- **Seasonal centered-window endpoints**: `rolling(L=3, center=True)` leaves
  the first lead (`L=0`) and the last lead undefined (NaN) — there is no
  full 3-month window to average there. These simply don't draw (matplotlib
  skips NaN points), so seasonal lines are one point shorter at each end
  than their monthly counterparts, and the start/end marker is placed at the
  first valid lead (`L=1`) instead of `L=0`. `_set_xaxis` additionally drops
  the corresponding two outermost ticks (see 4a) so no season-label tick is
  left dangling past where any line reaches.
- **`decode_times=False` then manual calendar patch**: the store's `calendar`
  attribute is the non-CF string `'360'`; `_decode_cf` patches it to
  `'360_day'` before calling `xr.decode_cf`. If a future store build already
  writes CF-compliant attrs, this patch becomes a no-op (safe either way).
- **Fewer than 2 init times available, or `--init-date` resolves to `now_idx
  == 0`**: `prev_idx = now_idx - 1` would be `-1`, silently wrapping to the
  *last* init rather than raising — acceptable for now since the store
  always has a long hindcast history and `--init-date` is expected to target
  a recent-ish init, not the very first one on record; not handled
  defensively.
- **Models with shorter lead count** (`NASA-GEOSS2S`, `NCEP-CFSv2`; see
  `config.N_LEADS_PLOT` docstring) are NaN-padded to 12 leads in the store,
  so no special-casing is needed in the plotting code — NaNs simply stop the
  line early.
- **A future model with incomplete 1991-2020 hindcast coverage joining
  `avail_models`**: `_synthetic_plume` raises `ValueError` naming the
  incomplete model(s) rather than silently computing the historical MMM
  error from a shrinking model set partway through the hindcast — the
  1991-2020 period was deliberately chosen because all 7 current models
  (`GFDL-SPEAR`, the newest, starts exactly 1991) have complete coverage
  there (see §4b step 1).
- **Spread-synthetic's 10th/90th percentile at seasonal endpoint leads**: the
  seasonal variant's `synthetic` array is NaN at the two rolling-mean
  endpoint leads (see above); `synthetic.quantile([0.1, 0.9],
  dim='member')` on an all-NaN slice raises `RuntimeWarning: All-NaN slice
  encountered` (from `numpy.nanpercentile` under the hood) and correctly
  returns NaN — benign, expected, confirmed by comparing total
  `RuntimeWarning` counts before/after this feature (25 baseline vs. 27
  with Spread-synthetic added, for 2 seasonal figures x 1 warning each; no
  *new* degrees-of-freedom warnings from the covariance step itself).
- **Attrs leaking into the grid colorbar label**: `ds.ssta`/`ds.ssta_rel`
  inherit a stale `standard_name`/`units` from the raw store field (copied
  wholesale by `_n34_average`/`_tropics_average`), and multiplying by the
  scaling factor in `plot_grid` can additionally pull in *that* array's
  attrs. xarray's auto-labeling prefers `long_name`/`standard_name`
  over the array's `.name`, so a leftover attr silently mislabels the
  colorbar (observed: the relative-index grid showed
  `"sea_surface_temperature"` before this was fixed). `plot_grid` now
  replaces `da.attrs` outright (`da.attrs = {"units": "degC"}`) after
  renaming to `spec['prefix']`, rather than patching the existing dict.

## Verification Snippet

```python
# Run after changes to confirm key invariants
import sys
sys.path.insert(0, "scripts")
import config
import numpy as np

ds = config.load_nino34_verification()

assert set(ds.ssta.dims) >= {"model", "S", "M", "L"}, "ssta missing expected dims"
assert str(ds.target.dtype) == "object" or "datetime" in str(ds.target.dtype), \
    "target should be decoded to cftime"
assert hasattr(ds.indexes["S"], "to_datetimeindex"), "S index should be cftime-based"

two_clim_present = config.TWO_CLIM_GROUPS & set(ds.model.values)
assert two_clim_present, "expected at least one TWO_CLIM_GROUPS model present"

assert set(ds.ssta_rel.dims) == set(ds.ssta.dims), "ssta_rel dims should match ssta"
assert np.allclose(ds.ssta_rel.values, (ds.ssta - ds.ssta_trop).values, equal_nan=True), \
    "ssta_rel should equal ssta - ssta_trop exactly (unscaled)"

factor_monthly, factor_seasonal = config.rel_scaling_factor(ds)
assert set(factor_monthly.dims) == {"model", "month", "L"}, "factor_monthly should be (model, month, L)"
assert set(factor_monthly.month.values) == set(range(1, 13)), "factor_monthly should cover all 12 months"
assert set(factor_monthly.model.values) == set(ds.model.values), "factor_monthly should cover all models"
# Finite/positive on leads every model actually has (NASA-GEOSS2S, NCEP-CFSv2
# are NaN-padded past their real lead count — see Edge Cases in specs/skill.md
# / config.N_LEADS_PLOT docstring); NOT expected to exceed 1 uniformly, unlike
# the prior obs-only factor (this is model-dependent, no such bound applies).
# .min()/.max() skipna=True by default, so NaN-padded leads don't interfere.
assert float(factor_monthly.min()) > 0, "scaling factor should be positive wherever defined"
assert np.isfinite(float(factor_monthly.max())), "scaling factor should be finite wherever defined"

sys.path.insert(0, "scripts")
from latest_forecast import _resolve_init_idx
start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
assert _resolve_init_idx(start, None) == -1, "no init_date should resolve to the latest (-1)"
resolved = _resolve_init_idx(start, "2026-06-01")
assert start[resolved].year == 2026 and start[resolved].month == 6, \
    "2026-06-01 should resolve to the S index whose init lands in 2026-06"
try:
    _resolve_init_idx(start, "1900-01-01")
    raise AssertionError("out-of-range init_date should raise ValueError")
except ValueError:
    pass

for f in [
    "n34_monthly_grid.png",
    "n34_monthly_compare.png",
    "n34_monthly_spread.png",
    "n34_monthly_spread_synthetic.png",
    "n34_monthly_mean.png",
    "n34_seasonal_compare.png",
    "n34_seasonal_spread.png",
    "n34_seasonal_spread_synthetic.png",
    "n34_seasonal_mean.png",
    "n34r_monthly_grid.png",
    "n34r_monthly_compare.png",
    "n34r_monthly_spread.png",
    "n34r_monthly_spread_synthetic.png",
    "n34r_monthly_mean.png",
    "n34r_seasonal_compare.png",
    "n34r_seasonal_spread.png",
    "n34r_seasonal_spread_synthetic.png",
    "n34r_seasonal_mean.png",
    "latest_forecast_summary.md",
]:
    assert (config.PLOTS_DIR_LATEST_FORECAST / f).exists(), f"missing output {f}"

# Summary-table rank sanity: the current forecast's own historical MMM value
# must rank first among a same-value one-element pool, and the pool must
# span ANALYSIS_START_YEAR-present.
from latest_forecast import _historical_mmm, _rank_at_lead

now_idx = _resolve_init_idx(start, None)
spec_n34 = {"var": "ssta", "prefix": "n34", "name": "Nino 3.4", "scale": None}
avail = np.where(~np.isnan(ds[spec_n34["var"]].isel(S=now_idx, L=0).mean("M")))[0]
avail_models = ds.model.isel(model=avail).values
mmm_hist = _historical_mmm(ds, avail_models, spec_n34, seasonal=False)
assert int(mmm_hist.S.dt.year.min()) >= config.ANALYSIS_START_YEAR, \
    "historical MMM pool should be restricted to ANALYSIS_START_YEAR-present"
start_month_now = int(ds.S.isel(S=now_idx).dt.month)
current_S = ds.S.isel(S=now_idx)
current, rank, pool_size = _rank_at_lead(mmm_hist, start_month_now, current_S)
assert bool((rank >= 1).where(rank.notnull(), True).all()), "rank should be >= 1 wherever defined"
assert pool_size == int((mmm_hist.S.dt.month == start_month_now).sum()), \
    "pool size should equal the count of historical starts sharing the current init month"

# Synthetic-plume sanity: same-seed rerun reproduces byte-identical draws.
from latest_forecast import (
    N_SYNTHETIC_MEMBERS, SYNTHETIC_SEED, _synthetic_plume,
)
import xarray as xr

now_idx = _resolve_init_idx(start, None)
spec_n34 = {"var": "ssta", "prefix": "n34", "name": "Nino 3.4", "scale": None}
avail = np.where(~np.isnan(ds[spec_n34["var"]].isel(S=now_idx, L=0).mean("M")))[0]
mean_list = [ds[spec_n34["var"]].isel(S=now_idx, model=im).mean("M") for im in avail]
mmm = xr.concat(mean_list, dim="model").mean("model")  # matches plot_spread_synthetic's own mmm
synth1 = _synthetic_plume(ds, avail, spec_n34, False, now_idx, mmm)
synth2 = _synthetic_plume(ds, avail, spec_n34, False, now_idx, mmm)
assert synth1.sizes["member"] == N_SYNTHETIC_MEMBERS
xr.testing.assert_identical(synth1, synth2)  # fixed seed -> identical draws

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-07-06 | Initial script + `config.load_nino34_ssta()` written | ✓ |
| 2026-07-06 | Dropped duplicate `.pdf` output for the compare figure (PNG only) | ✓ |
| 2026-07-06 | Added seasonal (3-month running mean) variants of compare/spread/mean; renamed all outputs to `n34_{monthly,seasonal}_*` | ✓ |
| 2026-07-06 | Added relative Niño-3.4 (`n34r_*`): `config.TROPICS_LAT`, `config._tropics_average`, `config._forecast_anomaly` (refactored from the inline two-clim logic), `ds.ssta_trop`/`ds.ssta_rel` in `load_nino34_ssta()`, new `config.load_obs_scaling()` (ERSSTv5, 1991-2020); `latest_forecast.py` plot functions parametrized over an index spec; fixed a colorbar-attrs-leak bug in `plot_grid` (see Edge Cases) | ✓ |
| 2026-07-06 | Added disk caching to `config.load_nino34_ssta()`/`config.load_obs_scaling()` (`use_cache=True` default), `config.CACHE_DIR`, `config._store_fingerprint()`; ~15x speedup on repeat runs (~60s → ~4s) | ✓ |
| 2026-07-06 | Fixed color bug: `_plume_colors()` now returns one array (the U Miami-override array) used by Compare, Spread, and Mean alike — previously only Compare got the overridden colors, Spread/Mean fell back to the plain default cycle | ✓ |
| 2026-07-06 | Reworked x-axis tick/xlim discipline (`_tight_xlim`, `_set_xaxis` seasonal trimming): each plot computes its own ticks/xlim instead of Spread reusing Compare's (wider) xlim; Mean's tick array no longer has one extra tick past its last data point; seasonal variant drops the two outermost (NaN-endpoint) ticks. See Algorithm §4a | ✓ |
| 2026-07-06 | Added `--init-date YYYY-MM[-DD]` CLI argument (`_resolve_init_idx`, `now_idx`/`prev_idx`/`date_suffix` threaded through all plot functions in place of hardcoded `S=-1`/`S=-2`); output filenames get a `_<YYYY-MM>` suffix when used. Tested against `--init-date 2026-06-01`. See CLI Arguments, Algorithm §2 | ✓ |
| 2026-07-06 | `config.load_obs_scaling()`'s ERSSTv5 box-mean/anomaly logic factored out into a new shared helper `config._obs_index_anomalies()` (also used by the new `config.load_nino34_verification()`, see `specs/skill.md`). Behavior-preserving refactor only — `load_obs_scaling`'s public signature and return values are unchanged; verified bit-identical (`np.allclose`) against a pre-refactor fresh (`use_cache=False`) recompute. Not itself a change to this script's algorithm. | ✓ (noted only) |
| 2026-07-06 | **Revised the relative Niño-3.4 scaling factor formula.** Per the first author of L'Heureux, Tippett et al. (2024): "scale the model relative Niño-3.4 variance to match the observed 1991-2020 Niño-3.4 variance" — the denominator must be the *model's* `ssta_rel` std, not obs's, and it should depend on `(start month, lead)`, not target month alone. Removed `config.load_obs_scaling()` (and its cache files `cache/obs_scaling_1991_2020.{nc,fingerprint.json}`); added `config.rel_scaling_factor(ds)`, which reads `ds.obsa`/`ds.ssta_rel` from `config.load_nino34_verification()` and computes `factor(model, month, L) = std_obs(n34) / std_model(n34r)` via `groupby('S.month')` over 1991-2020 starts, with the denominator pooling ensemble members as `sqrt(var('S').mean('M'))` (not the ensemble mean — see Constants). `main()` now calls `load_nino34_verification()` instead of `load_nino34_ssta()`; `_index_transform` now takes `(model, start_month)` instead of `target` and selects `factor.sel(model=..., month=...)`; `plot_grid` likewise keys on the selected init's start month. Added `scripts/rel_scaling_compare.py` (exploratory, non-production) contrasting the chosen member-pooled denominator against the (flawed) ensemble-mean alternative: anomaly correlation is identical (max diff ~1e-15, as expected — correlation is scale-invariant), and MSESS is higher for the chosen estimator (mean +0.052 across model/month/lead, 1991-2020), supporting the choice. See Algorithm §1a and Constants & Scientific Rationale. | ✓ |
| 2026-07-06 | **Documented the second member-based pooling estimator and added its comparison.** The per-member denominator (`sqrt(var('S').mean('M'))`) is one of two valid member-based estimators, equal in expectation but differing in finite samples — the other is grand-mean pooling (`sqrt(var(['S','M']))`, the pooled std over the flattened start x member sample), which additionally carries the between-member spread of the per-member time-means. No change to `config.rel_scaling_factor`'s formula — docstring only. `scripts/rel_scaling_compare.py` now also computes this grand-mean factor (C, local `_factor_grandmean`) alongside the existing ensemble-mean alternative (B): AC identical to A (max diff ~1e-15), mean MSESS(A) − MSESS(C) = -0.015 across model/month/lead (1991-2020) — small relative to the +0.052 A-vs-B gap, consistent with A and C being equal in expectation. New plots `n34r_scaling_msess_C_{start,target}.png`, `n34r_scaling_msess_diff_AC_{start,target}.png`; existing A-vs-B diff plots renamed `n34r_scaling_msess_diff_AB_{start,target}.png` (from unsuffixed `_diff_{start,target}.png`) to disambiguate. See Algorithm §1a and Constants & Scientific Rationale. | ✓ |
| 2026-07-07 | **Switched `config.rel_scaling_factor`'s denominator from per-member to grand-mean pooling.** Formula change: `sqrt(ref.ssta_rel.groupby('S.month').var('S').mean('M'))` → `sqrt(ref.ssta_rel.groupby('S.month').var(['S', 'M']))` (both `factor_monthly` and `factor_seasonal`). Motivated by `scripts/rel_scaling_compare.py`'s A-vs-C evidence (`n34r_scaling_msess_diff_AC_start.png`): MSESS is not worse for grand-mean (mean MSESS(per-member) − MSESS(grand-mean) = −0.015 across model/month/lead, i.e. grand-mean marginally better), and grand-mean is the simpler estimator to describe. `rel_scaling_compare.py` restructured accordingly: `A` = per-member (now the local alternative, `_factor_permember`, moved out of `config.py`), `B` = ensemble-mean (flawed, unchanged), `C` = `config.rel_scaling_factor` (now grand-mean, chosen/production). Added the third pairwise comparison, B-vs-C (`n34r_scaling_msess_diff_BC_{start,target}.png`, full-range color scale — not dominated by outliers the way A-vs-C was): mean MSESS(B) − MSESS(C) = −0.067, i.e. grand-mean beats the flawed ensemble-mean estimator by even more than per-member did. All production `n34r_*`/`n34r_seasonal_*` figures regenerated (`latest_forecast.py`) — scaling factor numeric range shifted slightly, 1991-2020 monthly range [0.60, 2.42] → [0.51, 2.30]. See Algorithm §1a and Constants & Scientific Rationale. | ✓ |
| 2026-07-07 | Output moved from `plots/` to `plots/latest_forecast/` (per-script subdirectory, `config.PLOTS_DIR_LATEST_FORECAST`) | ✓ |
| 2026-07-07 | `config.ERSSTV5_NC` default changed to the repo-local `OBS_DIR / "ERSSTv5.sst.mnmean.nc"` (see `specs/skill.md` same-date row for details). `rel_scaling_factor` 1991-2020 range verified unchanged ([0.51, 2.30] monthly); figures regenerated. | ✓ |
| 2026-07-22 | **Hardened the synthetic family's hindcast-period assumption.** Confirmed empirically (360/360 starts per model at zero lead) that all 7 current NMME models have complete 1991-2020 coverage, which is *why* that period was chosen for the error covariance, not just because it's the standard climatology window. Added an explicit check in `_synthetic_plume` that raises `ValueError` if any `avail_models` member has incomplete zero-lead coverage over that period, rather than letting `mean('model')`'s skipna behavior silently shrink the effective model set for part of the hindcast. No change to computed output (all models currently pass). See Algorithm §4b and Edge Cases. | ✓ |
| 2026-07-22 | **Added the synthetic error-covariance spread family** (`plot_spread_synthetic`, `_synthetic_plume`), implementing Barnston, Tippett, van den Dool & Unger (2015, *J. Appl. Meteor. Climatol.*, **54**, 1579–1595, https://doi.org/10.1175/JAMC-D-14-0188.1, Fig. 9 lower panels): 100 Gaussian scenarios drawn from the historical (1991-2020) MMM forecast-error covariance across leads, stratified by the current start month, added to the current MMM. Both n34 and n34r verify against observed absolute Niño-3.4 (`ds.obsa`), per direct author request — not `ds.obsa_rel`. New constants `N_SYNTHETIC_MEMBERS=100`, `SYNTHETIC_SEED=0`, `SYNTHETIC_MEMBER_COLOR`, `SYNTHETIC_MMM_COLOR`. New outputs `n34{,r}_{monthly,seasonal}_spread_synthetic.png` (4 new files; totals 7→9 per index, 14→18 overall). No `config.py` changes — reuses `load_nino34_verification()`'s existing `ssta`/`ssta_rel`/`obsa` and `rel_scaling_factor`. Verified: all 4 new figures render with the same MMM as the corresponding `*_spread.png`; RuntimeWarning count 25→27 (2 new, both benign all-NaN-slice from the seasonal quantile's NaN endpoint leads — no new degrees-of-freedom warnings from the covariance step). See Algorithm §4b, Constants, Edge Cases. | ✓ |
| 2026-07-29 | **Added the summary-table output** (`write_summary_tables`, `_historical_mmm`, `_rank_at_lead`): `plots/latest_forecast/latest_forecast_summary.md`, monthly + seasonal MMM anomaly tables for both indices, each value's rank (1 = highest) among all MMM forecasts issued in the same calendar start month, `ANALYSIS_START_YEAR`-present. Design decisions confirmed with the user: one combined table per kind (4 rows: `n34`/`n34r` anom/rank) rather than 4 separate tables; ranking pool uses the **fixed model set from the current forecast** applied across all historical years (matching the plotted MMM line), not each year's true historical roster. See Algorithm §5. | ✓ |
| 2026-07-09 | **Added a multi-model mean (MMM) line to Compare/Spread/Mean.** New `MMM_COLOR = "0.75"` constant; each function collects the per-model arrays it already computes into a list during the loop, then plots `xr.concat(..., dim="model").mean("model")` after the loop with no explicit `zorder` (renders on top, drawn last) and `label="MMM"`. Also standardized `ax.legend(ncol=...)` to `2` in all three (was 1/2/3). **Also tried, then reverted same-session:** adding an 8th "MMM" panel to the Grid facet (`plot_grid`) by broadcasting the same MMM series across the `M` coordinate into a uniform-color block — the user judged this a bad idea and asked it removed; Grid stays at 7 panels, 8th `col_wrap` slot empty, as before. See Algorithm §4 and Constants. All 12 line-plot figures (`n34_*`/`n34r_*` compare/spread/mean, monthly/seasonal) regenerated; Grid figures unchanged from pre-session. | ✓ |
