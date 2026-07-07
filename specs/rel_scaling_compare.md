# rel_scaling_compare.py — Behavioral Specification

> Last reviewed against code: 2026-07-07 (production switched to grand-mean pooling; added B-C comparison)

## Purpose

Provide data-based evidence for the member-pooling choice in
`config.rel_scaling_factor()` (see `specs/latest_forecast.md` §1a): the
relative Niño-3.4 scaling factor's denominator pools ensemble members via
the grand mean — the pooled std over the flattened (start, member) sample,
`sqrt(var(['S', 'M']))` (factor **C, chosen production formula since
2026-07-07**). This script provides evidence on all three pairwise
comparisons among the candidate denominators:

- **A vs B (per-member vs. ensemble-mean, flawed)**: A is the per-member
  temporal variance averaged across members (`sqrt(var('S').mean('M'))`,
  the original production choice, now the local alternative). B is the
  variance of the ensemble mean (`std('S')` of `ssta_rel.mean('M')`). B is
  flawed — it shrinks variance wherever forecast skill is low, mixing
  predictability into a quantity meant to capture only variance.
- **A vs C (per-member vs. grand-mean, equal in expectation)**: A and C are
  equal in expectation but differ in finite samples — C additionally
  carries the between-member spread of the per-member time-means (law of
  total variance), which A's per-member-then-average pooling does not.
  Neither is flawed by construction; this comparison quantifies how much
  the finite-sample distinction actually matters. **It motivated the
  2026-07-07 switch** (`n34r_scaling_msess_diff_AC_start.png`): MSESS is
  not worse for C — slightly better on average — and C is the simpler
  estimator to describe.
- **B vs C (ensemble-mean vs. grand-mean, flawed vs. chosen)**: added for
  completeness alongside A vs C, so the comparison against the flawed
  estimator (B) is available for both member-based candidates (A and C),
  not just the one that used to be production.

For all three axes, the script scales the same forecasts with each factor
and compares anomaly correlation (AC) and MSESS against a common
observational target, to make the practical impact of each choice visible.
Exploratory/diagnostic only — does not feed production figures or the
manuscript directly, though its A-vs-C evidence informed the
`config.rel_scaling_factor` formula change.

## Inputs

| File | Relevant columns | Filters applied |
|------|-------------------|-----------------|
| `$NMME_STORE_DIR/nmme_sst.zarr` (`config.STORE_SST`) | `sst(S, M, L, Y, X)`, all 7 model groups | Via `config.load_nino34_verification()` (which wraps `load_nino34_ssta()`); restricted to forecast starts `S` in `[CLIM_START_YEAR, CLIM_END_YEAR]` (1991-2020) |
| `config.ERSSTV5_NC` (ERSSTv5 monthly SST) | `sst(time, lat, lon)` | Via `load_nino34_verification()`'s `obsa`/`obsa_rel` (aligned onto the `(S, L)` grid), same 1991-2020 restriction |

### Caching

Uses `config.load_nino34_verification()` → `load_nino34_ssta()`'s existing
disk cache (`cache/nino34_ssta.nc`, see `specs/latest_forecast.md`). No
caching of its own — both scaling-factor variants and all skill scores are
cheap reductions over the small, already-cached verification Dataset.

## Outputs

| File | Contents | Format |
|------|----------|--------|
| `plots/rel_scaling_compare/n34r_scaling_msess_A_start.png` | MSESS (per-member-pooled denominator, alternative), by start month x lead, 8-panel (7 models + MMM) | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_A_target.png` | Same, relabeled to target month | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_B_start.png` | MSESS (ensemble-mean denominator, flawed), by start month x lead | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_B_target.png` | Same, relabeled to target month | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_C_start.png` | MSESS (grand-mean pooled denominator, **chosen/production**), by start month x lead | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_C_target.png` | Same, relabeled to target month | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_AB_start.png` | MSESS(A) − MSESS(B) (per-member scaling − ens. mean scaling), by start month x lead, full data-range diverging color scale | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_AB_target.png` | Same, relabeled to target month | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_AC_start.png` | MSESS(A) − MSESS(C) (per-member scaling − grand-mean scaling), by start month x lead, color scale clipped to ±0.03 (`cmap="RdBu_r"`, explicit `vmin`/`vmax`) so the small-magnitude signal in 6 of 8 panels is visible against SPEAR/GEOSS2S's much larger outlier magnitudes (which saturate/clip) | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_AC_target.png` | Same, relabeled to target month | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_BC_start.png` | MSESS(B) − MSESS(C) (ens. mean scaling − grand-mean scaling), by start month x lead, full data-range diverging color scale (not dominated by outliers the way A-C is — no clipping needed) | PNG, dpi=150 |
| `plots/rel_scaling_compare/n34r_scaling_msess_diff_BC_target.png` | Same, relabeled to target month | PNG, dpi=150 |

Console output also prints the max absolute AC difference for all three
pairs (A-vs-B, A-vs-C, B-vs-C — sanity checks, not plotted, since AC is not
saved to a figure), and the mean MSESS(A) − MSESS(B), MSESS(A) − MSESS(C),
and MSESS(B) − MSESS(C) across model/month/lead.

## Algorithm

### 1. Load and restrict

`ds = config.load_nino34_verification()`, then `ds_p = ds.sel(S=slice(period[0]-01-01, period[1]-12-01))`
with `period = (config.CLIM_START_YEAR, config.CLIM_END_YEAR)` (1991-2020) —
identical restriction pattern to `skill.py`'s `main()`.

### 2. Three scaling factors

- **A (per-member, alternative)**: `_factor_permember(ds, period)`, defined
  in this script (moved out of `config.py` on 2026-07-07 — it was the
  production formula until then; retained here purely as the comparison
  baseline that motivated the switch):
  ```
  num = ref.obsa.groupby('S.month').std('S')                          # (month, L)
  den = np.sqrt(ref.ssta_rel.groupby('S.month').var('S').mean('M'))   # (model, month, L)
  factor_A = num / den
  ```
- **B (ensemble-mean, flawed alternative)**: `_factor_ensmean(ds, period)`,
  defined in this script (deliberately kept out of `config.py` — it exists
  only to demonstrate why the member-based estimators, A and C, are
  preferred):
  ```
  num = ref.obsa.groupby('S.month').std('S')                       # (month, L)
  den = ref.ssta_rel.mean('M').groupby('S.month').std('S')          # (model, month, L)
  factor_B = num / den
  ```
- **C (grand-mean, chosen)**: `factor_C, _ = config.rel_scaling_factor(ds, period=period)`
  — the production factor as of 2026-07-07 (grand-mean-pooled denominator;
  see `specs/latest_forecast.md` §1a for the full formula).

  where `ref` is `ds` restricted to `period` via `ds.where(year in period, drop=True)`
  — same restriction logic as `config.rel_scaling_factor`, for both A and B.
  Only the monthly (not seasonal) variant is computed for A and B, since only
  the monthly forecast is compared here.

### 3. Scaled forecasts

`_scaled_forecast(ds_p, factor)` computes the ensemble-mean scaled
relative-index forecast for a given factor:
```
f = factor.sel(month=ds_p.S.dt.month)   # pointwise/vectorized sel by each start's calendar month → (model, S, L)
forecast = (ds_p.ssta_rel * f).mean('M')  # (model, S, L)
```
Applied with `factor_A` → `fcst_A`, `factor_B` → `fcst_B`, `factor_C` → `fcst_C`.

### 4. Verification target

The observed relative index scaled to match observed Niño-3.4 variance (the
paper's own obs-only RONI, identical for A, B, and C — see
`specs/latest_forecast.md` for why this is *not* the same formula as
`config.rel_scaling_factor`, which uses a model-relative denominator):
```
obs_scale = ds_p.obsa.groupby('S.month').std('S') / ds_p.obsa_rel.groupby('S.month').std('S')
obs_target = ds_p.obsa_rel * obs_scale.sel(month=ds_p.S.dt.month)
```

### 5. Skill comparison

Reuses `skill.py`'s own functions (`ac_by_start`, `msess_by_start`,
`add_mmm`, `to_target_month`, `plot_skill_grid`) rather than reimplementing
them, so the comparison's skill-score definitions are identical to the
production skill figures:

- `ac_A`, `ac_B`, `ac_C = ac_by_start(fcst_*, obs_target)` — **not** plotted;
  only `max(abs(ac_A - ac_B))`, `max(abs(ac_A - ac_C))`, and
  `max(abs(ac_B - ac_C))` are printed, as a sanity check that AC is
  invariant to the per-(model, month, L) scalar factor within a start-month
  group (it is a linear rescaling, and Pearson correlation is invariant to
  affine transforms of either input).
- `msess_A = msess_by_start(fcst_A, obs_target)`, then `add_mmm(...)` to
  append the multi-model-mean panel — same for `msess_B` and `msess_C`.
- `to_target_month()` relabels each to target month (`msess_A_target`,
  `msess_B_target`, `msess_C_target`).
- `msess_diff_ab = msess_A - msess_B` and `msess_diff_ab_target` — MSESS *is*
  not scale-invariant (unlike AC), so this difference is evidence that B is
  worse than A. `msess_diff_ac = msess_A - msess_C` and
  `msess_diff_ac_target` are the analogous evidence for how close C is to A
  (both being valid member-based estimators, equal in expectation) — this
  pair motivated the 2026-07-07 production switch. `msess_diff_bc =
  msess_B - msess_C` and `msess_diff_bc_target` complete the third pairing,
  evidence that B is worse than C too (added for completeness alongside
  A-vs-C).
- Plot titles for the diff figures are descriptive, not `MSESS(A) - MSESS(B)`
  formula notation — `"n34r MSESS, per-member scaling - ens. mean scaling"`
  for the A-B diff, `"n34r MSESS, per-member scaling - grand-mean scaling"`
  for the A-C diff, `"n34r MSESS, ens. mean scaling - grand-mean scaling"`
  for the B-C diff — identical between the start- and target-month framing
  of each (the framing itself is conveyed by the x-axis label, "start
  month" vs "target month").
- The A-B and B-C diffs (`msess_diff_ab*`, `msess_diff_bc*`) plot with
  `plot_skill_grid`'s default color scale (data-range-derived, diverging) —
  neither is dominated by outlier models enough to need clipping (verified
  visually before deciding not to clip B-C). The A-C diff
  (`msess_diff_ac*`) plots with an explicit clipped range,
  `ac_diff_vlim = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")`: SPEAR and
  GEOSS2S have much larger-magnitude values here (min ≈ −0.34, −0.40) than
  the other 5 models/MMM (≈ ±0.01-0.03), so the default data-range color
  scale would saturate all detail in those 6 panels. `cmap="RdBu_r"` must be
  passed explicitly alongside the clipped `vmin`/`vmax` — passing only
  `vmin`/`vmax` (even with `center=0`) makes xarray fall back to a
  sequential colormap; see Edge Cases.

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `period` | `(config.CLIM_START_YEAR, config.CLIM_END_YEAR)` = 1991-2020 | Same restriction as `config.rel_scaling_factor` and `skill.py`, so the comparison uses the same sample all three factors were fit on |
| Verification target | Observed RONI (`obsa_rel` scaled by the paper's own obs-only ratio) | A model-independent, denominator-choice-independent ground truth — using `ssta_rel` scaled by any of factor_A/B/C as its own target would trivially favor whichever factor produced it |
| `_factor_permember` kept out of `config.py` | — | It was the production formula before 2026-07-07; retained here as the comparison baseline that motivated the switch to grand-mean, not as a supported option going forward |
| `_factor_ensmean` kept out of `config.py` | — | It is a deliberately flawed alternative, used only for this comparison; keeping it here (not in `config.py`) avoids presenting it as a supported option |
| `ac_diff_vlim = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")` | Clip range for the A-C diff plots only (not A-B or B-C) | SPEAR/GEOSS2S's MSESS(A)-MSESS(C) magnitude (~10x the other models') would otherwise saturate the whole color scale via xarray's default data-range-derived limits, hiding the signal this comparison exists to show — see Algorithm §5. The A-B and B-C diffs don't have this problem (checked visually), so they keep the default data-range scale |
| Production switch (2026-07-07) | `config.rel_scaling_factor` denominator: per-member → grand-mean | Per this script's A-vs-C evidence: MSESS not worse for grand-mean (mean MSESS(A) − MSESS(C) = −0.015, i.e. grand-mean marginally better), and grand-mean is the simpler estimator to describe. See `specs/latest_forecast.md` §1a/Constants and its Synchronization Log for the full rationale |

## Edge Cases & Error Handling

- **Models with shorter lead count** (`NASA-GEOSS2S`, `NCEP-CFSv2`; see
  `config.N_LEADS_PLOT` docstring): NaN-padded past their real lead count in
  the store, so `msess_by_start`/`ac_by_start`'s masking (via
  `x.where(ok)`) naturally excludes those cells — no special-casing needed.
- **AC "identical" is only exact up to floating-point roundoff**: observed
  max `|AC_A - AC_B|` ~ 1e-15, max `|AC_A - AC_C|` ~ 1e-15, max
  `|AC_B - AC_C|` ~ 1e-15 (2026-07-07 run) — printed, not asserted, since
  this is a diagnostic script rather than a test.
- **MSESS(A) − MSESS(C) is not expected to be ~0** despite A and C being
  equal in expectation: they are fit on the same finite 1991-2020 sample, so
  their difference reflects genuine (if small) finite-sample estimator
  variance, not a bug — observed mean −0.015 across model/month/lead
  (2026-07-07 run), versus +0.052 for MSESS(A) − MSESS(B) and −0.067 for
  MSESS(B) − MSESS(C) (both the systematically flawed B being worse than
  either member-based estimator, A or C).
- **`skill.py`'s `plot_skill_grid(..., vmin=, vmax=)` silently drops the
  diverging colormap when both are given explicitly** — verified empirically
  (xarray, this project's pinned version): `center=0` alone does *not*
  restore it; only an explicit `cmap="RdBu_r"` does. `plot_skill_grid`
  therefore takes a `cmap` passthrough parameter (default `None`, i.e. no
  behavior change for existing full-data-range callers); any future
  clipped-range diff plot must pass `cmap="RdBu_r"` explicitly alongside
  `vmin`/`vmax`, or it will render in the sequential default (viridis)
  instead — this was caught visually, not by an assertion, since it's a
  cosmetic-not-numeric failure mode.
- **SPEAR/GEOSS2S render fully saturated (solid clipped color) in the A-C
  diff plots**, not their true extreme values — by design, since the clip
  exists specifically to make the other 6 panels legible. Their actual
  magnitude is available from the unclipped `msess_diff_ac`/`_target`
  DataArrays (not just the plots) and from the printed
  `mean MSESS(A) - MSESS(C)` console line.

## Verification Snippet

```python
# Run after changes to confirm key invariants
import sys
sys.path.insert(0, "scripts")
import config
import numpy as np
from rel_scaling_compare import _factor_permember, _factor_ensmean, _scaled_forecast

ds = config.load_nino34_verification()
period = (config.CLIM_START_YEAR, config.CLIM_END_YEAR)

factor_A = _factor_permember(ds, period=period)
factor_B = _factor_ensmean(ds, period=period)
factor_C, _ = config.rel_scaling_factor(ds, period=period)
for f in (factor_A, factor_B, factor_C):
    assert set(f.dims) == {"model", "month", "L"}, \
        "all three factors should be (model, month, L)"
# A and C should be close on average (equal in expectation) — a coarse
# sanity check that the roles haven't been mixed up. Use the mean, not max,
# absolute difference: a handful of sparse-data cells (e.g. GEOSS2S's short
# lead count) can produce large single-cell differences even though A and C
# are close everywhere else (observed mean ~0.03, max ~0.7, 2026-07-07 run).
assert float(np.abs(factor_A - factor_C).mean()) < 0.2, \
    "per-member (A) and grand-mean (C) factors should be close on average"

period_sel = f"{period[0]}-01-01", f"{period[1]}-12-01"
ds_p = ds.sel(S=slice(*period_sel))
fcst_A = _scaled_forecast(ds_p, factor_A)
fcst_B = _scaled_forecast(ds_p, factor_B)
fcst_C = _scaled_forecast(ds_p, factor_C)
for f in (fcst_A, fcst_B, fcst_C):
    assert set(f.dims) == {"model", "S", "L"}, \
        "all three scaled forecasts should be (model, S, L)"

for f in [
    "n34r_scaling_msess_A_start.png",
    "n34r_scaling_msess_A_target.png",
    "n34r_scaling_msess_B_start.png",
    "n34r_scaling_msess_B_target.png",
    "n34r_scaling_msess_C_start.png",
    "n34r_scaling_msess_C_target.png",
    "n34r_scaling_msess_diff_AB_start.png",
    "n34r_scaling_msess_diff_AB_target.png",
    "n34r_scaling_msess_diff_AC_start.png",
    "n34r_scaling_msess_diff_AC_target.png",
    "n34r_scaling_msess_diff_BC_start.png",
    "n34r_scaling_msess_diff_BC_target.png",
]:
    assert (config.PLOTS_DIR_REL_SCALING_COMPARE / f).exists(), f"missing output {f}"

print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-07-06 | Initial `scripts/rel_scaling_compare.py` written, comparing `config.rel_scaling_factor`'s member-pooled denominator against an ensemble-mean alternative (via AC/MSESS, 1991-2020). Result: AC identical (max diff ~1e-15), MSESS higher for the member-pooled factor (mean +0.05 across model/month/lead) — supports the choice documented in `config.rel_scaling_factor`'s docstring and `specs/latest_forecast.md` §1a/Constants. | ✓ |
| 2026-07-06 | Added a second comparison axis: `_factor_grandmean` (factor C, grand-mean pooling `sqrt(var(['S','M']))`) alongside the existing ensemble-mean alternative (B). C is equal in expectation to the chosen per-member estimator (A) but differs in finite samples. New outputs `n34r_scaling_msess_C_{start,target}.png`, `n34r_scaling_msess_diff_AC_{start,target}.png`; renamed the existing A-vs-B diff outputs from unsuffixed `_diff_{start,target}.png` to `_diff_AB_{start,target}.png` to disambiguate. Result (1991-2020 run): AC identical to A (max diff ~1e-15), mean MSESS(A) − MSESS(C) = −0.015 across model/month/lead — small relative to the +0.052 A-vs-B gap, consistent with A and C being equal in expectation. Also expanded the two-member-method note in `config.rel_scaling_factor`'s docstring (no formula change) and `specs/latest_forecast.md` §1a/Constants. | ✓ |
| 2026-07-07 | Cosmetic-only follow-up, no score/formula changes. (1) A-B/A-C diff plot titles: `_start` variants now use the same descriptive text as `_target` (`"n34r MSESS, member scaling - ens. mean scaling"` / `"...per-member scaling - grand-mean scaling"`) instead of `MSESS(A) - MSESS(B)`/`MSESS(A) - MSESS(C)` formula notation. (2) A-C diff plots (`msess_diff_AC_{start,target}.png`) now use a clipped color scale, `vmin=-0.03, vmax=0.03, cmap="RdBu_r"`, since SPEAR/GEOSS2S's much larger magnitude there (min ≈ −0.34, −0.40) was saturating the default data-range scale and hiding the ≈ ±0.01-0.03 signal in the other 6 panels; `skill.py`'s `plot_skill_grid` gained a `cmap` passthrough param for this (default `None`, no behavior change elsewhere). Discovered along the way: passing explicit `vmin`/`vmax` to `plot_skill_grid` disables xarray's automatic diverging-colormap detection even with `center=0` — `cmap="RdBu_r"` must be passed explicitly (see Edge Cases). | ✓ |
| 2026-07-07 | **Switched production to grand-mean pooling; restructured A/B/C.** Based on the A-vs-C evidence above (MSESS not worse for grand-mean, simpler to describe), `config.rel_scaling_factor`'s denominator formula changed from per-member to grand-mean pooling (see `specs/latest_forecast.md` Synchronization Log for the formula diff). This script's roles flipped to match: **A** = per-member (`_factor_permember`, new local function — moved out of `config.py`, was `config.rel_scaling_factor` before this change), **B** = ensemble-mean (`_factor_ensmean`, unchanged), **C** = `config.rel_scaling_factor` (now grand-mean, chosen/production — was the local `_factor_grandmean`, now removed since it duplicates the production formula). Added the third pairwise comparison, B-vs-C, for completeness alongside A-vs-B and A-vs-C: `msess_diff_BC_{start,target}.png` (full data-range color scale, not clipped — checked visually, not dominated by outliers the way A-C is). Result (1991-2020 run): mean MSESS(B) − MSESS(C) = −0.067 (grand-mean beats the flawed ensemble-mean estimator by even more than per-member did, +0.052 for A-vs-B). All production `n34r_*` figures regenerated via `latest_forecast.py` (scaling factor range shifted [0.60, 2.42] → [0.51, 2.30]). | ✓ |
| 2026-07-07 | Output moved from `plots/` to `plots/rel_scaling_compare/` (per-script subdirectory, `config.PLOTS_DIR_REL_SCALING_COMPARE`) | ✓ |
| 2026-07-07 | `config.ERSSTV5_NC` default changed to the repo-local `OBS_DIR / "ERSSTv5.sst.mnmean.nc"` (see `specs/skill.md` same-date row for details). All three mean MSESS pairwise differences verified bit-stable (+0.0521 / −0.0149 / −0.0670); figures regenerated. | ✓ |
