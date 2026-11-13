# forecast_flashlight.py — Behavioral Specification

> Last reviewed against code: 2026-11-13

## Purpose

An animated monitoring figure: does each month's NMME relative Niño-3.4
(n34r) seasonal forecast follow the observed path? The user's framing: "the
obs are a fixed path and the forecast each month is the beam of a
flashlight, hopefully following the path." Axes are fixed in calendar time
and value across frames, which is the deliberate departure from the
reference it was modeled on (a user-supplied NMME ENSO Outlook RONI GIF whose
x-axis is labeled relative to each issue date, so the observed line moves
between frames). Exploratory/monitoring, not a manuscript figure.

## Inputs

| File | Relevant variables | Filters applied |
|------|--------------------|-----------------|
| `config.STORE_SST` via `config.load_nino34_verification()` | `ssta_rel(model, S, M, L)`, `obsa(S, L)` | Inits `--first-init` … `--last-init` (default 2025-11 … 2026-11); per init, models with non-NaN lead-0 forecast (`latest_forecast._available_models`) |
| `config.ERSSTV5_NC` via `config._obs_index_anomalies()` | `n34_anom(time)`, `rel_anom(time)` (1991-2020 monthly climatology) | Plotted only within the fixed x-window |
| `config.rel_scaling_factor(ds)` | seasonal factor `(model, month, L)` | — |

## Outputs

| File | Contents | Format |
|------|----------|--------|
| `plots/forecast_flashlight/n34r_seasonal_<first>_<last>.gif` | One frame per init: current beam (MMM + 5–95/15–85/25–75% bands), fixed observed path, dotted line at issue month, "Issued <Mon YYYY>" box | GIF, 1300×700 (13×7 in, dpi 100), 1 s/frame, last frame held 3 s, loops |
| `plots/forecast_flashlight/n34r_seasonal_<first>_<last>_trail.gif` (`--trail`) | As above plus every earlier init's MMM line and 5–95% band in faded gray, beneath the current beam | GIF, same |

## Algorithm

1. **Load.** `ds = config.load_nino34_verification()`; resolve the first/last
   init to positional S indices with `latest_forecast._resolve_init_idx`
   (calendar-month match); build the n34r spec with
   `config.rel_scaling_factor(ds)`.
2. **Beam, per init** (`beam()`):
   - MMM: for each available model, `latest_forecast._index_transform(...,
     seasonal=True)` (3-month centered running mean over L, then the
     model/start-month/lead n34r factor), mean over M, then mean over models.
     Identical to the MMM in `latest_forecast`'s n34r seasonal figures.
   - Spread: `valid_L, cov = latest_forecast._lead_error_cov(ds, avail,
     spec, True, now_idx)` — 1991-2020 MMM forecast-error covariance across
     leads, stratified by the init's start month, errors against observed
     **absolute** Niño-3.4 (`ds.obsa`), per that function's convention.
     `sd_L = sqrt(diag(cov))`.
   - Bands: for each q ∈ {0.05, 0.15, 0.25, 0.75, 0.85, 0.95},
     `q_L = MMM_L + Φ⁻¹(q)·sd_L` — the closed-form marginals of the synthetic
     plume (Barnston et al. 2015), as in `plot_strength_probabilities`, not
     Monte Carlo.
   - Target time of lead L = init month + L months (center month of the
     3-month season). Only `valid_L` (L = 1…10) is kept. NMME `L` labels are
     0.5, 1.5, …, so the whole-month offset is `int(L)`.
   - Returned Dataset: `mean`, `sd` (per-lead hindcast-error sd, used by
     `forecast_flashlight_video.py`), and the `q<pct>` band bounds.
3. **Observed path** (`obs_n34r_seasonal()`): 3-month centered running mean
   of ERSSTv5 `n34_anom` and `rel_anom`; scale per center calendar month
   m: `s_m = std(n34_anom_seas | m) / std(rel_anom_seas | m)` over
   1991-2020; path = `rel_anom_seas · s_m`. This is the obs-only RONI-style
   scaling (L'Heureux et al. 2024), the same target as
   `rel_scaling_compare.py` §4 but computed on the seasonal series.
4. **Fixed axes.** x: one month before the first beam's first target through
   the last beam's last target, monthly ticks labeled by season (year shown
   on the first tick and on Jan/Jul). y: min/max over all beams' 5%/95%
   bounds and the obs in the x-window, floored/ceiled to 0.5 °C.
5. **Frames.** Each frame clears and redraws: (`--trail` only) gray 5–95%
   band + gray MMM line for every earlier init; the current beam's three
   bands (outermost first) and MMM; the observed path; a dotted vertical
   line at the issue month; the issue-date box. Last frame repeated
   `LAST_FRAME_HOLD` times (Pillow merges duplicates into one longer frame).

## Constants & Scientific Rationale

| Name | Value | Why |
|------|-------|-----|
| `BAND_QUANTILES` | (0.05,0.95), (0.15,0.85), (0.25,0.75) | Match the reference NMME ENSO Outlook GIF's 90/70/50% ranges |
| Spread method | Gaussian MMM ± hindcast-error sd | User's choice (2026-10-19) over pooled-member quantiles: the calibrated synthetic spread already used in `latest_forecast` |
| Index | n34r seasonal (ERSSTv5 obs) | User's choice (2026-10-19) over CPC RONI obs or standard Niño-3.4: keeps obs and forecast on the project's own definition |
| `--first-init` / `--last-init` defaults | 2025-11 / 2026-10 | User request (initially Feb 2026 to match the reference; moved to Nov 2025 the same day) |
| Trail default | off | User chose "pure flashlight"; the faded trail was added afterward as an option rather than replacing the default, at user request |
| `TRAIL_COLOR`, `TRAIL_BAND_ALPHA`, `TRAIL_LINE_ALPHA` | "0.5", 0.11, 0.5 | Cosmetic; tuned with the user (gray, then "a little darker" from 0.07/0.35) |
| `FRAME_MS`, `LAST_FRAME_HOLD` | 1000 ms, 3 | Cosmetic |
| Climatology / scaling period | 1991-2020 | `config.CLIM_START_YEAR`–`CLIM_END_YEAR`, project-wide |

## Edge Cases & Error Handling

- **Gap between path and beam.** The seasonal beam starts at lead 1 (lead 0
  is NaN after the centered running mean), and the obs path ends at the
  last complete season (needs the following month's obs). For the Oct 2026
  issue this leaves ASO and SON uncovered. Not filled — see Open Items.
- **Model set varies by frame.** GFDL-SPEAR enters the MMM only when it has
  that init (2026: Jan, Apr, Jul, Oct), per the project's use-when-available
  policy; `beam()` prints the model list per init.
- **Positional init range.** Inits are `range(i0, i1 + 1)` over S positions;
  valid because the store's S axis is a contiguous monthly sequence
  (all-NaN starts remain as coordinates rather than gaps).
- **Duplicate last frames.** Pillow collapses identical consecutive frames
  into one frame with summed duration, so the GIF reports n_inits frames,
  the last one 3000 ms.
- **Overlapping trail bands** darken where earlier forecasts agree; this is
  intended (shows consensus), not an artifact to suppress.

## Open Items

1. **Bridge the path-to-beam gap?** Decision: leave the gap (current, honest
   about what the seasonal forecast defines) vs. CPC-style blending of
   observed months with lead-0/lead-1 to fill the first seasons. Blending is
   a method change and would need its own spec entry.
2. **Spread reference vs. plotted path.** Band widths come from errors
   against observed *absolute* Niño-3.4 (`_lead_error_cov`'s convention),
   while the plotted path is the *scaled relative* index. Decision: accept
   (consistent with `latest_forecast`) vs. verify against the scaled
   relative obs for this figure. Matters if the figure goes public.

## Verification Snippet

```python
# Run from scripts/: mamba run -n pangeo-local python -c "<this>"
import numpy as np
import config
from latest_forecast import _resolve_init_idx
from forecast_flashlight import beam, obs_n34r_seasonal, BAND_QUANTILES

ds = config.load_nino34_verification()
start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
spec = {"var": "ssta_rel", "prefix": "n34r", "name": "relative Nino 3.4",
        "scale": config.rel_scaling_factor(ds)}
b = beam(ds, start, spec, _resolve_init_idx(start, "2026-10"))
assert b.sizes["L"] == 10, b.sizes                     # leads 1..10 survive the seasonal mean
qs = [f"q{round(q * 100):02d}" for pair in BAND_QUANTILES for q in pair]
for v in ["mean", *qs]:
    assert b[v].notnull().all(), v
assert ((b.q05 <= b.q15) & (b.q15 <= b.q25) & (b.q25 <= b["mean"])
        & (b["mean"] <= b.q75) & (b.q75 <= b.q85) & (b.q85 <= b.q95)).all()
assert np.allclose(b["mean"] - b.q25, b.q75 - b["mean"])  # symmetric Gaussian bands
obs = obs_n34r_seasonal()
assert obs.sel(time=slice("2025-10", "2026-08")).notnull().all()
print("ok")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-10-19 | Initial `scripts/forecast_flashlight.py`: fixed-axis animated n34r seasonal forecast GIF (Gaussian synthetic-spread bands, scaled ERSSTv5 relative-index obs path), modeled on a user-supplied NMME ENSO Outlook RONI GIF. Added `--trail` (faded gray earlier MMMs + 5–95% bands, separate `_trail` file); default first init changed Feb 2026 → Nov 2025 the same day. New `config.PLOTS_DIR_FORECAST_FLASHLIGHT`. | ✓ |
| 2026-10-19 | `beam()` now also returns `sd` (per-lead hindcast-error sd) for `forecast_flashlight_video.py`; comment documenting `int(L)` on NMME's 0.5-offset lead labels. No change to GIF output. | ✓ |
| 2026-11-13 | Default `--last-init` 2026-10 → 2026-11 (new NMME init). Regenerated both GIFs (now 13 frames); removed the stale `_2025-11_2026-10` outputs. No method change. | ✓ |
