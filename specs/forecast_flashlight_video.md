# forecast_flashlight_video.py — Behavioral Specification

> Last reviewed against code: 2026-11-13

## Purpose

A literal-flashlight "night walk" video companion to
`forecast_flashlight.py`: a walker carries a flashlight along the observed
relative Niño-3.4 (n34r) path, and each month's NMME forecast is the beam,
with brightness proportional to forecast probability density. The observed
path ahead is hidden in the dark and appears only where the light falls on
it, so forecast error shows up as the path leaving the beam — no separate
scoring. Exploratory/outreach prototype, not a manuscript figure.

## Inputs

| File | Relevant variables | Filters applied |
|------|--------------------|-----------------|
| `config.load_nino34_verification()` (via `forecast_flashlight.beam`) | `ssta_rel`, `obsa` | Inits `--first-init` … `--last-init` (default 2025-05 … 2026-11); available models per init |
| `config.ERSSTV5_NC` (via `forecast_flashlight.obs_n34r_seasonal`) | `n34_anom`, `rel_anom` | Path drawn within the fixed x-window |

All forecast and observed quantities come from `forecast_flashlight.py`'s
`beam()` (MMM `mean`, per-lead hindcast-error `sd`) and
`obs_n34r_seasonal()`; see `specs/forecast_flashlight.md` Algorithm §2–3.
This script adds no statistics, only rendering.

## Outputs

| File | Contents | Format |
|------|----------|--------|
| `plots/forecast_flashlight/n34r_seasonal_<first>_<last>_nightwalk.mp4` | Night-walk animation | MP4 (H.264, yuv420p, CRF 20), 1280×720, 30 fps; ≈ (n_inits−1)·1 s + 3 s (default 19 inits: 631 frames, 21.0 s) |
| `--stills T1,T2,… --stills-dir DIR` | Single frames at issue-time values T (0 = first issue; fractional = interpolated) instead of the video | PNG `still_<T>.png` |

## Algorithm

Time is a continuous issue coordinate τ ∈ [0, n−1]; τ = k is the k-th
real issue.

1. **Month index.** All x positions are month indices (year·12 + month − 1);
   seasons are placed at their center month.
2. **Beam knots per issue** (`beam_knots`): source at x = init − 2 (center of
   the last season observed at issue time), y = observed n34r there,
   sd = `SD_SOURCE`; then leads 1…10 at x = target center month, y = MMM,
   sd = hindcast-error sd. Every issue has the same knot structure (source
   + 10 leads, source 3 months before lead 1), so knots interpolate
   elementwise between issues.
3. **Beam image** (`render`): for each image column x in the beam's span,
   μ(x), σ(x) by linear interpolation of the knots;
   brightness = (1 − exp(−φ(y; μ, σ)/`DENSITY_SCALE`)) × end-fade, where φ
   is the normal pdf and end-fade ramps to 0 over the last `END_FADE_MONTHS`.
   `DENSITY_SCALE` is fixed for the whole video, so each column carries the
   same total probability and wider (more uncertain) columns are dimmer.
   The source-to-lead-1 neck is a visual connector, not a forecast.
4. **Interpolation between issues.** At τ = k + f, knots =
   (1−f)·knots_k + f·knots_{k+1} (linear, constant speed — no easing, no
   pauses); the image is re-rendered from the interpolated knots.
5. **Afterglow.** For each earlier issue j with age a = τ − j > 0:
   w_j = `GLOW_GAIN` · min(1, a/`GLOW_RISE`) · (`GLOW_FLOOR` + (1 −
   `GLOW_FLOOR`)·exp(−a/`GLOW_TAU`)); glow = pixelwise max_j(w_j · img_j).
   Light = clip(current + glow, 0, 1).
6. **Observed path.** Behind the walker (x ≤ source x): solid gray line with
   season dots, plus a partial segment to the walker's interpolated
   position. Ahead: the full observed path, densely sampled
   (`PATH_SAMPLES_PER_MONTH`), drawn in `LIT_PATH_RGB` with per-segment
   opacity clip(light(x, y)/`LIT_FULL`, 0, 1)^`LIT_GAMMA`, light sampled
   from the current frame's light image (nearest pixel).
7. **Flashlight.** A 4-part polygon glyph defined in screen pixels (lens at
   origin, pointing +x), scaled by `FLASHLIGHT_SCALE`, rotated by the
   screen-space angle from the source to the lead-1 knot (so the angle is
   correct despite the axes' aspect), translated to the source, and mapped
   back to data coordinates each frame.
8. **Issue label.** "Issued <Mon YYYY>" of the nearest real issue
   round(τ), alpha 1 − 0.9·|τ − round(τ)|.
9. **Fixed axes.** x from (first source − 3 months) to the last issue's
   last target; y from min/max of the in-window obs and all beams' 5%/95%
   bounds, floored/ceiled to 0.5 °C. Year labels only on Jan/Jul ticks.
10. **Frames.** τ = linspace(0, n−1, (n−1)·`SECONDS_PER_ISSUE`·`FPS` + 1),
    then `FINAL_HOLD_S`·`FPS` copies of τ = n−1.

## Constants & Scientific Rationale

All values below are display/animation choices made with the user while
prototyping on 2026-10-19; none affects the forecast statistics.

| Name | Value | Why |
|------|-------|-----|
| `--first-init` / `--last-init` | 2025-05 / 2026-10 | User request (start May 2025; covers the weak 2025 La Niña into the 2026 El Niño) |
| `SECONDS_PER_ISSUE`, `FINAL_HOLD_S`, `FPS` | 1.0, 3.0, 30 | User asked for smooth progress instead of pauses; final hold so the last forecast can be read |
| `DENSITY_SCALE` | 0.9 | Tone-map scale; fixed across frames so dimming reflects uncertainty |
| `SD_SOURCE` | 0.04 °C | Lens width; visual only |
| `END_FADE_MONTHS` | 0.6 | Softens the beam's far edge |
| `GLOW_GAIN`, `GLOW_TAU`, `GLOW_FLOOR`, `GLOW_RISE` | 0.5, 2 months, 0.15, 0.5 months | User asked for an afterglow that dims but never goes to zero; max-combination (not sum) chosen because a summed floor saturated by the end of the run |
| `LIT_FULL`, `LIT_GAMMA` | 0.8, 1.5 | Path ahead must stay dark under the dim afterglow floor and under the beam's faint tails; 0.45/1.0 lit the path almost everywhere |
| `FLASHLIGHT_SCALE` | 1.2 | User request (75% of 1.6) |
| `PATH_SAMPLES_PER_MONTH` | 12 | Smooth opacity gradient along the path |
| `NX`, `NY` | 900, 500 | Image resolution; bilinear-interpolated to 1280×720 |

## Edge Cases & Error Handling

- **NMME lead labels** are 0.5, 1.5, …; leads are located by target month
  (the `time` coord from `beam()`), never by `L` label.
- **Neck is not forecast.** For the Oct 2026 issue it spans ASO–SON and
  climbs from obs 1.8 to MMM ≈ 3.0 at OND; drawn as a narrow bright cone.
  See Open Items.
- **Path ahead at the end.** Observations end at the last walker position
  (init − 2 for the last issue), so the final frames have no lit path ahead.
- **Afterglow max vs. sum.** Pixelwise max loses the "consensus brightening"
  that the GIF's summed gray trail shows; accepted to keep the current beam
  legible.
- **Model count varies by issue** (SPEAR when available), inherited from
  `beam()`.
- **ffmpeg** comes from Homebrew (`/opt/homebrew/bin/ffmpeg`), not the conda
  env; matplotlib's `FFMpegWriter` finds it on PATH.

## Open Items

1. **Neck rendering.** Decision: keep the straight-line neck from the hand
   to lead 1 (current), dim/dash it to mark it as non-forecast, or bridge
   with CPC-style observed/lead-0 blending (method change). Same question as
   `specs/forecast_flashlight.md` Open Item 1.
2. **Hard beam ends.** Each beam ends at lead 10; the afterglow shows
   vertical edges. Decision: longer `END_FADE_MONTHS` vs. accept.
3. **Flashlight aim.** Aims at lead 1 only (keeps glyph and neck aligned).
   Alternative: mean of the first few leads (steadier, but misaligned with
   the neck).
4. **Prototype status.** Not yet reviewed for public release; inherits
   `specs/forecast_flashlight.md` Open Item 2 (spread verified against
   absolute, not relative, Niño-3.4).

## Verification Snippet

```python
# Run from scripts/: mamba run -n pangeo-local python -c "<this>"
import numpy as np
import config
from latest_forecast import _resolve_init_idx
from forecast_flashlight import beam, obs_n34r_seasonal
from forecast_flashlight_video import beam_knots, render, _mi

ds = config.load_nino34_verification()
start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
spec = {"var": "ssta_rel", "prefix": "n34r", "name": "relative Nino 3.4",
        "scale": config.rel_scaling_factor(ds)}
obs = obs_n34r_seasonal()
obs_mi = {_mi(t): float(v) for t, v in zip(obs.time.values, obs.values)}
ks = []
for d in ["2025-05", "2026-10"]:
    i = _resolve_init_idx(start, d)
    x, mu, sd = beam_knots(beam(ds, start, spec, i), start[i], obs_mi)
    assert len(x) == 11 and np.all(np.diff(x) > 0) and np.all(sd > 0)
    assert x[1] - x[0] == 3                      # source = init-2, lead 1 = init+1
    assert np.isfinite(mu).all()
    ks.append((x, mu, sd))
gx = np.linspace(ks[0][0][0] - 3, ks[1][0][-1], 400)
gy = np.linspace(-2, 4, 300)
img = render(*ks[1], gx, gy)
assert img.min() >= 0 and img.max() <= 1
assert (img[:, gx < ks[1][0][0]] == 0).all()     # nothing lit behind the source
print("ok")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| 2026-10-19 | Initial prototype `scripts/forecast_flashlight_video.py`. Iterated with user the same day: hit-dot scoring (lead-3 PIT classes + tally) replaced by the lit-path-ahead encoding; afterglow given a nonzero floor and switched from sum to max; pauses replaced by constant-speed linear interpolation; default first init 2025-11 → 2025-05; rotated flashlight glyph added (scale 1.6 → 1.2). | ✓ |
| 2026-11-13 | Default `--last-init` 2026-10 → 2026-11 (new NMME init). Regenerated the MP4 (631 frames, 21.0 s); removed the stale `_2025-05_2026-10` output. No rendering change. | ✓ |
