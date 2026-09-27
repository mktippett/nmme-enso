# Computing the Relative Niño-3.4 Index in Forecast Models

This note documents how this project computes the **relative Niño-3.4 index**
for NMME forecast models — including the variance-scaling step, where the
published observational recipe does not directly carry over to models and a
choice must be made. The choice made here, and the evidence behind it, are
written down because we could not find them written down anywhere else.

## 1. Why a relative index

ENSO is conventionally monitored with the Niño-3.4 index: the SST anomaly
averaged over 5°S–5°N, 170°W–120°W. In a warming climate this becomes
ambiguous — the anomaly mixes ENSO variability with the background warming
trend, and event classification drifts as the climatology period is updated.

[Van Oldenborgh et al. (2021)](https://doi.org/10.1029/2021GL095041) proposed
subtracting the tropical-mean SST anomaly, and
[L'Heureux et al. (2024)](https://doi.org/10.1175/JCLI-D-23-0406.1)
developed the operational form used by NOAA CPC, the **relative Niño-3.4
index** (also called RONI):

$$\text{RONI} = \big[\,\text{Niño-3.4 anom} - \langle\text{tropical SST anom}\rangle_{20S-20N}\,\big] \times \frac{\sigma(\text{Niño-3.4 anom})}{\sigma(\text{Niño-3.4 anom} - \langle\text{tropical anom}\rangle)}$$

The subtraction removes the common warming signal; the scaling factor restores
the familiar amplitude of the standard index, so that thresholds like ±0.5 °C
keep their usual meaning. In observations both standard deviations are
computed from the observational record per calendar month, giving a single
scalar per month (≈ 1.1–1.2).

## 2. The question this project had to answer

The paper defines the scaling factor **in observations**. To plot and verify
model forecasts of the relative index on the same footing, the factor must be
defined **for a model forecast** — and the observational recipe does not
transfer verbatim, for two reasons:

1. **Models have their own variance biases.** If the model's relative-index
   variance differs from the observed (some NMME models' ENSO is far too
   energetic), scaling by an obs-only ratio leaves that bias in place.
2. **Forecast variance depends on start month and lead.** A model's
   climatology, drift, and ENSO amplitude all vary with initialization month
   and forecast lead, so a single per-calendar-month number is not stratified
   finely enough.

Per direct guidance from the paper's first author, the principle adopted here
is: *scale the model's relative Niño-3.4 variance to match the observed
1991–2020 Niño-3.4 variance.* Concretely
(`config.rel_scaling_factor`):

$$f(\text{model}, \text{start month}, \text{lead}) = \frac{\sigma_{\text{obs}}(\text{Niño-3.4 anom})}{\sigma_{\text{model}}(\text{relative Niño-3.4 anom})}$$

with both standard deviations stratified by start month and lead, over
1991–2020 forecast starts. This extends the paper in two deliberate ways:
the denominator is the **model's own** relative-index spread (making the
factor model-dependent and bias-correcting), and the stratification is by
**(start month, lead)** rather than calendar month alone. Computed factors
for the seven NMME models range from about 0.5 to 2.4 — i.e., some
model/month/lead combinations need their relative index *halved* to match
observed variance, which is exactly the bias an obs-only factor would ignore.

The model's forecast anomalies themselves (numerator of the relative index)
are computed against the model's own hindcast climatology — per model, per
start month, ensemble-mean seasonal cycle over 1991–2020 targets, with a
split-climatology treatment for three models whose hindcasts have a known
configuration discontinuity in 1999 (see `specs/latest_forecast.md`).

## 3. The subtle part: how to pool ensemble members

The denominator $\sigma_{\text{model}}$ must be estimated from an ensemble
hindcast with dimensions (start `S`, member `M`). There are three natural
estimators, and they are not equivalent:

| | Estimator | Formula (xarray) |
|---|---|---|
| **A** | per-member | `sqrt(x.groupby("S.month").var("S").mean("M"))` |
| **B** | ensemble-mean | `x.mean("M").groupby("S.month").std("S")` |
| **C** | grand-mean *(chosen)* | `sqrt(x.groupby("S.month").var(["S", "M"]))` |

**B is flawed by construction.** Averaging members before taking the variance
shrinks the variance wherever the ensemble disagrees — i.e., wherever skill
is low. A scaling factor built on B therefore mixes *predictability* into a
quantity that is meant to capture only *variance*, contaminating the very
comparison the scaling is supposed to normalize. The damage grows with lead,
as members decorrelate.

**A and C are both valid member-pooling estimators**, equal in expectation
but different in finite samples: C is the pooled standard deviation over the
flattened (start × member) sample and so additionally carries the
between-member spread of the per-member time means (law of total variance),
which A averages out.

### Evidence

`scripts/rel_scaling_compare.py` verifies all three pairwise comparisons
against a common observational target (the paper's own RONI), over 1991–2020:

- **Anomaly correlation is identical across A, B, C** (max difference
  ~10⁻¹⁵). This is expected — correlation is invariant to a positive scalar
  factor within each (model, start-month, lead) stratum — and serves as a
  sanity check that the three runs differ only in the scaling.
- **MSESS separates them.** Mean squared error skill *is* amplitude-sensitive:
  - MSESS(C) − MSESS(B) = **+0.067** averaged over model × month × lead —
    the grand-mean estimator beats the flawed ensemble-mean estimator
    decisively, and the gap grows with lead, consistent with the
    skill-contamination mechanism.
  - MSESS(C) − MSESS(A) = **+0.015** — the two valid estimators are nearly
    equivalent in practice, with grand-mean marginally ahead.

![MSESS, ensemble-mean minus grand-mean scaling](../plots/rel_scaling_compare/n34r_scaling_msess_diff_BC_start.png)

*MSESS(B) − MSESS(C) by start month and lead: predominantly negative,
increasingly so at long lead — the ensemble-mean denominator degrades the
scaled forecast exactly where members have decorrelated.*

**C (grand-mean pooling) is the production choice** — not worse than A on
the evidence, and the simpler estimator to state ("pool all members and
starts together") and to implement.

### Why it matters in practice

A concrete illustration from the July 2026 initialization (computed before
the §5 ocean-mask change, so the relative-index values differ slightly
from the current definition), during a strong
developing El Niño: peak ensemble-mean forecast anomalies across the seven
NMME models span **+2.6 to +4.9 °C** in the standard Niño-3.4 index — the
upper end driven by models with known excessive ENSO variance — but
**+2.3 to +3.7 °C** in the scaled relative index, with the across-model
spread of the peak reduced from 0.71 to 0.50 °C. The variance matching does
real work: it makes the models' amplitudes comparable to each other and to
observations before any forecast interpretation happens.

## 4. Recipe summary

Given an NMME-style hindcast/forecast dataset with dims
(model, start `S`, member `M`, lead `L`):

1. **Box averages**: cosine-latitude-weighted Niño-3.4 (5°S–5°N, 170°W–120°W)
   and tropical-mean (20°S–20°N, all longitudes) SST. The tropical mean is
   taken over a **common ocean mask**, the same cells for every model
   (§5).
2. **Model anomalies**: subtract the model's own ensemble-mean seasonal
   cycle, per model and start month, over 1991–2020 (respecting any known
   hindcast discontinuities). Apply the identical procedure to both box
   averages.
3. **Unscaled relative index**: `n34r_raw = n34_anom − trop_anom`.
4. **Scaling factor**, per (model, start month, lead), over 1991–2020 starts:

   ```python
   num = obs_n34_anom.groupby("S.month").std("S")             # obs, target-aligned
   den = np.sqrt(n34r_raw.groupby("S.month").var(["S", "M"])) # grand-mean pooling
   factor = num / den
   ```

5. **Scaled relative index**: `n34r = n34r_raw * factor`, selecting the
   factor by the forecast's model, start month, and lead.
6. If forecasts are smoothed (e.g., 3-month running mean over lead for
   seasonal display), compute a **separate factor from the smoothed fields**
   — running-meaning changes both variances.

Implementation: `config.rel_scaling_factor` (factor),
`config.load_nino34_verification` (anomalies and obs alignment),
`scripts/latest_forecast.py` (application), `scripts/rel_scaling_compare.py`
(evidence). Behavioral details and edge cases: `specs/latest_forecast.md`,
`specs/rel_scaling_compare.md`.

## 5. The tropical-mean ocean mask

*Decided 2026-09-27. Implementation: `config.tropics_ocean_mask`,
`config.TROPICS_MASK_EXCLUDE_GROUPS`.*

### Definition

The 20°S–20°N mean is taken over the cells that are **ocean in every NMME
model that supplies a land mask**. A cell counts as ocean in a model if the
model has data there at any start, member, or lead; the mask is the
intersection over those models (currently COLA-RSMAS-CCSM4,
COLA-RSMAS-CESM1, GFDL-SPEAR, NASA-GEOSS2S, NCEP-CFSv2), and the same
10,859 cells (1° grid) are used for all seven models. The observed relative
index uses ERSSTv5's own land mask, which is a separate grid and was not
harmonized with this one.

### Why a mask was needed

Before this change, the tropical mean relied on land being NaN. That holds
for five models (land fraction 0.23–0.26 of 20°S–20°N cells) but **not for
the two CanSIPS models** (CanSIPS-IC4-CanESM5, CanSIPS-IC4-GEM52NEMO), whose
land cells carry smooth SST-like values, e.g. 24–26 °C over the Amazon
(5°S, 300°E), with no NaN anywhere in 20°S–20°N. Their tropical means
therefore included about 3,400 land cells. This is a defect in the
source data (it is how the CanSIPS fields arrive from the archive), recorded
as an open item for NMME-zarr.

We also considered the NMME land–sea mask distributed with the IRI Data
Library (`SOURCES/.Models/.NMME/.LSMASK/.land`). It matches none of the
models' ocean masks exactly (e.g. 87–161 cells of GEOS/CFSv2 data lie on
its land, and 104–486 of its ocean cells are NaN in CFSv2/SPEAR), and its
provenance is undocumented, so it was not used.

### Effect

Relative to the previous (each model's own NaN pattern) calculation, over
1991–2020:

| | CCSM4 | CESM1 | CanESM5 | GEM5.2 | SPEAR | GEOS | CFSv2 |
|---|---|---|---|---|---|---|---|
| Scaling factor change, median % (max) | 0.6 (1.0) | 0.7 (1.0) | 3.0 (5.4) | 5.9 (9.2) | 0.1 (0.1) | 0.5 (0.9) | 0.6 (1.3) |
| n34r AC change, max over lead | .001 | .000 | .011 | .028 | .000 | .001 | .000 |

The change is material only for the two CanSIPS models; for the others it
is under about 1%.

### Known cold coastal cells, deliberately not masked

CCSM4 and CESM1 each have about 330 cells *inside* the mask (334 and 327)
whose long-term mean is more than 3 °C below the median of the other models
at the same cell; the coldest in-mask climatology is about 11.5 °C.
(Another ~450 such cells per model, some as cold as ~3 °C, are already
outside the mask because another model treats them as land.) About 95% of them are adjacent
to land in the common mask (e.g. NW Australia, Madagascar, Gulf of Tonkin).
The values are identical in the IRI source (CCSM4 at 19°S, 122°E reads about
10 °C next to 24 °C and 29 °C cells), so they are not an NMME-zarr fetch
error. Their behavior is consistent with coastal cells that blend in part
of an adjacent land value: the anomalies still track the other models
(median correlation 0.92 with the median of SPEAR, GEOS and CFSv2) but are
damped by about 0.7×, close to the ratio of the climatologies. The fit has a
nonzero intercept (median 2.2 °C), so this is not a simple zero-filled blend.

They are left in because they have **no model-specific effect**. Two
alternatives were evaluated against the mask above:

| Alternative | Cells | Scaling factor change vs. chosen mask, median % |
|---|---|---|
| Also drop any cell with a cross-model check failure (> 3 °C below the other models' median) | 10,499 | 0.5–0.7, all models |
| Also drop a one-cell coastal ring | 9,894 | 1.2–1.9, all models |

In both cases CCSM4 and CESM1 move no more than the unaffected models, so
either alternative would redefine the tropical mean for every model rather
than correct CCSM4/CESM1, at the cost of one more choice to justify. The
cold offset itself is removed by the climatology; what remains is a
damped anomaly on about 3% of the area.

**Why a 3 °C threshold:** at 3 °C the cross-model check isolates exactly this
coastal defect (0–9 flagged cells in any other model). At 1 °C it instead
picks up basin-wide model biases in the open ocean (about 800 cells for
CESM1, about 2,000 for GEM5.2), so the lower threshold is not a reliable
contamination test. A neighbor check within a single model is weaker still:
it misses contamination two cells deep.

## References

- L'Heureux, M. L., M. K. Tippett, et al., 2024: A Relative Sea Surface
  Temperature Index for Classifying ENSO Events in a Changing Climate.
  *J. Climate*, **37**, 1197–1211,
  [doi:10.1175/JCLI-D-23-0406.1](https://doi.org/10.1175/JCLI-D-23-0406.1).
- van Oldenborgh, G. J., et al., 2021: Defining El Niño indices in a warming
  climate. *Environ. Res. Lett.*, **16**, 044003,
  [doi:10.1088/1748-9326/abe9ed](https://doi.org/10.1088/1748-9326/abe9ed).
