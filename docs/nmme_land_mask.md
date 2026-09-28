# Land Masking in the NMME SST Fields

This note records what we found about how the NMME monthly SST hindcasts and
forecasts treat land and coastal cells, and how this project handles it. It
is written for readers outside this project, including the model
maintainers. Findings date from 2026-09-27; data are the IRI Data Library
NMME `MONTHLY/.sst` fields (1° grid, 181 × 360), as mirrored locally by
`~/claude/NMME-zarr`. All values quoted below were checked to be identical
between the local mirror and the IRI Data Library.

## Summary

1. **The two CanSIPS models do not mask land.** CanSIPS-IC4-CanESM5 and
   CanSIPS-IC4-GEM52NEMO carry SST-like values over land rather than
   missing values. The other five models set land to missing.
2. **The IRI Data Library NMME land mask is not the models' mask.** It
   nearly matches the two COLA-RSMAS models and matches none of the others,
   and its provenance is undocumented. We do not use it.
3. **COLA-RSMAS-CCSM4 and COLA-RSMAS-CESM1 have cold coastal cells.** These
   are cells that are ocean in every masked model, with a long-term mean
   about 6 °C (up to 16 °C) colder than the other models. The pattern is
   the same in both models, at every lead, and in every decade.
4. **In this project the land problem only affects the tropical-mean SST**
   and hence the relative Niño-3.4 index (n34r). The Niño-3.4 box is all
   ocean in every model, so the standard index is unaffected.
5. **Strategy:** the tropical mean uses a common ocean mask, the cells that
   are ocean in all five land-masked models. The cold coastal cells are
   left in, because removing them has no model-specific effect.
6. **Impact:** the n34r scaling factor changes by 3–9% for the CanSIPS
   models and by 1% or less for the others. n34r anomaly correlation
   changes by at most 0.03.

## 1. CanSIPS models: land is not masked

| Model | Land fraction (missing values) in 20°S–20°N | Value at 5°S, 300°E (Amazon) |
|---|---|---|
| COLA-RSMAS-CCSM4 | 0.233 | missing |
| COLA-RSMAS-CESM1 | 0.233 | missing |
| GFDL-SPEAR | 0.263 | missing |
| NASA-GEOSS2S | 0.235 | missing |
| NCEP-CFSv2 | 0.227 | missing |
| **CanSIPS-IC4-CanESM5** | **0.000** | **25.7 °C** |
| **CanSIPS-IC4-GEM52NEMO** | **0.000** | **26.3 °C** |

(Start January 2010, member 1, lead 0.5.) The CanSIPS land values are
smooth and SST-like, and their 20°S–20°N minimum is 18–20 °C, with no
missing values anywhere in the band. A spatial average that
relies on missing values to exclude land therefore includes several
thousand land cells for these two models, but not for the others.

The masked models also do not share a coastline. Each masks a slightly
different set of cells (10,883–11,404 valid cells in 20°S–20°N), and
GFDL-SPEAR's mask varies by about 8 cells between samples.

## 2. The IRI Data Library NMME land mask

The IRI Data Library distributes an NMME land–sea mask at
`SOURCES/.Models/.NMME/.LSMASK/.land` (0/1, same 1° grid, no missing
values; fetchable as `.../LSMASK/.land/data.nc`). We compared it with each
model's own mask (20°S–20°N):

| Model | Model has data where LSMASK says land | Model missing where LSMASK says ocean |
|---|---|---|
| COLA-RSMAS-CCSM4 | 0 | 21 |
| COLA-RSMAS-CESM1 | 0 | 21 |
| GFDL-SPEAR | 14 | 486 |
| NASA-GEOSS2S | 87 | 144 |
| NCEP-CFSv2 | 161 | 104 |
| CanSIPS (both) | 3,413 | 0 |

The mask agrees with the two COLA-RSMAS models to within 21 cells, and was
probably what was applied to them. It is not the mask of any other model,
and we found no documentation of where it comes from. Using it would also
not remove the cold coastal cells in §3, which it treats as ocean.

## 3. Cold coastal cells in COLA-RSMAS-CCSM4 and COLA-RSMAS-CESM1

### What we see

Within 20°S–20°N, and inside the common ocean mask of §4, CCSM4 has 334
cells and CESM1 has 327 whose long-term mean is more than 3 °C below the
median of the other models at the same cell. **310 of these cells are
shared by the two models.** About 95% are next to a cell that at least one
model treats as land. The long-term mean is over all starts, members, and
leads.

Averaged over the affected cells, the offset compared with the median of
GFDL-SPEAR, NASA-GEOSS2S and NCEP-CFSv2 is about **−6 °C**:

| | By lead (0.5 → 8.5 months) | By start decade (1990s → 2020s) |
|---|---|---|
| CCSM4 | −5.8 → −6.2 °C | −5.9 → −6.3 °C |
| CESM1 | −6.0 → −6.2 °C | −6.0 → −6.2 °C |

The offset is essentially constant from the first lead, which points to a
fixed property of the grid or its processing (for example, how coastal
ocean cells were regridded to 1°) rather than to model drift.

The affected cells are spread along tropical coastlines worldwide, with the
most in the Maritime Continent, New Guinea and northern Australia
(105–150°E), and others along Africa, the Indian Ocean coasts, Central
America and the Caribbean, and South America. The worst cells in CCSM4
(long-term means; "others" is the median of the three reference models):

| Location | CCSM4 | Others | Difference |
|---|---|---|---|
| 10°S, 143°E (Torres Strait) | 11.5 °C | 27.5 °C | −15.8 °C |
| 12°N, 124°E (Philippines) | 13.6 °C | 28.5 °C | −14.9 °C |
| 13°N, 125°E (Philippines) | 15.1 °C | 28.8 °C | −13.4 °C |
| 5°N, 356°E (Gulf of Guinea) | 13.8 °C | 27.3 °C | −13.4 °C |
| 4°S, 114°E (Java Sea) | 16.4 °C | 29.3 °C | −12.7 °C |
| 14°S, 141°E (Gulf of Carpentaria) | 15.6 °C | 28.2 °C | −12.6 °C |

A single example row that reproduces directly from the IRI Data Library
(start 1 Jan 2010, member 1, lead 0.5, 19°S, 110–125°E):

```
.../.NMME/.COLA-RSMAS-CCSM4/.MONTHLY/.sst/S/(1 Jan 2010)VALUE/M/1/VALUE/L/0.5/VALUE/Y/-19/VALUE/X/(110)(125)RANGEEDGES/data.nc
```

| Longitude | 119°E | 120°E | 121°E | 122°E | 123°E |
|---|---|---|---|---|---|
| CCSM4 | 29.6 | 29.2 | 24.1 | **9.8** | missing |
| CESM1 | 29.4 | 28.9 | 23.7 | **9.5** | missing |

In addition to the cells above, each model has roughly 450 more such cold
cells, some as cold as about 3 °C. These cells fall outside the common
ocean mask because another model treats them as land. A single-start
snapshot (January 2010, lead 0.5, ensemble mean) over 40°S–40°N finds about
1,100–1,200 cells per model more than 3 °C colder than the reference
models, about 430 of them poleward of 20°. So the problem is not confined
to the tropics, but we have only examined the tropics in detail.

### What the anomalies do

Anomalies at the cold cells still track the other models. At lead 0.5, the
CCSM4 ensemble-mean series at the flagged cells correlates with the
reference-model median at 0.92 (median over cells; 10th percentile 0.73).
Its amplitude is damped: the regression slope is about 0.7 (median),
close to the ratio of the climatologies (0.81; median difference 0.08).
This is consistent with the coastal cells blending the ocean value with a
fraction of something colder, perhaps an adjacent land value. It is not a
simple blend with zero, since the fitted intercept is nonzero (median
2.2 °C). We have not determined the mechanism.

### How we detect them

A cross-model check on long-term means: flag a cell in model *m* if its
mean is more than 3 °C below the median of the other models' means at that
cell. At 3 °C this isolates the problem. CCSM4 and CESM1 have about 330
flagged cells each, and every other model has 0–9. At 1 °C the check
instead picks up large-scale model biases in the open ocean (about 800
cells for CESM1 and about 2,000 for GEM5.2-NEMO), so it is not a reliable
test for contamination. A spatial check within one model (a cell compared
with its own neighbors) is weaker: it misses contamination two cells deep,
where the neighbors are also cold.

## 4. What this project does

### Where it matters

The land problem matters here only for the **tropical-mean SST**
(20°S–20°N, all longitudes), which is subtracted from Niño-3.4 to form the
relative Niño-3.4 index (n34r; L'Heureux et al. 2024; see
[relative_nino34.md](relative_nino34.md)). The standard Niño-3.4 box
(5°S–5°N, 170°W–120°W) is entirely ocean in every model, so the standard
index and all its figures are unchanged (checked: byte-identical on
rerun).

### Strategy: a common ocean mask

The tropical mean is taken over the cells that are **ocean in every model
that masks land**, meaning all five models except CanSIPS. A cell counts
as ocean in a model if the model has data there at any start, member, or
lead. The mask is the intersection over the five models: 10,859 cells in
20°S–20°N, the same for all seven models (`config.tropics_ocean_mask`).

- This is part of the index definition, not only a repair. It puts every
  model on the same sample, and we would keep it even if CanSIPS land were
  masked at the source.
- The code checks its assumptions and stops rather than adapting silently:
  the CanSIPS land fraction must be exactly 0, and every other model's
  20°S–20°N land fraction must lie in 0.18–0.30.
- The observed relative index uses ERSSTv5's own land mask on its own grid.
  We did not harmonize the two.

### Cold coastal cells: left in, deliberately

We evaluated two ways of removing them, measured against the chosen mask
(change in the n34r variance-scaling factor, median over start month and
lead, 1991–2020):

| Alternative | Cells | Factor change |
|---|---|---|
| Also drop cells flagged by the 3 °C cross-model check | 10,499 | 0.5–0.7%, every model |
| Also drop a one-cell coastal ring | 9,894 | 1.2–1.9%, every model |

In both cases CCSM4 and CESM1 change no more than the models without the
problem, so removing the cells would redefine the tropical mean for every
model without correcting anything specific to COLA-RSMAS. The cold offset
itself is removed when the climatology is subtracted. What remains is a
damped anomaly on about 3% of the area.

## 5. Impact

Relative to the previous calculation (each model's own missing-value
pattern), over 1991–2020:

| | CCSM4 | CESM1 | CanESM5 | GEM5.2 | SPEAR | GEOS | CFSv2 |
|---|---|---|---|---|---|---|---|
| Scaling factor change, median % (max) | 0.6 (1.0) | 0.7 (1.0) | 3.0 (5.4) | 5.9 (9.2) | 0.1 (0.1) | 0.5 (0.9) | 0.6 (1.3) |
| Tropical-mean anomaly, RMS change (°C) | 0.007 | 0.009 | 0.028 | 0.059 | 0.001 | 0.007 | 0.008 |
| n34r anomaly correlation change, max over lead | 0.001 | 0.000 | 0.011 | 0.028 | 0.000 | 0.001 | 0.000 |

For reference, the tropical-mean anomaly standard deviation is 0.19–0.29 °C
across models. For the September 2026 forecast, the n34r multi-model mean
rose by 0.02–0.04 °C at every lead. The comparisons of scaling-factor
estimators in [relative_nino34.md](relative_nino34.md) §3 are unchanged to
three decimals.

A sibling project, `~/claude/enso-indices-us`, checked the cold cells
against the ENSO Longitude Index (ELI), which thresholds absolute SST. It
uses a common mask built the same way but requires data at two sampled
starts, giving 10,815 cells. About 48 cold cells per COLA-RSMAS model fall
in its 5°S–5°N, 115–282°E band. Excluding them shifts the COLA-RSMAS
models' climatological ELI about 1.5° west. The DJF ELI correlation
changes by 0.01 or less, and MSESS by 0.02 or less, except CCSM4's March
start, where it drops from 0.36 to 0.32. That project also left the cells
in.

## 6. Open items

- **CanSIPS land in the source store.** Whether `~/claude/NMME-zarr` should
  set CanSIPS land to missing at build time, and from which mask, and add
  a build check that fails when a model's 20°S–20°N land fraction is 0.
  Other projects that read the store are not protected by this project's
  mask.
- **Cold coastal cells: being reported to the COLA-RSMAS maintainers.** If they
  fix the cells upstream, the cross-model check above should find none.
  Rerun it and the §4 comparisons before relying on that.
- **Aligning the mask rule with enso-indices-us** ("data at any sample"
  here, "data at two sampled starts" there; 10,859 vs 10,815 cells). This
  is a consistency decision, not a results issue.
