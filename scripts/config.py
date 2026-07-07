"""
config.py — shared constants for NMME ENSO forecast skill / index diagnostics.

Import in each analysis script:
    import sys, config               # if scripts/ is on sys.path
    from config import STORE_SST, N34_LAT, nmme_groups

Store access
------------
The NMME zarr store lives in ~/claude/NMME-zarr/data/ and is maintained by
that project. Set NMME_STORE_DIR to override the default path, e.g.:
    export NMME_STORE_DIR=/data/nmme

Scope
-----
This project evaluates NMME skill at predicting the Niño-3.4 index itself
(anomaly correlation, RMSE, ensemble spread/reliability, event verification)
against ERSSTv5. Surface-temperature teleconnections are out of scope here —
see ~/claude/enso-t2m, which reads the same store for that purpose.

Store-derived helpers
---------------------
nmme_groups() and model_leads() read from the zarr store so the analysis
always reflects what is actually present — no separate model registry to
maintain.

Uniform sample discipline
--------------------------
All scripts use ANALYSIS_START_YEAR / CLIM_START_YEAR / CLIM_END_YEAR from
here rather than re-deriving a period filter locally.

Print-on-load pattern
----------------------
Store-derived helpers print what they found so silent failures are visible.
"""

import os
from pathlib import Path

ANALYSIS_DIR = Path(__file__).parent.parent  # nmme_enso/

# ---------------------------------------------------------------------------
# Store paths (env-var override with sensible default)
# ---------------------------------------------------------------------------
STORE_ROOT = Path(os.environ.get(
    "NMME_STORE_DIR",
    "/Users/tippett/claude/NMME-zarr/data",
))
STORE_SST = STORE_ROOT / "nmme_sst.zarr"

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
PLOTS_DIR = ANALYSIS_DIR / "plots"
OBS_DIR   = ANALYSIS_DIR / "observations"
CACHE_DIR = ANALYSIS_DIR / "cache"

# Each plotting script writes to its own plots/<script-name> subdirectory.
PLOTS_DIR_LATEST_FORECAST    = PLOTS_DIR / "latest_forecast"
PLOTS_DIR_SKILL              = PLOTS_DIR / "skill"
PLOTS_DIR_REL_SCALING_COMPARE = PLOTS_DIR / "rel_scaling_compare"

# ---------------------------------------------------------------------------
# Niño-3.4 region
# ---------------------------------------------------------------------------
N34_LAT = slice(-5, 5)
N34_LON = slice(190, 240)

# ---------------------------------------------------------------------------
# Tropical mean region — for the relative Niño-3.4 index (L'Heureux, Tippett
# et al. 2024, J. Climate). Entire tropics, all longitudes.
# ---------------------------------------------------------------------------
TROPICS_LAT = slice(-20, 20)

# ---------------------------------------------------------------------------
# Analysis / climatology period
# ---------------------------------------------------------------------------
ANALYSIS_START_YEAR = 1991
CLIM_START_YEAR = 1991
CLIM_END_YEAR   = 2020

# Float equivalent (months since 1960-01-01) for decode_times=False scripts
S_ANALYSIS_START = float((ANALYSIS_START_YEAR - 1960) * 12)  # 372.0

# Two-climatology models: split hindcast period to avoid future-data contamination
#   Period 1: before TWO_CLIM_SPLIT  → climatology from period-1 starts only
#   Period 2: TWO_CLIM_SPLIT onward  → climatology from 1999–CLIM_END_YEAR starts
TWO_CLIM_GROUPS = {"COLA-RSMAS-CCSM4", "COLA-RSMAS-CESM1", "NCEP-CFSv2"}
TWO_CLIM_SPLIT  = 1999   # period 2 starts January of this year

# Maximum lead count across all NMME models (L dimension size in the store).
# Models NASA-GEOSS2S (9) and NCEP-CFSv2 (10) are shorter and NaN-padded.
N_LEADS_PLOT = 12

# ERSSTv5 local file (override with ERSSTV5_NC env var) — Niño-3.4 verification
# reference. Default is the repo-local copy that the README Quickstart's
# curl -z command bootstraps and refreshes (PSL updates it monthly).
ERSSTV5_NC = Path(os.environ.get(
    "ERSSTV5_NC",
    str(OBS_DIR / "ERSSTv5.sst.mnmean.nc"),
))

# ---------------------------------------------------------------------------
# Short display names — used in figure titles and legends
# ---------------------------------------------------------------------------
SHORT_LABELS = {
    "COLA-RSMAS-CCSM4":      "CCSM4",
    "COLA-RSMAS-CESM1":      "CESM1",
    "CanSIPS-IC4-CanESM5":   "CanESM5",
    "CanSIPS-IC4-GEM52NEMO": "GEM5.2-NEMO",
    "GFDL-SPEAR":            "SPEAR",
    "NASA-GEOSS2S":          "GEOSS2S",
    "NCEP-CFSv2":            "CFSv2",
}


# ---------------------------------------------------------------------------
# Store-derived helpers — read from the zarr store at call time
# ---------------------------------------------------------------------------

def nmme_groups(store=None):
    """Sorted list of model group names present in the zarr store.

    Parameters
    ----------
    store : path-like or None
        Path to the zarr store; defaults to STORE_SST.

    Returns
    -------
    list[str]
    """
    store = Path(store or STORE_SST)
    import zarr as _zarr
    root = _zarr.open_group(str(store), mode="r")
    groups = sorted(root.group_keys())
    print(f"  nmme_groups ({store.name}): {groups}")
    return groups


def model_leads(store=None):
    """Dict mapping group name → L dimension size.

    Parameters
    ----------
    store : path-like or None
        Path to the zarr store; defaults to STORE_SST.

    Returns
    -------
    dict[str, int]
    """
    store = Path(store or STORE_SST)
    import zarr as _zarr
    root = _zarr.open_group(str(store), mode="r")
    leads = {}
    for g in root.group_keys():
        grp = root[g]
        if "sst" in grp:
            # (S, M, L, Y, X) — L is index 2
            leads[g] = grp["sst"].shape[2]
    print(f"  model_leads ({store.name}): {leads}")
    return leads


def short_label(group):
    """Short display name for an NMME model group."""
    return SHORT_LABELS.get(group, group.rsplit("-", 1)[-1])


def _store_fingerprint(store):
    """Cheap per-group metadata fingerprint: S size, last S value, and the
    'last_updated' attribute update_archive.py stamps on each group.

    Reads only zarr array metadata (shape) and the last S value — never the
    sst data itself — so this is fast even against a many-GB store. Used to
    invalidate the load_nino34_ssta() disk cache: the store is append-only
    along S (update_archive.py appends new starts and may overwrite the
    trailing RECHECK_TAIL starts in place), so any real change is visible in
    this fingerprint without rereading the historical grid data.
    """
    import zarr as _zarr
    root = _zarr.open_group(str(store), mode="r")
    fp = {}
    for g in sorted(root.group_keys()):
        grp = root[g]
        S = grp["S"]
        fp[g] = {
            "s_size": int(S.shape[0]),
            "s_last": float(S[-1]) if S.shape[0] else None,
            "last_updated": grp.attrs.get("last_updated"),
        }
    return fp


# ---------------------------------------------------------------------------
# Niño-3.4 index + forecast anomaly loader
# ---------------------------------------------------------------------------

def _n34_average(x):
    """Cosine-latitude-weighted Niño-3.4 box average over X, Y."""
    import numpy as np
    x = x.sortby("Y")
    weights = np.cos(np.deg2rad(x.Y))
    y = x.sel(X=N34_LON, Y=N34_LAT).weighted(weights).mean(["X", "Y"])
    y.attrs = x.attrs.copy()
    return y


def _tropics_average(x):
    """Cosine-latitude-weighted tropical-mean average over X, Y (20S-20N, all lon).

    Land points are NaN in the underlying grid; weighted().mean() skips NaNs
    by default for floating dtypes.
    """
    import numpy as np
    x = x.sortby("Y")
    weights = np.cos(np.deg2rad(x.Y))
    y = x.sel(Y=TROPICS_LAT).weighted(weights).mean(["X", "Y"])
    y.attrs = x.attrs.copy()
    return y


def _decode_cf(ds, time_var):
    """Decode a 360-day-calendar time variable to cftime.

    xarray writes calendar='360' for this store, but cftime/CF decoding
    requires the spelling '360_day'; patch before decoding. Decodes every
    CF-time variable in ds, not just time_var (e.g. also decodes 'target'
    if it carries matching units/calendar attrs).
    """
    import xarray as xr
    if ds[time_var].attrs.get("calendar") == "360":
        ds[time_var].attrs["calendar"] = "360_day"
    return xr.decode_cf(ds, decode_times=True)


def _forecast_anomaly(sst_da, target):
    """Two-climatology-scheme ensemble-mean anomaly for a box-averaged sst.

    sst_da has dims (model, S, M, L); target is the cftime valid-time
    DataArray (dims S, L) used to build the fixed-climatology mask. This is
    the anomaly logic shared by the Niño-3.4 box average and the tropical-
    mean box average (load_nino34_ssta computes both and differences them
    to form the unscaled relative index) — do not collapse the two
    climatology schemes into one, see load_nino34_ssta docstring.

    Returns
    -------
    xarray.DataArray
        Dims (model, S, M, L), same shape as sst_da.
    """
    import cftime
    import xarray as xr

    clim_start = cftime.Datetime360Day(CLIM_START_YEAR, 1, 16)
    clim_end = cftime.Datetime360Day(CLIM_END_YEAR, 12, 16)
    climo_1991_2020 = (target >= clim_start) & (target <= clim_end)

    two_clim_models = sorted(TWO_CLIM_GROUPS & set(sst_da.model.values))
    one_clim_models = [m for m in sst_da.model.values if m not in TWO_CLIM_GROUPS]

    split_start = f"{TWO_CLIM_SPLIT}-01-01"
    split_prev_end = f"{TWO_CLIM_SPLIT - 1}-12-01"
    clim_end_str = f"{CLIM_END_YEAR}-12-01"

    f1 = sst_da.sel(model=two_clim_models).sel(S=slice(None, split_prev_end))
    ssta1 = f1.groupby("S.month") - f1.mean("M").groupby("S.month").mean("S")

    f2 = sst_da.sel(model=two_clim_models).sel(S=slice(split_start, None))
    ssta2 = f2.groupby("S.month") - (
        f2.sel(S=slice(split_start, clim_end_str))
        .mean("M").groupby("S.month").mean("S")
    )

    ssta_two_clim = xr.concat([ssta1, ssta2], dim="S")

    f = sst_da.sel(model=one_clim_models)
    fc = f.where(climo_1991_2020)
    ssta_one_clim = f.groupby("S.month") - fc.mean("M").groupby("S.month").mean("S")

    return xr.concat([ssta_two_clim, ssta_one_clim], dim="model")


def load_nino34_ssta(store=None, use_cache=True):
    """Load the NMME Niño-3.4 index and compute forecast anomalies (ssta).

    Reads every model group in the SST zarr store, reduces to the cosine-
    weighted Niño-3.4 box average and the cosine-weighted tropical-mean box
    average (20S-20N, all lon), computes cftime target (valid) time, and
    removes a per-model, per-init-month ensemble-mean climatology from each
    to form ssta and ssta_trop. ssta_rel = ssta - ssta_trop is the unscaled
    relative Niño-3.4 index (L'Heureux, Tippett et al. 2024, J. Climate) —
    scaling to match the variance of ssta (against the model's own
    relative-index variance) is applied by the caller (see
    rel_scaling_factor), not computed here.

    Two climatology schemes are used, matching a known NMME hindcast
    discontinuity in TWO_CLIM_GROUPS (applied identically to the Niño-3.4
    and tropical-mean box averages via _forecast_anomaly):
      - TWO_CLIM_GROUPS: climatology split at TWO_CLIM_SPLIT on *init* time
        (S). Starts before the split use a climatology from that segment's
        own starts; starts from the split onward use a climatology computed
        over [TWO_CLIM_SPLIT, CLIM_END_YEAR] starts only (a fixed baseline
        applied to all later starts, including forecasts past CLIM_END_YEAR).
      - All other models: a single climatology defined on *target* (valid)
        time, CLIM_START_YEAR-CLIM_END_YEAR.

    In both cases the anomaly removed is the ensemble-mean (mean over M)
    seasonal cycle, stratified by init month (S.month). Do not collapse
    these two schemes into one — the split exists because those models'
    hindcasts have a real configuration change at TWO_CLIM_SPLIT.

    Caching
    -------
    The full-history box-average read is the expensive step (the store is
    tens of GB), but the historical portion is append-only and unchanged
    between successive runs — only new forecast starts (and the trailing
    RECHECK_TAIL=2 starts update_archive.py rechecks) ever change. When
    use_cache=True (default), the returned ds is cached to
    CACHE_DIR/nino34_ssta.nc, keyed on _store_fingerprint(store); a
    matching fingerprint on the next call skips recomputation entirely and
    loads the small cached dataset instead. Pass use_cache=False to force
    a fresh read (e.g. after suspecting the fingerprint missed a change).

    Parameters
    ----------
    store : path-like or None
        Path to the SST zarr store; defaults to STORE_SST.
    use_cache : bool
        Read/write the disk cache in CACHE_DIR; True by default.

    Returns
    -------
    xarray.Dataset
        Dims (model, S, M, L). Variables: sst (Niño-3.4 index), trop
        (tropical-mean index), target (cftime valid time = S + L), ssta
        (Niño-3.4 forecast anomaly), ssta_trop (tropical-mean forecast
        anomaly), ssta_rel (unscaled relative Niño-3.4 anomaly).
    """
    import json

    import numpy as np
    import xarray as xr

    store = Path(store or STORE_SST)
    fingerprint = _store_fingerprint(store)

    cache_nc = CACHE_DIR / "nino34_ssta.nc"
    cache_fp = CACHE_DIR / "nino34_ssta.fingerprint.json"
    if use_cache and cache_nc.exists() and cache_fp.exists():
        if json.loads(cache_fp.read_text()) == fingerprint:
            time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
            ds = xr.load_dataset(cache_nc, decode_times=time_coder)
            print(
                f"  load_nino34_ssta: cache hit ({cache_nc.name}) "
                f"models={list(ds.model.values)} shape={dict(ds.ssta.sizes)}"
            )
            return ds
        print(f"  load_nino34_ssta: cache stale ({cache_nc.name}), recomputing")

    ds_list = []
    for group in nmme_groups(store):
        g = xr.open_zarr(str(store), group=group, decode_times=False)
        sst_full = g.sst
        n34 = _n34_average(sst_full).astype(np.float64).compute()
        trop = _tropics_average(sst_full).astype(np.float64).compute()
        g["sst"] = n34
        g["trop"] = trop
        g = g.drop_vars(["X", "Y", "_filled", "T"])
        g.coords["model"] = group
        g = g.expand_dims("model")
        ds_list.append(g)
    ds = xr.merge(ds_list)
    ds.S.attrs["units"] = "months since 1960-01-01"
    ds.S.attrs["calendar"] = "360"

    # target (valid) time = init + lead — computed while S is still a raw
    # float, before decoding to cftime (order matters, see load_nino34_ssta
    # docstring / spec).
    ds["target"] = ds.L + ds.S
    ds.target.attrs["units"] = ds.S.attrs["units"]
    ds.target.attrs["calendar"] = "360_day"

    ds = _decode_cf(ds, "S")  # decodes both S and target to cftime

    ds["ssta"] = _forecast_anomaly(ds.sst, ds.target)
    ds.ssta.attrs["long_name"] = "ssta"

    ds["ssta_trop"] = _forecast_anomaly(ds.trop, ds.target)
    ds.ssta_trop.attrs["long_name"] = "ssta_trop"

    ds["ssta_rel"] = ds.ssta - ds.ssta_trop
    ds.ssta_rel.attrs["long_name"] = "ssta_rel (unscaled relative Nino-3.4)"

    two_clim_models = sorted(TWO_CLIM_GROUPS & set(ds.model.values))
    print(
        f"  load_nino34_ssta: models={list(ds.model.values)} "
        f"shape={dict(ds.ssta.sizes)} ssta_rel_shape={dict(ds.ssta_rel.sizes)} "
        f"two_clim={two_clim_models} split={TWO_CLIM_SPLIT} "
        f"climo={CLIM_START_YEAR}-{CLIM_END_YEAR}"
    )

    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        # S/target inherit encoding={'units': 'months since ...'} from
        # _decode_cf (the store's raw calendar); "months" is a valid cftime
        # decoding unit but not a fixed-duration CF unit, so the netCDF
        # encoder rejects it on write. Drop it and let xarray pick a
        # supported default (calendar is inferred from the cftime dtype).
        ds_cache = ds.copy()
        ds_cache.S.encoding = {}
        ds_cache.target.encoding = {}
        ds_cache.to_netcdf(cache_nc)
        cache_fp.write_text(json.dumps(fingerprint))
        print(f"  load_nino34_ssta: wrote cache ({cache_nc.name})")

    return ds


def _obs_index_anomalies(nc=None, clim_period=(CLIM_START_YEAR, CLIM_END_YEAR)):
    """Niño-3.4 and relative-index ERSSTv5 monthly anomalies, full record.

    Returns cosine-weighted Niño-3.4 (n34_anom) and relative-index
    (rel_anom = n34_anom - tropical-mean anom) anomalies over the *entire*
    ERSSTv5 record, computed against a monthly climatology from
    `clim_period` only (inclusive years) — i.e. the climatology window is
    narrower than the data returned. Used by load_nino34_verification
    (which needs full-record coverage to align against forecast targets
    that fall outside `clim_period`) and, via its obsa output, by
    rel_scaling_factor's numerator.

    Parameters
    ----------
    nc : path-like or None
        Path to the ERSSTv5 monthly-mean netCDF; defaults to ERSSTV5_NC.
    clim_period : tuple[int, int]
        (start_year, end_year), inclusive, defining the climatology window.

    Returns
    -------
    xarray.Dataset
        Dims (time,), Gregorian monthly, full ERSSTv5 record. Variables:
        n34_anom, rel_anom.
    """
    import numpy as np
    import xarray as xr

    nc = Path(nc or ERSSTV5_NC)
    start_year, end_year = clim_period
    ds = xr.open_dataset(nc)
    sst = ds.sst

    def _box_mean(x, lat, lon=None):
        x = x.sortby("lat")
        weights = np.cos(np.deg2rad(x.lat))
        sel = x.sel(lat=lat) if lon is None else x.sel(lat=lat, lon=lon)
        return sel.weighted(weights).mean(["lat", "lon"])

    n34_obs = _box_mean(sst, N34_LAT, N34_LON)
    trop_obs = _box_mean(sst, TROPICS_LAT)

    clim_sel = dict(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
    n34_clim = n34_obs.sel(**clim_sel).groupby("time.month").mean("time")
    trop_clim = trop_obs.sel(**clim_sel).groupby("time.month").mean("time")

    n34_anom = n34_obs.groupby("time.month") - n34_clim
    trop_anom = trop_obs.groupby("time.month") - trop_clim
    rel_anom = n34_anom - trop_anom

    return xr.Dataset({"n34_anom": n34_anom, "rel_anom": rel_anom})


def rel_scaling_factor(ds, period=(CLIM_START_YEAR, CLIM_END_YEAR)):
    """(model, month, L) factor scaling the *model* relative Niño-3.4 index
    variance to match the *observed* Niño-3.4 variance.

    Per the first author of L'Heureux, Tippett et al. (2024, J. Climate):
    "scale the model relative Niño-3.4 variance to match the observed
    1991-2020 Niño-3.4 variance." This differs from the paper's own
    obs-only ratio (std_obs(n34)/std_obs(n34-trop), a single per-calendar-
    month number, formerly config.load_obs_scaling) in two ways that are
    new extensions beyond both the paper and the prior implementation:
      - the denominator is the *model's* relative-index std, not obs's —
        so the factor is model-dependent.
      - the factor is stratified by (start month, lead L) rather than
        target month alone, since NMME climatologies/skill vary by lead.

    factor(model, month, L) = std_obs(n34) / std_model(n34r), both computed
    over forecast starts S in `period` (default CLIM_START_YEAR-
    CLIM_END_YEAR) via groupby('S.month').std('S') on the (S, L) grid from
    load_nino34_verification(). The numerator (ds.obsa) has no model
    dimension and broadcasts across it.

    Ensemble-member pooling (denominator)
    --------------------------------------
    The denominator pools members via the grand mean: sqrt(var(['S', 'M']))
    — the pooled standard deviation over the flattened (start, member)
    sample — rather than std('S') of the ensemble mean (ds.ssta_rel.mean(
    'M')). This is deliberate: the ensemble mean shrinks variance wherever
    skill is low, which would mix predictability into a quantity meant to
    capture only variance (a scaling factor derived from the ensemble-mean
    std would itself depend on skill, contaminating the comparison it's
    meant to normalize). See scripts/rel_scaling_compare.py, which confirms
    this in practice: AC is identical to the ensemble-mean alternative (as
    expected, correlation is scale-invariant), and MSESS is higher for the
    grand-mean-pooled factor (mean +0.05 across model/month/lead, 1991-2020)
    than for the ensemble-mean alternative.

    There are TWO member-based estimators, both of which pool members
    without the ensemble-mean's skill contamination — they are equal in
    expectation but differ in finite samples:
      - grand-mean (used here): sqrt(var(['S', 'M'])) — the pooled standard
        deviation over the flattened (start, member) sample.
      - per-member: sqrt(var('S').mean('M')) — the temporal variance is
        computed separately per member, then averaged across members.
    The grand-mean form additionally carries the between-member spread of
    the per-member time-means (law of total variance), which the per-member
    form averages out. **Changed 2026-07-07: switched from per-member to
    grand-mean.** Per-member was the original choice (see Synchronization
    Log in specs/rel_scaling_compare.md for its rationale at the time); the
    switch to grand-mean followed a review of
    scripts/rel_scaling_compare.py's AC-vs-per-member evidence
    (msess_diff_AC), which showed MSESS is not worse for grand-mean
    (slightly higher on average, -0.015 mean MSESS(per-member) -
    MSESS(grand-mean) across model/month/lead, 1991-2020) — combined with
    grand-mean being the simpler estimator to describe ("pool all members
    together", vs. "average each member's own variance across members").
    scripts/rel_scaling_compare.py now compares all three pairings
    (grand-mean vs. ensemble-mean, grand-mean vs. per-member, ensemble-mean
    vs. per-member) for ongoing evidence.

    A separate factor is computed for the seasonal (3-month centered
    running mean over L) variant, since running-meaning changes both
    variances.

    Parameters
    ----------
    ds : xarray.Dataset
        Output of load_nino34_verification() — must carry ssta_rel (model,
        S, M, L) and obsa (S, L).
    period : tuple[int, int]
        (start_year, end_year), inclusive, restricting forecast starts S
        for both the numerator and denominator std.

    Returns
    -------
    tuple[xarray.DataArray, xarray.DataArray]
        (factor_monthly, factor_seasonal), each dims (model, month, L).
    """
    import numpy as np

    start_year, end_year = period
    year = ds.S.dt.year
    ref = ds.where((year >= start_year) & (year <= end_year), drop=True)

    num = ref.obsa.groupby("S.month").std("S")  # (month, L) — obs, model-independent
    den = np.sqrt(ref.ssta_rel.groupby("S.month").var(["S", "M"]))  # (model, month, L)
    factor_monthly = num / den

    obsa_s = ref.obsa.rolling(L=3, center=True).mean()
    rel_s = ref.ssta_rel.rolling(L=3, center=True).mean()
    num_s = obsa_s.groupby("S.month").std("S")
    den_s = np.sqrt(rel_s.groupby("S.month").var(["S", "M"]))
    factor_seasonal = num_s / den_s

    print(
        f"  rel_scaling_factor ({start_year}-{end_year}): "
        f"shape={dict(factor_monthly.sizes)} "
        f"monthly range=[{float(factor_monthly.min()):.2f}, {float(factor_monthly.max()):.2f}] "
        f"seasonal range=[{float(factor_seasonal.min()):.2f}, {float(factor_seasonal.max()):.2f}]"
    )

    return factor_monthly, factor_seasonal


def load_nino34_verification(store=None, use_cache=True):
    """load_nino34_ssta() plus ERSSTv5 observational anomalies aligned onto
    the same (S, L) grid, for forecast-skill diagnostics (correlation,
    MSESS, ...) against observations.

    Alignment problem
    ------------------
    ds.target (from load_nino34_ssta) is a 360-day cftime valid time
    stamped on the 16th of each month (target = S + L; L=0.5 is the init
    month itself). ERSSTv5's `time` is Gregorian, stamped on the 1st, with
    no calendar attribute. A float/timestamp join between the two — e.g.
    building a matching 360-day mid-month obs time axis and doing an
    exact-float `.sel(T=S+L)` — is fragile: it depends on an assumed
    day-of-month convention on both sides, and floating-point `.sel`
    equality can silently miss or (worse) silently match the wrong month
    if the two conventions ever drift.

    Instead, both grids are collapsed to an integer (year, month) index —
    the only thing that is actually invariant across calendar and
    day-of-month convention — and the gather is an exact-integer `.sel`:
        ym = (year - 1960) * 12 + (month - 1)
    `.dt.year` / `.dt.month` work identically on cftime Datetime360Day and
    on Gregorian datetime64 (see xarray-gotchas), so this maps the 16th
    and the 1st of the same (year, month) to the same integer with no
    day-of-month assumption at all.

    Observational anomalies are computed over ERSSTv5's *full* record (see
    _obs_index_anomalies), not restricted to CLIM_START_YEAR-CLIM_END_YEAR,
    so that forecast targets falling outside that window (e.g. a
    CLIM_END_YEAR-ish start whose 11-month-lead target lands the following
    year) still have a defined obs value; only the *climatology* removed
    to form the anomaly is fixed to CLIM_START_YEAR-CLIM_END_YEAR. Targets
    beyond the end of the ERSSTv5 record itself (genuine future forecasts)
    get NaN via reindex — ac_by_start/msess_by_start in skill.py already
    drop these via their own ~isnan masking.

    Parameters
    ----------
    store : path-like or None
        Path to the SST zarr store; defaults to STORE_SST.
    use_cache : bool
        Passed through to load_nino34_ssta (this function adds no caching
        of its own — the obs alignment is cheap relative to the store read).

    Returns
    -------
    xarray.Dataset
        load_nino34_ssta()'s output, plus obsa and obsa_rel: ERSSTv5
        Niño-3.4 / relative-index anomalies on the (S, L) grid (no model
        dim — broadcasts against ssta's/ssta_rel's model dim).
    """
    import numpy as np

    ds = load_nino34_ssta(store, use_cache)
    obs = _obs_index_anomalies(clim_period=(CLIM_START_YEAR, CLIM_END_YEAR))

    ym_obs = (obs.time.dt.year - 1960) * 12 + (obs.time.dt.month - 1)
    obs = obs.assign_coords(ym=ym_obs).swap_dims({"time": "ym"})

    ym_tgt = (ds.target.dt.year - 1960) * 12 + (ds.target.dt.month - 1)
    ym_full = np.arange(int(ym_tgt.min()), int(ym_tgt.max()) + 1)
    obs = obs.reindex(ym=ym_full)  # NaN for targets past the obs record

    ds["obsa"] = obs.n34_anom.sel(ym=ym_tgt)
    ds.obsa.attrs["long_name"] = "obsa (ERSSTv5 Nino-3.4 anomaly, target-aligned)"
    ds["obsa_rel"] = obs.rel_anom.sel(ym=ym_tgt)
    ds.obsa_rel.attrs["long_name"] = "obsa_rel (ERSSTv5 relative-index anomaly, target-aligned)"

    print(
        f"  load_nino34_verification: obsa shape={dict(ds.obsa.sizes)} "
        f"nan_frac={float(ds.obsa.isnull().mean().values):.3f} "
        f"climo={CLIM_START_YEAR}-{CLIM_END_YEAR}"
    )

    return ds


if __name__ == "__main__":
    print("STORE_SST:", STORE_SST)
    if STORE_SST.exists():
        print("Groups:", nmme_groups())
        print("Leads:", model_leads())
    else:
        print("Store not found — set NMME_STORE_DIR")
    print("config.py OK")
