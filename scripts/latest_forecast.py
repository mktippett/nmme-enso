"""
latest_forecast.py — NMME Niño-3.4 latest-forecast plume figures.

Reads the local NMME SST zarr store (config.STORE_SST) via
config.load_nino34_verification(), and plots the most recent forecast
initialization (and the previous one, for comparison) as a set of
monitoring figures. These are exploratory/monitoring figures, not
manuscript figures: they carry the climatology note in the title rather
than travelling with a separate caption.

Two indices are plotted, each into its own subfolder (plots/latest_forecast/n34/, n34r/):
  - Niño-3.4 (ds.ssta): the standard box anomaly.
  - relative Niño-3.4 (ds.ssta_rel, scaled): Niño-3.4 anomaly minus the
    tropical-mean anomaly, rescaled by a (model, start-month, lead)-
    dependent factor (config.rel_scaling_factor) so the model's own
    relative-index variance matches the observed (ERSSTv5) Niño-3.4
    variance over 1991-2020 — see L'Heureux, Tippett et al. (2024,
    J. Climate, papers/) for the index's origin; the specific
    model-relative, (month, lead)-stratified factor here is an extension
    beyond the paper's obs-only ratio (see specs/latest_forecast.md).

Each plume figure (compare, spread, mean) is produced in two variants:
  - monthly: the raw monthly forecast anomaly.
  - seasonal: a centered 3-month running mean over lead, labeled by target
    season (DJF, JFM, ...) rather than by month. The smoothing is a display
    choice made here, not in config.load_nino34_verification() — the loader
    always returns monthly ssta/ssta_rel so other scripts get unsmoothed
    data. The
    relative index's scaling factor is also computed separately for the
    monthly and seasonal variants (running-meaning changes the variance).

A ninth figure, strength_probabilities (n34r/seasonal only, no monthly
variant), is a stacked-bar chart of ENSO strength-category probability by
target season, format after NOAA CPC's ENSO Strength Probabilities chart
(https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/strengths/).
It reuses _synthetic_plume (the same MMM + historical-error-covariance
draws as plot_spread_synthetic) at a higher draw count for smoother category
percentages — see plot_strength_probabilities.

By default, plots the latest available initialization (and the one before
it, for Compare). Pass --init-date to plot a specific past initialization
instead — matched by calendar month (day is ignored, since NMME inits land
on the 16th of each month after the 360-day-to-Gregorian conversion; see
specs/latest_forecast.md). Output filenames get a "_<YYYY-MM>" suffix in
that case, so they don't overwrite the latest-init monitoring figures.

Run:
    mamba run -n pangeo-local python scripts/latest_forecast.py
    mamba run -n pangeo-local python scripts/latest_forecast.py --init-date 2026-06-01
"""
import argparse

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.patches
import matplotlib.pyplot as plt

import config

MMM_COLOR = "0.75"

# Synthetic-plume family (plot_spread_synthetic): equally-likely scenarios drawn
# from a Gaussian with the historical MMM forecast-error covariance across leads,
# added to the current MMM (Barnston, Tippett, van den Dool & Unger 2015, JAMC,
# https://doi.org/10.1175/JAMC-D-14-0188.1, Fig. 9 lower panels).
N_SYNTHETIC_MEMBERS = 100
SYNTHETIC_SEED = 0
SYNTHETIC_MEMBER_COLOR = "#7a86c8"
SYNTHETIC_MMM_COLOR = "#1a1a1a"

# plot_strength_probabilities draws far more synthetic members than the
# plotting default above — smooth stacked-bar category percentages need a
# finer empirical distribution than a plume figure does (cf.
# kalshi_roni_pricing.py's analogous override for its bucket probabilities).
N_STRENGTH_DRAWS = 5000

# Single-letter month initials, index 0 = January (the notebook's "m_str").
SEASON_INITIALS = "JFMAMJJASOND"


def _season_label(month):
    """3-letter overlapping-season label centered on `month` (1-12), e.g. Jan -> DJF."""
    i = month - 1
    prev_m = SEASON_INITIALS[(i - 1) % 12]
    this_m = SEASON_INITIALS[i]
    next_m = SEASON_INITIALS[(i + 1) % 12]
    return prev_m + this_m + next_m


def _seasonal(da):
    """Centered 3-month running mean over lead. Endpoints (first/last lead) are NaN."""
    return da.rolling(L=3, center=True).mean()


def _plume_colors():
    """Default color cycle, with the first three overridden for plume plots."""
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    plume_colors = colors.copy()
    plume_colors[0] = "#005030"
    plume_colors[1] = "#F17221"
    plume_colors[2] = "#1f77b4"
    return plume_colors


def _resolve_init_idx(start, init_date):
    """Positional S index for `init_date` (matched by calendar year/month,
    since inits land on the 16th of each month), or -1 (latest) if None."""
    if init_date is None:
        return -1
    target = pd.Timestamp(init_date)
    matches = np.where((start.year == target.year) & (start.month == target.month))[0]
    if len(matches) == 0:
        raise ValueError(
            f"No initialization found for {init_date!r} "
            f"(available: {start[0]:%Y-%m} to {start[-1]:%Y-%m})"
        )
    return int(matches[0])


def _available_models(ds, spec, now_idx):
    """Model names with a non-NaN latest-init, zero-lead forecast."""
    ok = ds[spec["var"]].isel(S=now_idx, L=0).mean("M").notnull()
    return ds.model.where(ok, drop=True).values


def _index_transform(da, spec, seasonal, model, start_month):
    """Apply the seasonal rolling-mean and, for scaled indices, the
    (model, start-month, lead)-dependent variance-matching factor (from
    spec["scale"]).

    `da` must have `L` as its only remaining non-(M,) dim at the point of
    the call (e.g. dims (L,) after an S/model/M reduction, or (M, L)
    before the M reduction — the factor depends on L directly, not via
    target month, so it broadcasts either way). Order matches the
    original code: rolling-mean and the factor multiply are both linear in
    L and commute with a subsequent mean("M").
    """
    out = _seasonal(da) if seasonal else da
    scale = spec.get("scale")
    if scale is not None:
        factor_monthly, factor_seasonal = scale
        factor = factor_seasonal if seasonal else factor_monthly
        out = out * factor.sel(model=model, month=start_month)
    return out


def _synthetic_plume(ds, avail, spec, seasonal, now_idx, mmm):
    """Draw N_SYNTHETIC_MEMBERS Gaussian scenarios from the historical MMM
    forecast-error covariance across leads, centered on the current `mmm`.

    Methodology (Barnston, Tippett, van den Dool & Unger 2015, JAMC, Fig. 9
    lower panels; https://doi.org/10.1175/JAMC-D-14-0188.1): the MMM
    forecast error over the 1991-2020 hindcast, stratified by the current
    start month, has a 12x12 lead-by-lead covariance. Its diagonal is the
    per-lead error variance (~SEE^2); its off-diagonals are the lead-to-lead
    error correlation. Drawing from a multivariate normal with this
    covariance and adding to the current MMM produces a plume whose width
    and lead-to-lead coherence reflect actual out-of-sample skill, unlike
    the raw ensemble-member spread (plot_spread), which mixes model-specific
    dispersion with skill and is not calibrated.

    The historical MMM forecast and the error reference are built in the
    same transformed/scaled space as the plotted plume (seasonal rolling
    mean, then — for the n34r index — the (model, month, L) variance-
    matching factor), so the seasonal and n34r covariances are each
    computed in that figure's own final units rather than derived from a
    shared monthly draw. Both n34 and n34r verify against the observed
    *absolute* Niño-3.4 (ds.obsa), consistent with rel_scaling_factor
    (which calibrates the model relative index to observed absolute
    Niño-3.4 variance), not the observed relative index (ds.obsa_rel).

    The error is used demeaned (np.cov subtracts the sample mean per lead):
    any residual MMM conditional bias is intentionally not injected into
    the synthetic members — the plume stays centered on the current MMM.

    Returns
    -------
    xarray.DataArray
        dims (member, L), member = 0..N_SYNTHETIC_MEMBERS-1, L = mmm's full
        lead coordinate (NaN at leads dropped for missing data, e.g. the
        seasonal rolling mean's first/last lead).
    """
    avail_models = avail  # already model names (see _available_models)
    year = ds.S.dt.year
    hindcast = (year >= config.CLIM_START_YEAR) & (year <= config.CLIM_END_YEAR)

    # The hindcast period is fixed to CLIM_START_YEAR-CLIM_END_YEAR (1991-2020)
    # specifically because every NMME model has complete forecast coverage
    # there (verified 360/360 starts per model at zero lead) — mean("model")
    # below is skipna, so a model missing part of this window would silently
    # drop out of some historical starts rather than raising; this guards
    # that assumption instead of assuming it holds forever (e.g. if a future
    # model with a shorter hindcast record is added to the store).
    lead0 = ds[spec["var"]].sel(model=avail_models).isel(L=0).mean("M").where(hindcast, drop=True)
    incomplete = lead0.isnull().any("S")
    if bool(incomplete.any()):
        bad = [str(m) for m in incomplete.where(incomplete, drop=True).model.values]
        raise ValueError(
            f"_synthetic_plume: model(s) {bad} lack complete "
            f"{config.CLIM_START_YEAR}-{config.CLIM_END_YEAR} hindcast coverage; "
            "the historical MMM error assumes every avail model contributes at "
            "every hindcast start"
        )

    fc = ds[spec["var"]].sel(model=avail_models).mean("M").where(hindcast, drop=True)
    obs = ds.obsa.where(hindcast, drop=True)
    if seasonal:
        fc = fc.rolling(L=3, center=True).mean()
        obs = obs.rolling(L=3, center=True).mean()
    scale = spec.get("scale")
    if scale is not None:
        factor_monthly, factor_seasonal = scale
        factor = factor_seasonal if seasonal else factor_monthly
        fc = fc * factor.sel(model=avail_models, month=fc.S.dt.month)
    mmm_hist = fc.mean("model")  # (S, L)

    err = mmm_hist - obs  # (S, L)
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)
    err_month = err.where(err.S.dt.month == start_month_now, drop=True)  # (S~30, L)
    err_valid = err_month.dropna("L", how="any")  # drop leads with any missing sample
    valid_L = err_valid["L"]

    n_dropped = mmm.sizes["L"] - valid_L.sizes["L"]
    expected_dropped = 2 if seasonal else 0
    if n_dropped != expected_dropped:
        print(
            f"  warning: _synthetic_plume ({spec['prefix']}, "
            f"{'seasonal' if seasonal else 'monthly'}) dropped {n_dropped} "
            f"lead(s), expected {expected_dropped}"
        )

    cov = np.cov(err_valid.transpose("S", "L").values, rowvar=False)  # (nL, nL)
    rng = np.random.default_rng(SYNTHETIC_SEED)
    draws = rng.multivariate_normal(np.zeros(cov.shape[0]), cov, size=N_SYNTHETIC_MEMBERS)

    draws_da = xr.DataArray(
        draws,
        dims=("member", "L"),
        coords={"L": valid_L, "member": np.arange(N_SYNTHETIC_MEMBERS)},
    )
    synthetic = draws_da + mmm.sel(L=valid_L)
    return synthetic.reindex(L=mmm["L"])


def _fmt_init(date):
    """Nominal initialization date for display, e.g. 'Aug 1, 2026'."""
    return pd.Timestamp(date).strftime("%b %-d, %Y")


_INIT_BBOX = dict(boxstyle="round", facecolor="white", edgecolor="0.6", alpha=0.85)


def _init_textbox(ax, text):
    """Small boxed annotation (upper-left, axes fraction) carrying the
    initialization date, kept out of ax.set_title() so titles stay short."""
    ax.text(
        0.02, 0.98, text, transform=ax.transAxes, ha="left", va="top",
        fontsize=10, bbox=_INIT_BBOX,
    )


def _place_grid_init(fg, text):
    """Write the init-date box into the grid's first unused facet slot
    (col_wrap leaves one empty whenever the model count isn't a multiple of
    the wrap width); fall back to a suptitle if the grid is exactly full."""
    empty_idxs = [i for i, nd in enumerate(fg.name_dicts.flat) if nd is None]
    if not empty_idxs:
        fg.fig.suptitle(text, x=0.42, y=1.03, fontsize=14, fontweight="bold")
        return
    ax = fg.axs.flat[empty_idxs[0]]
    ax.set_visible(True)
    ax.set_axis_off()
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=13, fontweight="bold", bbox=_INIT_BBOX)


def _out_path(spec, kind, name, date_suffix, ext="png"):
    """Figure/table path: plots/latest_forecast/<prefix>/<kind>/<name><suffix>.<ext>."""
    d = config.PLOTS_DIR_LATEST_FORECAST / spec["prefix"] / kind
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{name}{date_suffix}.{ext}"


def _tight_xlim(ticks):
    """xlim spanning exactly the tick range, padded by half a time step on each side."""
    half_step = (ticks[1] - ticks[0]) / 2
    return (ticks[0] - half_step, ticks[-1] + half_step)


def _set_xaxis(fig, ax, ticks, seasonal):
    """Apply x-tick positions, either date or season-label formatting, and a
    tight xlim so no tick/label falls beyond the plotted data.

    In the seasonal variant, `_seasonal`'s centered rolling mean leaves the
    first and last lead NaN, so those two ticks never have a plotted point;
    drop them before setting ticks/xlim.
    """
    if seasonal:
        ticks = ticks[1:-1]
        ax.set_xticks(ticks)
        ax.set_xticklabels([_season_label(d.month) for d in ticks])
        ax.set_xlabel("Target season")
    else:
        ax.set_xticks(ticks)
        fig.autofmt_xdate()
        ax.set_xlabel("Valid date")
    ax.set_xlim(_tight_xlim(ticks))


def plot_grid(ds, start, spec, now_idx, date_suffix):
    """Per-model facet grid of the latest forecast anomaly (lead x member)."""
    da = ds[spec["var"]].isel(S=now_idx)
    scale = spec.get("scale")
    if scale is not None:
        start_month = int(ds.S.isel(S=now_idx).dt.month)
        da = da * scale[0].sel(month=start_month)  # monthly factor; grid has no seasonal variant

    # ds.ssta/ssta_rel carry a stale standard_name/units inherited from the
    # raw store field (and, for n34r, the multiply can pull in the ERSSTv5
    # scale factor's attrs too) — xarray's auto-labeling prefers long_name/
    # standard_name over the array name, so a leftover attr silently
    # mislabels the colorbar. Replace attrs outright rather than patching
    # them, so both index grids get a clean "prefix [units]" label.
    da = da.rename(spec["prefix"])
    da.attrs = {"units": "degC"}

    fg = da.plot(col="model", col_wrap=4)
    fig = fg.fig
    fig.set_facecolor("white")
    _place_grid_init(fg, f"Init:\n{_fmt_init(start[now_idx])}")

    out = _out_path(spec, "monthly", "grid", date_suffix)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_compare(ds, start, avail, model_colors, spec, now_idx, prev_idx, date_suffix, seasonal=False):
    """Ensemble-mean plume: latest init (solid) vs. previous init (dashed)."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))

    leads_prev = pd.date_range(start[prev_idx], periods=12, freq="MS")
    leads_now = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_prev = int(ds.S.isel(S=prev_idx).dt.month)
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    prev_list, now_list = [], []
    for model in avail:
        color = model_colors[model]

        prev = _index_transform(ds[var].isel(S=prev_idx).mean("M").sel(model=model), spec, seasonal, model, start_month_prev)
        ax.plot(leads_prev, prev, "--", lw=2, color=color, alpha=0.5)
        ax.plot(leads_prev[l0], prev.isel(L=l0), "o", lw=2, color=color, alpha=0.5)
        prev_list.append(prev)

        now = _index_transform(ds[var].isel(S=now_idx).mean("M").sel(model=model), spec, seasonal, model, start_month_now)
        ax.plot(leads_now, now, lw=3, color=color, label=config.short_label(model), alpha=0.75)
        ax.plot(leads_now[l0], now.isel(L=l0), "s", lw=3, color=color, alpha=0.75)
        now_list.append(now)

    mmm_prev = xr.concat(prev_list, dim="model").mean("model")
    mmm_now = xr.concat(now_list, dim="model").mean("model")
    ax.plot(leads_prev, mmm_prev, "--", lw=3, color=MMM_COLOR, alpha=0.9)
    ax.plot(leads_prev[l0], mmm_prev.isel(L=l0), "o", lw=3, color=MMM_COLOR, alpha=0.9)
    ax.plot(leads_now, mmm_now, lw=4, color=MMM_COLOR, label="MMM", alpha=0.9)
    ax.plot(leads_now[l0], mmm_now.isel(L=l0), "s", lw=4, color=MMM_COLOR, alpha=0.9)

    ticks = pd.date_range(start[prev_idx], periods=13, freq="MS")
    _set_xaxis(fig, ax, ticks, seasonal)
    kind_label = "seasonal (3-month running mean)" if seasonal else "monthly"
    ax.set_title(f"{spec['name']} {kind_label} anomaly (1991-2020 climatology mostly)")
    _init_textbox(ax, f"Init: {_fmt_init(start[prev_idx])} (dashed)\n→ {_fmt_init(start[now_idx])} (solid)")
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out_png = _out_path(spec, kind, "compare", date_suffix)
    fig.savefig(out_png, dpi=200, format="png")
    plt.close(fig)
    print(f"  wrote {out_png}")


def plot_spread(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=False):
    """All members (thin) + ensemble mean (thick) per model, latest init only."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for model in avail:
        color = model_colors[model]
        members = _index_transform(ds[var].isel(S=now_idx).sel(model=model), spec, seasonal, model, start_month_now)
        mean = members.mean("M")

        ax.plot(leads, members.T, lw=1.5, color=color, alpha=0.35)
        ax.plot(leads, mean, lw=4, color=color, label=config.short_label(model), alpha=0.85)
        ax.plot(leads[l0], mean.isel(L=l0), "s", lw=3, color=color)
        mean_list.append(mean)

    mmm = xr.concat(mean_list, dim="model").mean("model")
    ax.plot(leads, mmm, lw=5, color=MMM_COLOR, label="MMM", alpha=0.9)
    ax.plot(leads[l0], mmm.isel(L=l0), "s", lw=4, color=MMM_COLOR)

    ticks = pd.date_range(start[now_idx], periods=12, freq="MS")
    _set_xaxis(fig, ax, ticks, seasonal)
    kind_label = "seasonal (3-month running mean)" if seasonal else "monthly"
    ax.set_title(f"NMME forecast {spec['name']} {kind_label} anomaly (1991-2020 climatology)")
    _init_textbox(ax, f"Init: {_fmt_init(start[now_idx])}")
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out = _out_path(spec, kind, "spread", date_suffix)
    fig.savefig(out, dpi=200, format="png")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_spread_synthetic(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=False):
    """Calibrated alternative to plot_spread: N_SYNTHETIC_MEMBERS Gaussian
    scenarios drawn from the historical MMM forecast-error covariance across
    leads, added to the current MMM (see _synthetic_plume), plus a 10th/90th
    percentile envelope. Latest init only."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for model in avail:
        members = _index_transform(ds[var].isel(S=now_idx).sel(model=model), spec, seasonal, model, start_month_now)
        mean_list.append(members.mean("M"))
    mmm = xr.concat(mean_list, dim="model").mean("model")

    synthetic = _synthetic_plume(ds, avail, spec, seasonal, now_idx, mmm)
    pct = synthetic.quantile([0.1, 0.9], dim="member")
    lo, hi = pct.sel(quantile=0.1), pct.sel(quantile=0.9)

    ax.plot(leads, synthetic.transpose("L", "member"), lw=1, color=SYNTHETIC_MEMBER_COLOR, alpha=0.25)
    ax.plot([], [], lw=1.5, color=SYNTHETIC_MEMBER_COLOR, alpha=0.6, label=f"synthetic members (n={N_SYNTHETIC_MEMBERS})")
    ax.plot(leads, lo, "--", lw=2, color=SYNTHETIC_MMM_COLOR, alpha=0.8, label="10th/90th percentile")
    ax.plot(leads, hi, "--", lw=2, color=SYNTHETIC_MMM_COLOR, alpha=0.8)
    ax.plot(leads, mmm, lw=5, color=SYNTHETIC_MMM_COLOR, label="MMM", alpha=0.9)
    ax.plot(leads[l0], mmm.isel(L=l0), "s", lw=4, color=SYNTHETIC_MMM_COLOR)

    ticks = pd.date_range(start[now_idx], periods=12, freq="MS")
    _set_xaxis(fig, ax, ticks, seasonal)
    kind_label = "seasonal (3-month running mean)" if seasonal else "monthly"
    ax.set_title(
        f"NMME forecast {spec['name']} {kind_label} anomaly — synthetic plume "
        f"(MMM + 1991-2020 error covariance)"
    )
    _init_textbox(ax, f"Init: {_fmt_init(start[now_idx])}")
    ax.legend(ncol=1)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out = _out_path(spec, kind, "spread_synthetic", date_suffix)
    fig.savefig(out, dpi=200, format="png")
    plt.close(fig)
    print(f"  wrote {out}")


STRENGTH_RED_EDGE = "#ff0000"
STRENGTH_BLUE_EDGE = "#0000ff"
STRENGTH_NEUTRAL_EDGE = "#7c7c7c"


def _strength_categories():
    """The ENSO strength categories (NOAA CPC ENSO Strength Probabilities
    chart convention — colors sampled from a CPC chart screenshot), split
    into the chart's 3 bar groups. Each entry is (label, lo, hi, lo_closed,
    hi_closed, facecolor); la_nina and el_nino lists are ascending (weakest
    first — the bottom of each stacked bar), neutral is a single category.
    Boundaries are exhaustive and non-overlapping across all three groups —
    El Niño categories are lower-inclusive/upper-exclusive, La Niña
    categories are lower-exclusive/upper-inclusive, Neutral is open on both
    sides — so every threshold value (±0.5, ±1.0, ±1.5, ±2.0, 3.0) belongs
    to exactly one category.

    El Niño adds a Super El Niño category (index >= 3.0) above Very Strong
    (now capped at < 3.0) — a project-specific extension beyond the CPC
    chart, which stops at Very Strong (index >= 2.0); confirmed with the
    user 2026-09-07.
    """
    return {
        "la_nina": [
            ("Weak La Niña\n-1.0°C < index ≤ -0.5°C", -1.0, -0.5, False, True, "#dbeaff"),
            ("Moderate La Niña\n-1.5°C < index ≤ -1.0°C", -1.5, -1.0, False, True, "#9fc2ff"),
            ("Strong La Niña\n-2.0°C < index ≤ -1.5°C", -2.0, -1.5, False, True, "#4d88ff"),
            ("Very Strong La Niña\nindex ≤ -2.0°C", -np.inf, -2.0, False, True, "#0033cc"),
        ],
        "neutral": ("Neutral\n-0.5°C < index < 0.5°C", -0.5, 0.5, False, False, "#d3d3d3"),
        "el_nino": [
            ("Weak El Niño\n0.5°C ≤ index < 1.0°C", 0.5, 1.0, True, False, "#ffe5e5"),
            ("Moderate El Niño\n1.0°C ≤ index < 1.5°C", 1.0, 1.5, True, False, "#ffb3b3"),
            ("Strong El Niño\n1.5°C ≤ index < 2.0°C", 1.5, 2.0, True, False, "#ff6666"),
            ("Very Strong El Niño\n2.0°C ≤ index < 3.0°C", 2.0, 3.0, True, False, "#990000"),
            ("Super El Niño\nindex ≥ 3.0°C", 3.0, np.inf, True, False, "#660000"),
        ],
    }


def _in_category(x, lo, hi, lo_closed, hi_closed):
    left = x >= lo if lo_closed else x > lo
    right = x <= hi if hi_closed else x < hi
    return left & right


def plot_strength_probabilities(ds, start, avail, spec, now_idx, date_suffix):
    """Grouped/stacked-bar ENSO strength-category probability by target
    season — n34r/seasonal only, format after NOAA CPC's ENSO Strength
    Probabilities chart (see module docstring): 3 bars per season (La Niña,
    Neutral, El Niño), the La Niña and El Niño bars internally stacked by
    strength category. Category probabilities are the empirical fraction of
    _synthetic_plume draws (at N_STRENGTH_DRAWS, not the plotting-default
    N_SYNTHETIC_MEMBERS) falling in each _strength_categories bucket, per
    season.
    """
    var = spec["var"]
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for model in avail:
        members = _index_transform(ds[var].isel(S=now_idx).sel(model=model), spec, True, model, start_month_now)
        mean_list.append(members.mean("M"))
    mmm = xr.concat(mean_list, dim="model").mean("model")

    global N_SYNTHETIC_MEMBERS
    orig_n = N_SYNTHETIC_MEMBERS
    N_SYNTHETIC_MEMBERS = N_STRENGTH_DRAWS
    try:
        synthetic = _synthetic_plume(ds, avail, spec, True, now_idx, mmm)  # (member, L)
    finally:
        N_SYNTHETIC_MEMBERS = orig_n

    # Seasonal rolling mean leaves the first/last lead NaN (see _set_xaxis);
    # _synthetic_plume's expected_dropped=2 check already relies on that
    # being exactly the first and last lead (as does write_summary_tables).
    leads_seasonal = leads[1:-1]
    valid = synthetic.isel(L=slice(1, -1))  # (member, L=10)

    categories = _strength_categories()
    x = np.arange(len(leads_seasonal))
    bar_width, offset = 0.26, 0.27

    fig, ax = plt.subplots(figsize=(13, 7))

    bottom = np.zeros(len(leads_seasonal))
    for label, lo, hi, lo_closed, hi_closed, color in categories["la_nina"]:
        pct = (_in_category(valid, lo, hi, lo_closed, hi_closed).mean("member") * 100).values
        ax.bar(x - offset, pct, bottom=bottom, width=bar_width, color=color, edgecolor=STRENGTH_BLUE_EDGE, linewidth=1)
        bottom += pct

    label, lo, hi, lo_closed, hi_closed, color = categories["neutral"]
    pct = (_in_category(valid, lo, hi, lo_closed, hi_closed).mean("member") * 100).values
    ax.bar(x, pct, width=bar_width, color=color, edgecolor=STRENGTH_NEUTRAL_EDGE, linewidth=1)

    bottom = np.zeros(len(leads_seasonal))
    for label, lo, hi, lo_closed, hi_closed, color in categories["el_nino"]:
        pct = (_in_category(valid, lo, hi, lo_closed, hi_closed).mean("member") * 100).values
        ax.bar(x + offset, pct, bottom=bottom, width=bar_width, color=color, edgecolor=STRENGTH_RED_EDGE, linewidth=1)
        bottom += pct

    ax.set_xticks(x)
    ax.set_xticklabels([_season_label(d.month) for d in leads_seasonal])
    ax.set_xlabel("Season")
    ax.set_ylabel("Percent Chance (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.grid(axis="y", visible=True, alpha=0.3)

    fig.suptitle(f"ENSO Strength Probabilities (issued {start[now_idx]:%B %Y})", x=0.44, y=1.0, fontsize=18, fontweight="bold")
    ax.set_title("Based on the NMME MMM and historical performance", fontsize=12)

    # Legend order matches the CPC chart: El Niño strongest-to-weakest, then
    # Neutral, then La Niña weakest-to-strongest — built explicitly (rather
    # than from plotting order) since only the El Niño side needs reversing.
    legend_entries = (
        [(lbl, c, STRENGTH_RED_EDGE) for lbl, *_, c in reversed(categories["el_nino"])]
        + [(categories["neutral"][0], categories["neutral"][-1], STRENGTH_NEUTRAL_EDGE)]
        + [(lbl, c, STRENGTH_BLUE_EDGE) for lbl, *_, c in categories["la_nina"]]
    )
    handles = [
        matplotlib.patches.Patch(facecolor=c, edgecolor=e, label=lbl, linewidth=1)
        for lbl, c, e in legend_entries
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.02), fontsize=8, frameon=True)

    fig.set_facecolor("white")
    plt.tight_layout()

    out = _out_path(spec, "seasonal", "strength_probabilities", date_suffix)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def _historical_mmm(ds, avail, spec, seasonal):
    """Fixed-model-set MMM forecast anomaly for every start S from
    config.ANALYSIS_START_YEAR to the present.

    Used as the ranking pool for the latest-forecast summary tables
    (write_summary_tables). `avail` is the model set available at the
    *current* latest forecast (see _available_models); using that same
    fixed set across all of history keeps this MMM numerically identical
    to the one already plotted in plot_mean/plot_compare for the current
    forecast (rather than reproducing each historical year's actual,
    varying NMME model roster).

    Returns
    -------
    xarray.DataArray
        dims (S, L), restricted to S.dt.year >= config.ANALYSIS_START_YEAR.
    """
    fc = ds[spec["var"]].sel(model=avail).mean("M")  # (model, S, L)
    if seasonal:
        fc = fc.rolling(L=3, center=True).mean()
    scale = spec.get("scale")
    if scale is not None:
        factor_monthly, factor_seasonal = scale
        factor = factor_seasonal if seasonal else factor_monthly
        fc = fc * factor.sel(model=avail, month=fc.S.dt.month)
    mmm = fc.mean("model")  # (S, L)
    year = ds.S.dt.year
    return mmm.where(year >= config.ANALYSIS_START_YEAR, drop=True)


def _rank_at_lead(mmm_hist, start_month, current_S):
    """Rank (1 = highest) of the current forecast's anomaly at each lead
    among all historical MMM forecasts issued in the same calendar start
    month, config.ANALYSIS_START_YEAR-present.

    Returns
    -------
    tuple[xarray.DataArray, xarray.DataArray, int]
        (current values, ranks), each dims (L,) — rank is NaN wherever the
        current value itself is NaN (e.g. seasonal running-mean endpoints)
        — and the ranking pool size (number of years).
    """
    pool = mmm_hist.where(mmm_hist.S.dt.month == start_month, drop=True)  # (S', L)
    current = pool.sel(S=current_S)
    rank = (pool > current).sum("S") + 1
    rank = rank.where(current.notnull())
    return current, rank, pool.sizes["S"]


def write_summary_tables(ds, start, now_idx, index_specs, avail_by_prefix, date_suffix):
    """Markdown table of the latest-forecast MMM anomaly (monthly and
    seasonal target periods), with each value's rank (1 = highest) among
    all MMM forecasts issued in the same calendar start month,
    config.ANALYSIS_START_YEAR-present — see specs/latest_forecast.md §5.
    """
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)
    current_S = ds.S.isel(S=now_idx)
    start_date = start[now_idx]
    n_leads = ds.sizes["L"]

    leads_monthly = pd.date_range(start_date, periods=n_leads, freq="MS")
    leads_seasonal = leads_monthly[1:-1]

    lines = [
        f"# Latest-forecast MMM summary — start {start_date:%Y-%m}",
        "",
        f"Multi-model-mean (MMM) anomaly for the {n_leads}-lead forecast "
        f"initialized {start_date:%B %Y}, for the standard Niño-3.4 index "
        f"(n34) and the relative Niño-3.4 index (n34r). Below each anomaly, "
        f"its rank (1 = highest) among all MMM forecasts issued in "
        f"{start_date:%B}, {config.ANALYSIS_START_YEAR}-present, same "
        f"target period. MMM uses the fixed model set available for the "
        f"current forecast, applied consistently across all ranked years "
        f"(same MMM definition as the plume figures).",
        "",
    ]

    for seasonal, title in ((False, "Monthly"), (True, "Seasonal (3-month running mean)")):
        col_labels = (
            [_season_label(d.month) for d in leads_seasonal]
            if seasonal
            else [f"{d:%b %Y}" for d in leads_monthly]
        )

        rows = []
        for spec in index_specs:
            avail = avail_by_prefix[spec["prefix"]]
            mmm_hist = _historical_mmm(ds, avail, spec, seasonal)
            current, rank, pool_size = _rank_at_lead(mmm_hist, start_month_now, current_S)
            if seasonal:
                current = current.isel(L=slice(1, -1))
                rank = rank.isel(L=slice(1, -1))
            anom_row = ["–" if not np.isfinite(v) else f"{v:.2f}" for v in current.values]
            rank_row = ["–" if not np.isfinite(r) else f"{int(r)}" for r in rank.values]
            rows.append((f"{spec['prefix']} anom", anom_row))
            rows.append((f"{spec['prefix']} rank (n={pool_size})", rank_row))

        lines.append(f"## {title}")
        lines.append("")
        lines.append("| | " + " | ".join(col_labels) + " |")
        lines.append("|---|" + "---|" * len(col_labels))
        for label, vals in rows:
            lines.append(f"| {label} | " + " | ".join(vals) + " |")
        lines.append("")

    out = config.PLOTS_DIR_LATEST_FORECAST / f"latest_forecast_summary{date_suffix}.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out}")


def plot_mean(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=False):
    """Ensemble-mean-only plume, latest init."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for model in avail:
        color = model_colors[model]
        mean = _index_transform(ds[var].isel(S=now_idx).mean("M").sel(model=model), spec, seasonal, model, start_month_now)
        ax.plot(leads, mean, lw=3, color=color, label=config.short_label(model), alpha=0.75)
        ax.plot(leads[l0], mean.isel(L=l0), "s", lw=3, color=color, alpha=0.75)
        mean_list.append(mean)

    mmm = xr.concat(mean_list, dim="model").mean("model")
    ax.plot(leads, mmm, lw=4, color=MMM_COLOR, label="MMM", alpha=0.9)
    ax.plot(leads[l0], mmm.isel(L=l0), "s", lw=4, color=MMM_COLOR)

    ticks = pd.date_range(start[now_idx], periods=12, freq="MS")
    _set_xaxis(fig, ax, ticks, seasonal)
    kind_label = "seasonal (3-month running mean)" if seasonal else "monthly"
    ax.set_title(f"{spec['name']} {kind_label} anomaly (1991-2020 climatology)")
    _init_textbox(ax, f"Init: {_fmt_init(start[now_idx])}")
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out = _out_path(spec, kind, "mean", date_suffix)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--init-date",
        default=None,
        help="Plot this initialization (YYYY-MM or YYYY-MM-DD), matched by calendar "
        "month, instead of the latest available. Output filenames get a "
        "_<YYYY-MM> suffix.",
    )
    args = parser.parse_args()

    ds = config.load_nino34_verification()
    start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
    now_idx = _resolve_init_idx(start, args.init_date)
    prev_idx = now_idx - 1
    date_suffix = f"_{start[now_idx]:%Y-%m}" if args.init_date is not None else ""
    colors = _plume_colors()
    # Keyed by model name (not position) so color assignment survives the
    # per-index_spec `avail` subset/order — see _available_models.
    model_colors = {name: colors[i] for i, name in enumerate(ds.model.values)}
    factor_monthly, factor_seasonal = config.rel_scaling_factor(ds)

    config.PLOTS_DIR_LATEST_FORECAST.mkdir(parents=True, exist_ok=True)

    index_specs = [
        {"var": "ssta", "prefix": "n34", "name": "Nino 3.4", "scale": None},
        {
            "var": "ssta_rel",
            "prefix": "n34r",
            "name": "relative Nino 3.4",
            "scale": (factor_monthly, factor_seasonal),
        },
    ]

    avail_by_prefix = {}
    for spec in index_specs:
        avail = _available_models(ds, spec, now_idx)
        avail_by_prefix[spec["prefix"]] = avail
        print(
            f"  available models (init {start[now_idx]:%Y-%m}, {spec['prefix']}): "
            f"{[config.short_label(m) for m in avail]}"
        )

        plot_grid(ds, start, spec, now_idx, date_suffix)
        for seasonal in (False, True):
            plot_compare(ds, start, avail, model_colors, spec, now_idx, prev_idx, date_suffix, seasonal=seasonal)
            plot_spread(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=seasonal)
            plot_spread_synthetic(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=seasonal)
            plot_mean(ds, start, avail, model_colors, spec, now_idx, date_suffix, seasonal=seasonal)

        if spec["prefix"] == "n34r":
            plot_strength_probabilities(ds, start, avail, spec, now_idx, date_suffix)

    write_summary_tables(ds, start, now_idx, index_specs, avail_by_prefix, date_suffix)


if __name__ == "__main__":
    main()
