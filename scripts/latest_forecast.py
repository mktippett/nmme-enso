"""
latest_forecast.py — NMME Niño-3.4 latest-forecast plume figures.

Reads the local NMME SST zarr store (config.STORE_SST) via
config.load_nino34_verification(), and plots the most recent forecast
initialization (and the previous one, for comparison) as a set of
monitoring figures. These are exploratory/monitoring figures, not
manuscript figures: they carry the climatology note in the title rather
than travelling with a separate caption.

Two indices are plotted, each as its own figure set (n34_* / n34r_*):
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
import matplotlib.pyplot as plt

import config

MMM_COLOR = "0.75"

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
    """Positional model indices with a non-NaN latest-init, zero-lead forecast."""
    return np.where(~np.isnan(ds[spec["var"]].isel(S=now_idx, L=0).mean("M")))[0]


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


def plot_grid(ds, spec, now_idx, date_suffix):
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
    start_date = str(ds.S.isel(S=now_idx).values)[:10]
    fig.suptitle(f"Start = {start_date}", x=0.42, y=1.03, fontsize=14, fontweight="bold")

    out = config.PLOTS_DIR_LATEST_FORECAST / f"{spec['prefix']}_monthly_grid{date_suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_compare(ds, start, avail, plume_colors, spec, now_idx, prev_idx, date_suffix, seasonal=False):
    """Ensemble-mean plume: latest init (solid) vs. previous init (dashed)."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))

    leads_prev = pd.date_range(start[prev_idx], periods=12, freq="MS")
    leads_now = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_prev = int(ds.S.isel(S=prev_idx).dt.month)
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    prev_list, now_list = [], []
    for im in avail:
        model = ds.model.isel(model=im).item()
        color = plume_colors[im]

        prev = _index_transform(ds[var].isel(S=prev_idx).mean("M").isel(model=im), spec, seasonal, model, start_month_prev)
        ax.plot(leads_prev, prev, "--", lw=2, color=color, alpha=0.5)
        ax.plot(leads_prev[l0], prev.isel(L=l0), "o", lw=2, color=color, alpha=0.5)
        prev_list.append(prev)

        now = _index_transform(ds[var].isel(S=now_idx).mean("M").isel(model=im), spec, seasonal, model, start_month_now)
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
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out_png = config.PLOTS_DIR_LATEST_FORECAST / f"{spec['prefix']}_{kind}_compare{date_suffix}.png"
    fig.savefig(out_png, dpi=200, format="png")
    plt.close(fig)
    print(f"  wrote {out_png}")


def plot_spread(ds, start, avail, colors, spec, now_idx, date_suffix, seasonal=False):
    """All members (thin) + ensemble mean (thick) per model, latest init only."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for im in avail:
        model = ds.model.isel(model=im).item()
        color = colors[im]
        members = _index_transform(ds[var].isel(S=now_idx).isel(model=im), spec, seasonal, model, start_month_now)
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
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")
    plt.tight_layout()

    kind = "seasonal" if seasonal else "monthly"
    out = config.PLOTS_DIR_LATEST_FORECAST / f"{spec['prefix']}_{kind}_spread{date_suffix}.png"
    fig.savefig(out, dpi=200, format="png")
    plt.close(fig)
    print(f"  wrote {out}")


def plot_mean(ds, start, avail, colors, spec, now_idx, date_suffix, seasonal=False):
    """Ensemble-mean-only plume, latest init."""
    l0 = 1 if seasonal else 0
    var = spec["var"]

    fig, ax = plt.subplots(figsize=(12, 7))
    leads = pd.date_range(start[now_idx], periods=12, freq="MS")
    start_month_now = int(ds.S.isel(S=now_idx).dt.month)

    mean_list = []
    for im in avail:
        model = ds.model.isel(model=im).item()
        color = colors[im]
        mean = _index_transform(ds[var].isel(S=now_idx).mean("M").isel(model=im), spec, seasonal, model, start_month_now)
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
    ax.legend(ncol=2)
    ax.grid(visible=True)
    fig.set_facecolor("white")

    kind = "seasonal" if seasonal else "monthly"
    out = config.PLOTS_DIR_LATEST_FORECAST / f"{spec['prefix']}_{kind}_mean{date_suffix}.png"
    fig.savefig(out)
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

    for spec in index_specs:
        avail = _available_models(ds, spec, now_idx)
        print(
            f"  available models (init {start[now_idx]:%Y-%m}, {spec['prefix']}): "
            f"{[config.short_label(ds.model.isel(model=i).item()) for i in avail]}"
        )

        plot_grid(ds, spec, now_idx, date_suffix)
        for seasonal in (False, True):
            plot_compare(ds, start, avail, colors, spec, now_idx, prev_idx, date_suffix, seasonal=seasonal)
            plot_spread(ds, start, avail, colors, spec, now_idx, date_suffix, seasonal=seasonal)
            plot_mean(ds, start, avail, colors, spec, now_idx, date_suffix, seasonal=seasonal)


if __name__ == "__main__":
    main()
