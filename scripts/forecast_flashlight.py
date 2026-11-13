"""
forecast_flashlight.py — animated n34r seasonal forecast "flashlight" GIF.

One frame per NMME initialization (default Nov 2025-Oct 2026). Axes are fixed in
calendar time and value across all frames, so the observed path stays put and
each month's forecast is a beam cast forward from its own init date — the
question the animation asks is whether the beam follows the path.

  - Observed path (fixed, every frame): ERSSTv5 relative Niño-3.4, 3-month
    running mean, scaled per center calendar month by the obs-only std ratio
    std(n34_anom)/std(rel_anom) over 1991-2020 (L'Heureux et al. 2024's
    RONI-style scaling; same obs target as rel_scaling_compare.py, but on the
    seasonal series).
  - Beam (current init only): n34r seasonal MMM (latest_forecast's
    model-relative scaled ssta_rel, 3-month running mean over lead) with
    5/15/25/75/85/95% bands. Bands are the analytic quantiles of
    Normal(MMM, sqrt(diag cov)), cov = latest_forecast._lead_error_cov's
    1991-2020 MMM error covariance for that init's start month — i.e. the
    marginals of the synthetic plume (Barnston et al. 2015), computed in
    closed form as in plot_strength_probabilities rather than by sampling.
  - Trail (--trail only): earlier inits' MMM lines and 5-95% bands, faded
    gray (current beam drawn on top), so past
    beams stay visible as the animation advances. Default is pure flashlight
    (no trail); the trail variant writes to a separate "_trail" file.

Exploratory figure (no caption travels with it): method details are in the
title/subtitle.

Run:
    mamba run -n pangeo-local python scripts/forecast_flashlight.py
    mamba run -n pangeo-local python scripts/forecast_flashlight.py --trail
    mamba run -n pangeo-local python scripts/forecast_flashlight.py --first-init 2025-06 --last-init 2026-10
"""
import argparse

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import norm

import config
from latest_forecast import (
    _available_models,
    _index_transform,
    _lead_error_cov,
    _resolve_init_idx,
    _season_label,
)

# Beam bands, (lower, upper) quantile pairs, outermost first so inner bands draw on top.
BAND_QUANTILES = [(0.05, 0.95), (0.15, 0.85), (0.25, 0.75)]
BAND_COLOR = "royalblue"
BAND_ALPHAS = [0.15, 0.30, 0.50]
MEAN_COLOR = "navy"
OBS_COLOR = "black"
TRAIL_COLOR = "0.5"  # --trail: earlier inits, gray
TRAIL_BAND_ALPHA = 0.11  # earlier inits' 5-95% band; overlaps accumulate
TRAIL_LINE_ALPHA = 0.5  # earlier inits' MMM line
TRAIL_LW = 1.2
FRAME_MS = 1000  # per-frame duration in the GIF
LAST_FRAME_HOLD = 3  # last frame shown this many times longer


def obs_n34r_seasonal(period=(config.CLIM_START_YEAR, config.CLIM_END_YEAR)):
    """Observed scaled relative Niño-3.4, 3-month running mean, monthly time
    axis labeled by the center month. Scale = std(n34_anom)/std(rel_anom) per
    center calendar month over `period`, both computed on the seasonal series."""
    obs = config._obs_index_anomalies().rolling(time=3, center=True).mean()
    clim = obs.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))
    scale = clim.n34_anom.groupby("time.month").std("time") / clim.rel_anom.groupby("time.month").std("time")
    out = (obs.rel_anom.groupby("time.month") * scale).drop_vars("month", errors="ignore").dropna("time")
    print(
        f"  obs n34r seasonal: {out.time.dt.strftime('%Y-%m').values[0]} to "
        f"{out.time.dt.strftime('%Y-%m').values[-1]}, scale range "
        f"{float(scale.min()):.2f}-{float(scale.max()):.2f}"
    )
    return out


def beam(ds, start, spec, now_idx):
    """Seasonal n34r MMM and band quantiles for one init, on valid leads only.

    Returns a Dataset with dims (L,) and a `time` coord (target center month):
    mean, sd (hindcast-error sd per lead), plus q<pct> variables for each
    BAND_QUANTILES bound.
    """
    avail = _available_models(ds, spec, now_idx)
    start_month = int(ds.S.isel(S=now_idx).dt.month)
    mmm = xr.concat(
        [
            _index_transform(ds[spec["var"]].isel(S=now_idx).sel(model=m), spec, True, m, start_month).mean("M")
            for m in avail
        ],
        dim="model",
    ).mean("model")

    valid_L, cov = _lead_error_cov(ds, avail, spec, True, now_idx)
    sd = xr.DataArray(np.sqrt(np.diag(cov)), dims="L", coords={"L": valid_L})
    mean = mmm.sel(L=valid_L)

    out = xr.Dataset({"mean": mean, "sd": sd})
    for q in sorted({q for pair in BAND_QUANTILES for q in pair}):
        out[f"q{round(q * 100):02d}"] = mean + norm.ppf(q) * sd
    # NMME L labels are 0.5, 1.5, ...; int() gives the whole-month offset of the target center.
    target = pd.DatetimeIndex([start[now_idx] + pd.DateOffset(months=int(lead)) for lead in valid_L.values])
    out = out.assign_coords(time=("L", target.to_period("M").to_timestamp()))
    print(
        f"  init {start[now_idx]:%Y-%m}: models {[config.short_label(m) for m in avail]}, "
        f"targets {target[0]:%Y-%m} to {target[-1]:%Y-%m}"
    )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--first-init", default="2025-11", help="First initialization (YYYY-MM).")
    parser.add_argument("--last-init", default="2026-11", help="Last initialization (YYYY-MM).")
    parser.add_argument("--trail", action="store_true", help="Keep earlier inits' MMM lines and 5-95%% bands as faded gray traces.")
    args = parser.parse_args()

    ds = config.load_nino34_verification()
    start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
    i0 = _resolve_init_idx(start, args.first_init)
    i1 = _resolve_init_idx(start, args.last_init)
    factor_monthly, factor_seasonal = config.rel_scaling_factor(ds)
    spec = {"var": "ssta_rel", "prefix": "n34r", "name": "relative Nino 3.4", "scale": (factor_monthly, factor_seasonal)}

    beams = [beam(ds, start, spec, i) for i in range(i0, i1 + 1)]
    obs = obs_n34r_seasonal()

    # Fixed axes: x from one month before the first beam's first target (so
    # the first init's own season is visible) to the last beam's last target;
    # y covers every frame's outer band and the obs within that window, in 0.5 steps.
    x_first = pd.Timestamp(beams[0].time.values[0]) - pd.DateOffset(months=1)
    x_last = pd.Timestamp(beams[-1].time.values[-1])
    obs_win = obs.sel(time=slice(x_first, x_last))
    lo_q, hi_q = BAND_QUANTILES[0]
    y_lo = min(float(obs_win.min()), *(float(b[f"q{round(lo_q * 100):02d}"].min()) for b in beams))
    y_hi = max(float(obs_win.max()), *(float(b[f"q{round(hi_q * 100):02d}"].max()) for b in beams))
    ylim = (np.floor(y_lo * 2) / 2, np.ceil(y_hi * 2) / 2)
    ticks = pd.date_range(x_first, x_last, freq="MS")
    half = pd.Timedelta(days=15)

    fig, ax = plt.subplots(figsize=(13, 7))
    fig.set_facecolor("white")

    def draw(k):
        ax.clear()
        b = beams[k]
        t = b.time.values
        if args.trail:
            lo_q, hi_q = BAND_QUANTILES[0]
            for j, prev in enumerate(beams[:k]):
                ax.fill_between(
                    prev.time.values, prev[f"q{round(lo_q * 100):02d}"], prev[f"q{round(hi_q * 100):02d}"],
                    color=TRAIL_COLOR, alpha=TRAIL_BAND_ALPHA, lw=0,
                )
                ax.plot(
                    prev.time.values, prev["mean"], color=TRAIL_COLOR, lw=TRAIL_LW, alpha=TRAIL_LINE_ALPHA,
                    label="Earlier forecasts (MMM, 5-95%)" if j == 0 else None,
                )
        for (lo, hi), alpha in zip(BAND_QUANTILES, BAND_ALPHAS):
            ax.fill_between(
                t, b[f"q{round(lo * 100):02d}"], b[f"q{round(hi * 100):02d}"],
                color=BAND_COLOR, alpha=alpha, lw=0,
                label=f"{round(lo * 100)}-{round(hi * 100)}%",
            )
        ax.plot(t, b["mean"], color=MEAN_COLOR, lw=2, label="NMME MMM")
        ax.plot(obs_win.time, obs_win, "o-", color=OBS_COLOR, lw=3.5, ms=6, label="Observed (ERSSTv5)")
        init = start[i0 + k]
        ax.axvline(pd.Timestamp(init.year, init.month, 1), color="0.4", ls=":", lw=1.5)

        ax.axhline(0, color="0.3", lw=0.8)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{_season_label(d.month)}\n{d.year}" if d.month in (1, 7) or d == ticks[0] else _season_label(d.month) for d in ticks])
        ax.set_xlim(ticks[0] - half, ticks[-1] + half)
        ax.set_ylim(ylim)
        ax.set_yticks(np.arange(ylim[0], ylim[1] + 0.25, 0.5))
        ax.set_xlabel("Target season")
        ax.set_ylabel("n34r (°C)")
        ax.grid(axis="y", ls="--", color="0.7")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(
            "NMME relative Niño-3.4 seasonal forecast vs. observed\n"
            "Bands: Gaussian MMM ± 1991-2020 hindcast error (start-month-stratified); obs scaled by 1991-2020 std ratio",
            fontsize=11,
        )
        ax.text(
            0.02, 0.97, f"Issued {init:%b %Y}", transform=ax.transAxes, ha="left", va="top",
            fontsize=16, fontweight="bold",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.6", alpha=0.85),
        )
        ax.legend(loc="upper right", fontsize=10)
        fig.tight_layout()

    frames = list(range(len(beams))) + [len(beams) - 1] * (LAST_FRAME_HOLD - 1)
    anim = FuncAnimation(fig, draw, frames=frames)
    config.PLOTS_DIR_FORECAST_FLASHLIGHT.mkdir(parents=True, exist_ok=True)
    out = config.PLOTS_DIR_FORECAST_FLASHLIGHT / f"n34r_seasonal_{start[i0]:%Y-%m}_{start[i1]:%Y-%m}{'_trail' if args.trail else ''}.gif"
    anim.save(out, writer=PillowWriter(fps=1000 / FRAME_MS), dpi=100)
    plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
