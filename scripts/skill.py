"""
skill.py — NMME Niño-3.4 forecast-skill heatmaps (correlation, MSESS).

Reads config.load_nino34_verification() (NMME ensemble-mean Niño-3.4
forecast anomalies aligned against ERSSTv5 observational anomalies on the
same (S, L) grid — see that function's docstring for the alignment
approach) and computes two skill scores, each stratified by calendar month
and lead, restricted to starts in CLIM_START_YEAR-CLIM_END_YEAR
(1991-2020):

  - Anomaly correlation (AC): per-(month, lead) Pearson correlation between
    the ensemble-mean forecast and ERSSTv5.
  - MSESS: 1 - MSE(forecast)/MSE(zero-anomaly reference), per (month, lead).

Both are computed twice: once stratified by *start* month (the month the
forecast was initialized) and once relabeled to *target* month (the month
being predicted, = start + lead). A multi-model-mean (MMM) panel — the
ensemble mean averaged again across models — is added as an 8th panel so
the facet grid is a clean 2x4 rather than 7-with-a-blank-slot, and lets the
MMM's skill be compared directly against each individual model's.

These are exploratory/diagnostic figures (no manuscript caption), so the
calculation variant (start vs. target month, 1991-2020 climatology, MMM
definition) is baked into each suptitle rather than deferred to a caption.

Run:
    mamba run -n pangeo-local python scripts/skill.py
"""
import calendar

import matplotlib
matplotlib.use("Agg")

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import config

MONTH_ABBR = [calendar.month_abbr[m] for m in range(1, 13)]


# ---------------------------------------------------------------------------
# Skill functions (ported from the user's notebook; not modified)
# ---------------------------------------------------------------------------

def ac_by_start(x, y):
    """Anomaly correlation between x and y, stratified by S.month.

    Computes per-start-month means over values where both x and y are not
    missing (masking is needed since models have different lead counts /
    forecast lengths), then the Pearson correlation per start month.
    """
    ok = ~np.isnan(x) & ~np.isnan(y)
    xa = x.where(ok).groupby("S.month") - x.where(ok).groupby("S.month").mean("S")
    ya = y.where(ok).groupby("S.month") - y.where(ok).groupby("S.month").mean("S")
    c = (xa * ya).groupby("S.month").mean("S") / xa.groupby("S.month").std("S") / ya.groupby("S.month").std("S")
    c.attrs["long_name"] = "correlation"
    c.month.attrs["long_name"] = "start month"
    return c


def msess_by_start(f, o):
    """MSESS skill score (vs. a zero-anomaly reference forecast), stratified
    by S.month.

    Assumes f and o are already anomalies (no mean removal here). Uses
    .mean("S") rather than .sum() so a (month, lead) cell is NaN only if
    all its samples are missing (e.g. a model with a shorter lead count),
    not silently treated as zero.
    """
    msess = 1.0 - ((f - o) ** 2).groupby("S.month").mean("S") / (o ** 2).groupby("S.month").mean("S")
    msess.attrs["long_name"] = "MSESS"
    msess.month.attrs["long_name"] = "start month"
    return msess


def to_target_month(sk):
    """Relabel a start-month-stratified skill array to target month
    (= start month + lead), by rolling each lead's month axis.

    Cleaner than an in-place per-lead .loc/roll loop: builds each lead's
    relabeled slice functionally and concatenates, rather than mutating a
    copy in place. sk.L carries half-integer leads (0.5, 1.5, ...); L=0.5
    is the init month itself (shift 0), L=1.5 shifts by 1 month, etc.
    """
    parts = [
        sk.sel(L=l).roll(month=int(round(float(l) - 0.5)), roll_coords=False)
        for l in sk.L.values
    ]
    out = xr.concat(parts, dim="L")
    out.month.attrs["long_name"] = "target month"
    return out


# ---------------------------------------------------------------------------
# Multi-model mean
# ---------------------------------------------------------------------------

def add_mmm(skill_by_model, forecast_mean, obs, skill_fn):
    """Append a multi-model-mean (MMM) panel to a per-model skill array.

    forecast_mean is the ensemble-mean forecast still carrying a model dim
    (dims model, S, L); MMM further averages across model (equal-weight)
    before computing its own skill, tagged model="MMM", and concatenated
    onto skill_by_model along "model".
    """
    mmm_forecast = forecast_mean.mean("model")
    mmm_skill = skill_fn(mmm_forecast, obs).expand_dims(model=["MMM"])
    return xr.concat([skill_by_model, mmm_skill], dim="model")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _display_labels(da):
    """Model coord relabeled to short display names ('MMM' passes through)."""
    labels = [m if m == "MMM" else config.short_label(m) for m in da.model.values]
    return da.assign_coords(model=labels)


def plot_skill_grid(da, suptitle, cbar_label, out_path, vmin=None, vmax=None, cmap=None):
    """8-panel (2x4) facet heatmap: month (x) by lead (y), one panel per
    model plus MMM.

    cmap is passed through because xarray's automatic diverging-colormap
    detection (data straddles zero -> RdBu_r) is disabled whenever vmin and
    vmax are BOTH given explicitly (center=0 alone does not restore it,
    verified empirically) — pass cmap="RdBu_r" explicitly for a difference
    plot with a clipped vmin/vmax.
    """
    da = _display_labels(da)
    fg = da.plot(
        col="model", col_wrap=4, x="month",
        vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={"label": cbar_label},
    )
    fig = fg.fig
    for ax in fg.axs.flat:
        ax.set_xticks(np.arange(1, 13))
        ax.set_xticklabels(MONTH_ABBR, rotation=90)
    fg.set_titles(template="{value}")
    fg.set_ylabels("lead (months)")
    fig.suptitle(suptitle, y=1.03, fontsize=14, fontweight="bold")
    fig.set_facecolor("white")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    config.PLOTS_DIR_SKILL.mkdir(parents=True, exist_ok=True)

    ds = config.load_nino34_verification()
    period = f"{config.CLIM_START_YEAR}-01-01", f"{config.CLIM_END_YEAR}-12-01"
    ds_p = ds.sel(S=slice(*period))
    print(f"  skill.py: restricted to S={period[0]}..{period[1]}, shape={dict(ds_p.ssta.sizes)}")

    forecast_mean = ds_p.ssta.mean("M")  # dims (model, S, L)
    obs = ds_p.obsa                      # dims (S, L)

    ac = ac_by_start(forecast_mean, obs)
    ac = add_mmm(ac, forecast_mean, obs, ac_by_start)

    msess = msess_by_start(forecast_mean, obs)
    msess = add_mmm(msess, forecast_mean, obs, msess_by_start)

    ac_target = to_target_month(ac)
    msess_target = to_target_month(msess)

    period_label = f"{config.CLIM_START_YEAR}-{config.CLIM_END_YEAR}"

    plot_skill_grid(
        ac, f"Niño-3.4 anomaly correlation by start month ({period_label})",
        "correlation", config.PLOTS_DIR_SKILL / "n34_ac_start.png", vmin=0.0, vmax=1.0,
    )
    plot_skill_grid(
        ac_target, f"Niño-3.4 anomaly correlation by target month ({period_label})",
        "correlation", config.PLOTS_DIR_SKILL / "n34_ac_target.png", vmin=0.0, vmax=1.0,
    )
    plot_skill_grid(
        msess, f"Niño-3.4 MSESS by start month ({period_label})",
        "MSESS", config.PLOTS_DIR_SKILL / "n34_msess_start.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_target, f"Niño-3.4 MSESS by target month ({period_label})",
        "MSESS", config.PLOTS_DIR_SKILL / "n34_msess_target.png", vmax=1.0,
    )


if __name__ == "__main__":
    main()
