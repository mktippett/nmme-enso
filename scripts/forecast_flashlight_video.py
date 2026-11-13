"""
forecast_flashlight_video.py — "night walk" MP4 of n34r seasonal forecasts.

PROTOTYPE. A literal-flashlight companion to forecast_flashlight.py, reusing
its beam() (n34r seasonal MMM + Gaussian 1991-2020 hindcast-error spread) and
obs_n34r_seasonal() (scaled ERSSTv5 relative Niño-3.4). Fixed axes, dark
background:

  - Beam = probability density. Each target-season column is the forecast
    Normal(MMM, sd) pdf, tone-mapped to brightness with one fixed scale for
    the whole video, so every column carries the same total light and the
    beam dims as it widens with lead.
  - Walker. The light source sits on the observed path at the last season
    observed when the forecast was issued (center month = init - 2). The
    beam's neck from the source to the first forecast season (lead 1) is a
    visual connector — mean/sd linearly interpolated from (obs, SD_SOURCE) —
    not a forecast. A flashlight glyph sits with its lens on the source,
    rotated (in screen space, so the angle looks right despite the axes'
    month/°C aspect) to point at the lead-1 forecast mean, i.e. along the
    beam's neck.
  - Night walk, lit path. The observed path behind the walker is drawn
    solid. Ahead of the walker the (already observed) path is dark except
    where light falls on it: its opacity is the light intensity at that
    point (current beam + afterglow), (light / LIT_FULL) ** LIT_GAMMA. Where the
    forecast was right the path shows up inside the beam; where it missed,
    the path stays dark — forecast error is shown by the light itself, with
    no separate scoring.
  - Afterglow. Earlier beams leave a glow that decays (e-folding GLOW_TAU
    months) to a dim floor (GLOW_FLOOR of its initial strength) rather than
    to zero, so every earlier forecast stays faintly visible. Earlier beams
    combine by pixelwise max, not sum: with a nonzero floor, summing a year
    of overlapping beams saturates and swamps the current beam.

Motion is continuous at constant speed (SECONDS_PER_ISSUE per issue, no
pauses except FINAL_HOLD_S at the end): beam mean/sd/position are linearly
interpolated from one issue's beam to the next. Those in-between frames are
animation, not forecasts (stated in the subtitle); the "Issued <Mon YYYY>"
box names the nearest real issue and is fully opaque only on it.

Run:
    mamba run -n pangeo-local python scripts/forecast_flashlight_video.py
    mamba run -n pangeo-local python scripts/forecast_flashlight_video.py --stills 0,5,5.5,11 --stills-dir /tmp/x
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FFMpegWriter
from matplotlib.collections import LineCollection
from scipy.stats import norm

import config
from latest_forecast import _resolve_init_idx, _season_label
from forecast_flashlight import beam, obs_n34r_seasonal

FPS = 30
SECONDS_PER_ISSUE = 1.0  # constant-speed transition from one issue to the next
FINAL_HOLD_S = 3.0
NX, NY = 900, 500  # density image resolution
SD_SOURCE = 0.04  # beam half-width (°C) at the walker's hand
DENSITY_SCALE = 0.9  # brightness = 1 - exp(-pdf / DENSITY_SCALE), fixed for all frames
END_FADE_MONTHS = 0.6  # soft edge at the beam's far end
GLOW_TAU = 2.0  # afterglow e-folding time, months
GLOW_GAIN = 0.5
GLOW_FLOOR = 0.15  # afterglow never decays below this fraction of GLOW_GAIN
GLOW_RISE = 0.5  # months for a departing beam to hand over to afterglow
LIT_FULL = 0.8  # light intensity at which the path ahead is fully opaque
LIT_GAMMA = 1.5  # >1 keeps the path dark under dim light (incl. the afterglow floor)
LIT_PATH_RGB = (0.38, 0.84, 1.0)  # path ahead of the walker, where lit
PATH_SAMPLES_PER_MONTH = 12
BG = "#06080e"
LIGHT_RGB = np.array([1.0, 0.90, 0.62])
PATH_COLOR = "0.80"
# Flashlight glyph in screen pixels (before FLASHLIGHT_SCALE), lens at the origin pointing +x.
FLASHLIGHT_SCALE = 1.2
FLASHLIGHT_PARTS = {
    "body": ([(-12, 5), (-44, 5), (-44, -5), (-12, -5)], "#4a5160"),
    "head": ([(0, 9), (-12, 5), (-12, -5), (0, -9)], "#6b7385"),
    "tail": ([(-44, 4), (-48, 4), (-48, -4), (-44, -4)], "#2c313b"),
    "lens": ([(0, 8), (-2.5, 8), (-2.5, -8), (0, -8)], "#fff6d6"),
}
FG = "0.85"


def _mi(ts):
    """Month index (months since year 0) of a timestamp."""
    ts = pd.Timestamp(ts)
    return ts.year * 12 + ts.month - 1


def _mi_to_ts(m):
    return pd.Timestamp(int(m) // 12, int(m) % 12 + 1, 1)


def beam_knots(b, init, obs_mi):
    """Knot arrays (x in month index, mu, sd) for one issue: the walker at
    init-2 on the observed path, then leads 1..10 of the forecast."""
    src = _mi(init) - 2
    x = np.concatenate([[src], [_mi(t) for t in b.time.values]])
    mu = np.concatenate([[obs_mi[src]], b["mean"].values])
    sd = np.concatenate([[SD_SOURCE], b["sd"].values])
    return x.astype(float), mu, sd


def render(x, mu, sd, gx, gy):
    """Brightness image (NY, NX) for a beam given its knots."""
    inside = (gx >= x[0]) & (gx <= x[-1])
    m = np.interp(gx, x, mu)
    s = np.interp(gx, x, sd)
    pdf = norm.pdf(gy[:, None], m[None, :], s[None, :])
    end_fade = np.clip((x[-1] - gx) / END_FADE_MONTHS, 0, 1)
    return (1 - np.exp(-pdf / DENSITY_SCALE)) * (inside * end_fade)[None, :]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--first-init", default="2025-05")
    parser.add_argument("--last-init", default="2026-11")
    parser.add_argument("--stills", default=None, help="Comma-separated issue-time values (0 = first issue; 5.5 = halfway 6th->7th) to write as PNGs instead of the video.")
    parser.add_argument("--stills-dir", default=None)
    args = parser.parse_args()

    ds = config.load_nino34_verification()
    start = ds.indexes["S"].to_datetimeindex(time_unit="ms")
    i0 = _resolve_init_idx(start, args.first_init)
    i1 = _resolve_init_idx(start, args.last_init)
    spec = {"var": "ssta_rel", "prefix": "n34r", "name": "relative Nino 3.4", "scale": config.rel_scaling_factor(ds)}
    inits = [start[i] for i in range(i0, i1 + 1)]
    beams = [beam(ds, start, spec, i) for i in range(i0, i1 + 1)]

    obs = obs_n34r_seasonal()
    obs_mi = {_mi(t): float(v) for t, v in zip(obs.time.values, obs.values)}
    knots = [beam_knots(b, init, obs_mi) for b, init in zip(beams, inits)]
    n = len(knots)

    # Fixed axes (same y rule as forecast_flashlight.py, 0.5 °C steps).
    x_lo = knots[0][0][0] - 3
    x_hi = knots[-1][0][-1]
    obs_win = [v for m, v in obs_mi.items() if x_lo <= m <= x_hi]
    y_lo = min(*obs_win, *(float(b.q05.min()) for b in beams))
    y_hi = max(*obs_win, *(float(b.q95.max()) for b in beams))
    ylim = (np.floor(y_lo * 2) / 2, np.ceil(y_hi * 2) / 2)
    gx = np.linspace(x_lo - 0.5, x_hi + 0.5, NX)
    gy = np.linspace(ylim[0], ylim[1], NY)
    imgs = [render(*k, gx, gy) for k in knots]

    # Dense samples along the full observed path in the window, for the lit path ahead.
    obs_months = np.array(sorted(m for m in obs_mi if x_lo <= m <= x_hi), dtype=float)
    path_x = np.linspace(obs_months[0], obs_months[-1], int((obs_months[-1] - obs_months[0]) * PATH_SAMPLES_PER_MONTH) + 1)
    path_y = np.interp(path_x, obs_months, [obs_mi[m] for m in obs_months])
    path_ix = np.clip(np.searchsorted(gx, path_x), 0, NX - 1)
    path_iy = np.clip(np.searchsorted(gy, path_y), 0, NY - 1)
    path_segs = np.stack([np.column_stack([path_x[:-1], path_y[:-1]]), np.column_stack([path_x[1:], path_y[1:]])], axis=1)

    def state(tau):
        k = min(int(np.floor(tau)), n - 1)
        f = tau - k
        if k == n - 1 or f == 0:
            x, mu, sd = knots[k]
            cur = imgs[k]
        else:
            x, mu, sd = (a * (1 - f) + c * f for a, c in zip(knots[k], knots[k + 1]))
            cur = render(x, mu, sd, gx, gy)
        glow = np.zeros_like(cur)
        for j in range(n):
            age = tau - j
            if age > 0:
                decay = GLOW_FLOOR + (1 - GLOW_FLOOR) * np.exp(-age / GLOW_TAU)
                glow = np.maximum(glow, GLOW_GAIN * min(1.0, age / GLOW_RISE) * decay * imgs[j])
        return (x[0], mu[0]), (x[1], mu[1]), np.clip(cur + glow, 0, 1)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100, facecolor=BG)
    ax = fig.add_axes([0.07, 0.12, 0.90, 0.78], facecolor=BG)
    im = ax.imshow(np.zeros((NY, NX, 3)), extent=(gx[0], gx[-1], gy[0], gy[-1]), origin="lower", aspect="auto", interpolation="bilinear", zorder=1)
    ax.axhline(0, color="0.35", lw=0.8, zorder=2)
    (path_line,) = ax.plot([], [], "-", color=PATH_COLOR, lw=2.5, zorder=3)
    (path_dots,) = ax.plot([], [], "o", color=PATH_COLOR, ms=4, zorder=3)
    walker_halo = [ax.plot([], [], "o", color="#fff3c4", ms=ms, alpha=a, zorder=5)[0] for ms, a in [(28, 0.08), (18, 0.18), (10, 0.5)]]
    flashlight = {
        name: ax.add_patch(Polygon(np.array(xy, dtype=float), closed=True, facecolor=color, edgecolor="0.75", lw=0.6, zorder=6))
        for name, (xy, color) in FLASHLIGHT_PARTS.items()
    }

    def place_flashlight(src, aim):
        to_px = ax.transData
        p0, p1 = to_px.transform(src), to_px.transform(aim)
        theta = np.arctan2(p1[1] - p0[1], p1[0] - p0[0])
        rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        for name, (xy, _) in FLASHLIGHT_PARTS.items():
            flashlight[name].set_xy(to_px.inverted().transform(FLASHLIGHT_SCALE * np.asarray(xy, dtype=float) @ rot.T + p0))
    lit_path = LineCollection(path_segs, linewidths=2.5, zorder=4, capstyle="butt")
    ax.add_collection(lit_path)

    ticks = np.arange(np.ceil(x_lo), x_hi + 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{_season_label(_mi_to_ts(t).month)}\n{_mi_to_ts(t).year}" if _mi_to_ts(t).month in (1, 7) else _season_label(_mi_to_ts(t).month) for t in ticks], color=FG, fontsize=9)
    ax.set_xlim(gx[0], gx[-1])
    ax.set_ylim(ylim)
    ax.set_yticks(np.arange(ylim[0], ylim[1] + 0.25, 0.5))
    ax.tick_params(colors=FG)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("0.4")
    ax.set_ylabel("n34r (°C)", color=FG)
    fig.text(0.07, 0.962, "NMME relative Niño-3.4 seasonal forecasts — a night walk", color="white", fontsize=15, fontweight="bold")
    fig.text(0.07, 0.925, "Beam brightness ∝ forecast probability density (MMM ± 1991-2020 hindcast error). Observed path (ERSSTv5) ahead is visible only where lit.\nMotion between monthly issues is interpolated.", color=FG, fontsize=9.5, va="center", linespacing=1.4)
    issue_box = ax.text(0.015, 0.97, "", transform=ax.transAxes, ha="left", va="top", fontsize=16, fontweight="bold", color="white",
                        bbox=dict(boxstyle="round", facecolor="#1b2030", edgecolor="0.5"))
    fig.text(0.07, 0.03, "━ walked path", color=PATH_COLOR, fontsize=10, va="bottom")
    fig.text(0.20, 0.03, "━ observed path ahead, where the light reaches it", color=LIT_PATH_RGB, fontsize=10, va="bottom")

    def draw(tau):
        (src_x, src_y), aim, light = state(tau)
        im.set_data(np.asarray(matplotlib.colors.to_rgb(BG))[None, None, :] * (1 - light[..., None]) + LIGHT_RGB * light[..., None])

        months = sorted(m for m in obs_mi if x_lo <= m <= src_x)
        px = months + ([src_x] if src_x > months[-1] else [])
        py = [obs_mi[m] for m in months] + ([src_y] if src_x > months[-1] else [])
        path_line.set_data(px, py)
        path_dots.set_data(months, [obs_mi[m] for m in months])
        for h in walker_halo:
            h.set_data([src_x], [src_y])
        place_flashlight((src_x, src_y), aim)

        alpha = np.clip(light[path_iy, path_ix] / LIT_FULL, 0, 1) ** LIT_GAMMA * (path_x > src_x)
        rgba = np.zeros((len(path_segs), 4))
        rgba[:, :3] = LIT_PATH_RGB
        rgba[:, 3] = alpha[1:]
        lit_path.set_color(rgba)

        nearest = min(int(round(tau)), n - 1)
        issue_box.set_text(f"Issued {inits[nearest]:%b %Y}")
        issue_box.set_alpha(1 - 0.9 * abs(tau - nearest))

    if args.stills:
        outdir = Path(args.stills_dir or config.PLOTS_DIR_FORECAST_FLASHLIGHT / "stills")
        outdir.mkdir(parents=True, exist_ok=True)
        for tau in (float(t) for t in args.stills.split(",")):
            draw(tau)
            out = outdir / f"still_{tau:05.2f}.png"
            fig.savefig(out, facecolor=BG)
            print(f"  wrote {out}")
        return

    taus = list(np.linspace(0, n - 1, int((n - 1) * SECONDS_PER_ISSUE * FPS) + 1))
    taus += [float(n - 1)] * int(FINAL_HOLD_S * FPS)

    config.PLOTS_DIR_FORECAST_FLASHLIGHT.mkdir(parents=True, exist_ok=True)
    out = config.PLOTS_DIR_FORECAST_FLASHLIGHT / f"n34r_seasonal_{inits[0]:%Y-%m}_{inits[-1]:%Y-%m}_nightwalk.mp4"
    writer = FFMpegWriter(fps=FPS, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "20"])
    with writer.saving(fig, out, dpi=100):
        for tau in taus:
            draw(tau)
            writer.grab_frame(facecolor=BG)
    plt.close(fig)
    print(f"  wrote {out} ({len(taus)} frames, {len(taus) / FPS:.1f} s)")


if __name__ == "__main__":
    main()
