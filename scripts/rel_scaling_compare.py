"""
rel_scaling_compare.py — evidence for the relative Niño-3.4 scaling factor's
ensemble-member-pooling choice.

config.rel_scaling_factor scales the model relative-index (ssta_rel)
variance to the observed Niño-3.4 variance using a denominator that pools
members via the grand mean, sqrt(var(['S', 'M'])) — the pooled standard
deviation over the flattened (start, member) sample (C, chosen production
formula as of 2026-07-07; see config.rel_scaling_factor's docstring for the
switch from the original per-member choice). This script provides evidence
on all three pairwise comparison axes among the candidate denominators:

  - A = per-member pooling, sqrt(var('S').mean('M')) — the original choice,
    kept here (not in config.py) as the alternative that motivated the
    switch.
  - B = ensemble-mean pooling, std('S') of ssta_rel.mean('M') — flawed: it
    shrinks variance wherever skill is low, mixing predictability into a
    quantity meant to capture only variance.
  - C = grand-mean pooling, sqrt(var(['S', 'M'])) — config.rel_scaling_factor,
    the chosen production formula.

A and C are equal *in expectation* but differ in finite samples: C
additionally carries the between-member spread of the per-member time-means
(law of total variance), whereas A averages out that term. Unlike B, neither
A nor C is flawed by construction — both are valid member-based estimators;
A vs C quantifies how much the finite-sample difference actually matters in
practice (it motivated the 2026-07-07 switch: MSESS was not worse — slightly
better on average — for C, which is also the simpler estimator to describe).

For all three axes the same verification target is used — an observed
relative-index scaled to match the observed Niño-3.4 variance (the paper's
own RONI).

  - Anomaly correlation is expected to be IDENTICAL across A, B, and C
    (correlation is invariant to a per-(model, month, L) scalar factor
    within a start-month group) — printed as a sanity check, not plotted.
  - MSESS is NOT scale-invariant and differs across the three. Plotted as
    (month x lead) heatmaps, by both start month and target month, matching
    skill.py's own figures, plus difference grids for all three pairs
    (A - B, A - C, B - C).

Exploratory/diagnostic only: this does not feed production figures.

Run:
    mamba run -n pangeo-local python scripts/rel_scaling_compare.py
"""
import numpy as np

import matplotlib
matplotlib.use("Agg")

import config
from skill import ac_by_start, msess_by_start, add_mmm, to_target_month, plot_skill_grid


def _factor_permember(ds, period):
    """Alternative (previously chosen, still valid) denominator: per-member
    temporal variance, averaged across members, then square-rooted,
    sqrt(var('S').mean('M')) — rather than config.rel_scaling_factor's
    grand-mean sqrt(var(['S', 'M'])). Equal in expectation to the grand-mean
    form; differs only in finite samples (see config.rel_scaling_factor's
    docstring). Kept out of config.py: it was the production choice before
    2026-07-07 and is retained here purely as the comparison baseline that
    motivated the switch."""
    start_year, end_year = period
    year = ds.S.dt.year
    ref = ds.where((year >= start_year) & (year <= end_year), drop=True)

    num = ref.obsa.groupby("S.month").std("S")
    den = np.sqrt(ref.ssta_rel.groupby("S.month").var("S").mean("M"))
    return num / den


def _factor_ensmean(ds, period):
    """Alternative (flawed) denominator: std('S') of the ensemble-mean
    ssta_rel, rather than a member-pooled estimator. Kept out of config.py
    deliberately — this exists only to demonstrate why the member-pooled
    estimators (A, C) are preferred (see config.rel_scaling_factor's
    docstring)."""
    start_year, end_year = period
    year = ds.S.dt.year
    ref = ds.where((year >= start_year) & (year <= end_year), drop=True)

    num = ref.obsa.groupby("S.month").std("S")
    den = ref.ssta_rel.mean("M").groupby("S.month").std("S")
    return num / den


def _scaled_forecast(ds_p, factor):
    """Ensemble-mean scaled relative-index forecast, dims (model, S, L),
    factor keyed by each start's calendar month (pointwise/vectorized sel)."""
    f = factor.sel(month=ds_p.S.dt.month)
    return (ds_p.ssta_rel * f).mean("M")


def main():
    config.PLOTS_DIR_REL_SCALING_COMPARE.mkdir(parents=True, exist_ok=True)

    ds = config.load_nino34_verification()
    period = (config.CLIM_START_YEAR, config.CLIM_END_YEAR)
    period_sel = f"{period[0]}-01-01", f"{period[1]}-12-01"
    ds_p = ds.sel(S=slice(*period_sel))
    print(
        f"  rel_scaling_compare: restricted to S={period_sel[0]}..{period_sel[1]}, "
        f"shape={dict(ds_p.ssta_rel.sizes)}"
    )

    factor_A = _factor_permember(ds, period=period)              # per-member pooling (alternative)
    factor_B = _factor_ensmean(ds, period=period)                # ensemble-mean pooling (flawed)
    factor_C, _ = config.rel_scaling_factor(ds, period=period)   # grand-mean pooling (chosen, production)

    fcst_A = _scaled_forecast(ds_p, factor_A)
    fcst_B = _scaled_forecast(ds_p, factor_B)
    fcst_C = _scaled_forecast(ds_p, factor_C)

    # Verification target: RONI — obsa_rel scaled to match obs Niño-3.4
    # variance (the paper's own obs-only ratio). Identical for A, B, and C.
    obs_num = ds_p.obsa.groupby("S.month").std("S")
    obs_den = ds_p.obsa_rel.groupby("S.month").std("S")
    obs_scale = obs_num / obs_den
    obs_target = ds_p.obsa_rel * obs_scale.sel(month=ds_p.S.dt.month)

    ac_A = ac_by_start(fcst_A, obs_target)
    ac_B = ac_by_start(fcst_B, obs_target)
    ac_C = ac_by_start(fcst_C, obs_target)
    max_ac_diff_ab = float(np.abs(ac_A - ac_B).max())
    max_ac_diff_ac = float(np.abs(ac_A - ac_C).max())
    max_ac_diff_bc = float(np.abs(ac_B - ac_C).max())
    print(
        f"  rel_scaling_compare: max |AC_A - AC_B| = {max_ac_diff_ab:.2e}, "
        f"max |AC_A - AC_C| = {max_ac_diff_ac:.2e}, "
        f"max |AC_B - AC_C| = {max_ac_diff_bc:.2e} "
        "(all expected ~0 — correlation is invariant to the scalar factor)"
    )

    msess_A = msess_by_start(fcst_A, obs_target)
    msess_A = add_mmm(msess_A, fcst_A, obs_target, msess_by_start)
    msess_B = msess_by_start(fcst_B, obs_target)
    msess_B = add_mmm(msess_B, fcst_B, obs_target, msess_by_start)
    msess_C = msess_by_start(fcst_C, obs_target)
    msess_C = add_mmm(msess_C, fcst_C, obs_target, msess_by_start)

    msess_A_target = to_target_month(msess_A)
    msess_B_target = to_target_month(msess_B)
    msess_C_target = to_target_month(msess_C)
    msess_diff_ab = msess_A - msess_B
    msess_diff_ab_target = msess_A_target - msess_B_target
    msess_diff_ac = msess_A - msess_C
    msess_diff_ac_target = msess_A_target - msess_C_target
    msess_diff_bc = msess_B - msess_C
    msess_diff_bc_target = msess_B_target - msess_C_target

    print(f"  rel_scaling_compare: mean MSESS(A) - MSESS(B) = {float(msess_diff_ab.mean()):.4f}")
    print(f"  rel_scaling_compare: mean MSESS(A) - MSESS(C) = {float(msess_diff_ac.mean()):.4f}")
    print(f"  rel_scaling_compare: mean MSESS(B) - MSESS(C) = {float(msess_diff_bc.mean()):.4f}")

    period_label = f"{period[0]}-{period[1]}"

    plot_skill_grid(
        msess_A, f"n34r MSESS, per-member pooled denominator (A, alternative), by start month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_A_start.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_A_target, f"n34r MSESS, per-member pooled denominator (A, alternative), by target month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_A_target.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_B, f"n34r MSESS, ensemble-mean denominator (B, flawed), by start month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_B_start.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_B_target, f"n34r MSESS, ensemble-mean denominator (B, flawed), by target month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_B_target.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_C, f"n34r MSESS, grand-mean pooled denominator (C, chosen), by start month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_C_start.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_C_target, f"n34r MSESS, grand-mean pooled denominator (C, chosen), by target month ({period_label})",
        "MSESS", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_C_target.png", vmax=1.0,
    )
    plot_skill_grid(
        msess_diff_ab, f"n34r MSESS, per-member scaling - ens. mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_AB_start.png",
    )
    plot_skill_grid(
        msess_diff_ab_target, f"n34r MSESS, per-member scaling - ens. mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_AB_target.png",
    )
    # A-C diff color range: SPEAR and GEOSS2S have large negative outliers
    # (min ~ -0.34, -0.40) that would otherwise saturate a full-range
    # diverging colormap and hide the much smaller (~+-0.03) signal in the
    # other five models/MMM. Clip to +-0.03 so those panels are legible;
    # SPEAR/GEOSS2S render fully saturated (clipped) rather than showing
    # their true extreme values — acceptable since their magnitude is
    # already reported in console output (mean MSESS(A) - MSESS(C) above)
    # and in the unclipped msess_A_*/msess_C_* panels. cmap must be passed
    # explicitly alongside vmin/vmax — see plot_skill_grid's docstring.
    ac_diff_vlim = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")
    plot_skill_grid(
        msess_diff_ac, f"n34r MSESS, per-member scaling - grand-mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_AC_start.png", **ac_diff_vlim,
    )
    plot_skill_grid(
        msess_diff_ac_target, f"n34r MSESS, per-member scaling - grand-mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_AC_target.png", **ac_diff_vlim,
    )
    plot_skill_grid(
        msess_diff_bc, f"n34r MSESS, ens. mean scaling - grand-mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_BC_start.png",
    )
    plot_skill_grid(
        msess_diff_bc_target, f"n34r MSESS, ens. mean scaling - grand-mean scaling ({period_label})",
        "MSESS diff", config.PLOTS_DIR_REL_SCALING_COMPARE / "n34r_scaling_msess_diff_BC_target.png",
    )


if __name__ == "__main__":
    main()
