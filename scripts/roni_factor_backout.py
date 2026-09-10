"""
roni_factor_backout.py — back out CPC's RONI scaling factor from ERSSTv6.

CPC's Relative Oceanic Niño Index (RONI) page
(https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/roni/)
describes the index as:

    "3-month running mean of ERSST.v6 SST anomalies in the Niño 3.4 region
    (5°N-5°S, 120°-170°W) with average tropical mean (20°N-20°S) SST
    anomalies subtracted [both anomalies over a 1991-2020 base period]. The
    difference is then adjusted so the variance equals the original Niño
    3.4 index."

That fixes the anomaly recipe (box definitions, 1991-2020 base period,
3-month running mean) but not the scaling factor's numeric value, which
varies with the calendar month.

The current answer, and the order the sections below are built in (see
docs/CPC_RONI_scaling.md for the write-up):

 1. check_monthly_seasonal_consistency() — CPC's published monthly
    (Rnino34.ascii.txt) and seasonal (RONI.ascii.txt) products are mutually
    consistent, using no ERSSTv6 data at all. The factor is therefore
    applied to *monthly* values, with the seasonal product nothing more
    than their 3-month running mean.
 2. backout_monthly_factors_direct() — the primary estimate: each calendar
    month's factor by OLS through the origin against the published monthly
    series.
 3. rounding_convention_check() — the ~0.005 residual left by (2) is the
    publication convention: CPC's published values are floored to two
    decimals, not rounded to nearest. Refitting with a free intercept
    absorbs it and drops the residual to the quantization floor.
 4. reconstruct_from_factors() / reconstruction_tables() — the full round
    trip, ERSSTv6 to published values, monthly and seasonal.

backout_factors(), backout_factors_recent() and backout_monthly_factors()
are superseded estimators kept for comparison: the first two fit one factor
per *season* (an approximation, since the factor varies within a season),
the third recovers monthly factors from the seasonal product alone. They
quantify what the simpler models cost; they are not the recipe.

Exploratory/diagnostic only — not part of the production pipeline, not
referenced by any other script. Writes its printed report and the
underlying tables (CSV) to plots/roni_factor_backout/.

Run:
    mamba run -n pangeo-local python scripts/roni_factor_backout.py
"""
import numpy as np
import pandas as pd
import xarray as xr

import config

# Niño-3.4 box per CPC's RONI page (5N-5S, 120W-170W) — same box as
# config.N34_LAT / config.N34_LON (190-240E == 170W-120W); reused rather
# than restated to guarantee the two projects' Niño-3.4 boxes match.
N34_LAT, N34_LON = config.N34_LAT, config.N34_LON
# Tropical-mean box per the RONI page (20N-20S, all lon) — same as
# config.TROPICS_LAT.
TROPICS_LAT = config.TROPICS_LAT
BASE_START, BASE_END = config.CLIM_START_YEAR, config.CLIM_END_YEAR  # 1991, 2020

MONTH_LETTER = "JFMAMJJASOND"


def _season_label(center_month):
    """3-letter CPC season code for a center month (1-12), e.g. 1 -> 'DJF'."""
    idx = center_month - 1  # 0-based
    letters = [MONTH_LETTER[(idx + off) % 12] for off in (-1, 0, 1)]
    return "".join(letters)


def _box_mean(da, lat, lon=None, weighted=True):
    """Cosine-latitude-weighted box average by default. weighted=False gives
    a plain unweighted grid-cell mean — scientifically wrong (a 2-degree
    grid box near the pole covers far less area than one at the equator),
    kept here only to test whether CPC's own RONI product skips area
    weighting (per the user, 2026-09-10: "this is scientifically WRONG but
    we are trying to replicate").
    """
    da = da.sortby("lat")
    sel = da.sel(lat=lat) if lon is None else da.sel(lat=lat, lon=lon)
    if not weighted:
        return sel.mean(["lat", "lon"])
    weights = np.cos(np.deg2rad(da.lat))
    return sel.weighted(weights).mean(["lat", "lon"])


def _ersstv6_monthly_anomalies(nc=None, weighted=True):
    """Monthly ERSSTv6 Niño-3.4 and tropical-mean anomalies (1991-2020 base),
    no running mean. Returns xr.Dataset (dims: time) with n34_anom,
    trop_anom.
    """
    nc = config.ERSSTV6_NC if nc is None else nc
    ds = xr.open_dataset(nc)
    print(f"_ersstv6_monthly_anomalies: {nc}, "
          f"{ds.sst.time.values[0]} .. {ds.sst.time.values[-1]}, "
          f"{ds.sst.shape}, weighted={weighted}")

    n34 = _box_mean(ds.sst, N34_LAT, N34_LON, weighted=weighted)
    trop = _box_mean(ds.sst, TROPICS_LAT, weighted=weighted)

    clim_sel = dict(time=slice(f"{BASE_START}-01-01", f"{BASE_END}-12-31"))
    n34_clim = n34.sel(**clim_sel).groupby("time.month").mean("time")
    trop_clim = trop.sel(**clim_sel).groupby("time.month").mean("time")

    n34_anom = n34.groupby("time.month") - n34_clim
    trop_anom = trop.groupby("time.month") - trop_clim
    return xr.Dataset({"n34_anom": n34_anom, "trop_anom": trop_anom})


def load_ersstv6_3mo_anomalies(nc=None, weighted=True):
    """Monthly ERSSTv6 Niño-3.4 and tropical-mean anomalies (1991-2020 base),
    3-month running mean (centered), with CPC season/year labels.

    Returns
    -------
    pandas.DataFrame
        Columns: year, month (center month, 1-12), SEAS, YR, n34_3mo,
        trop_3mo, diff_3mo (= n34_3mo - trop_3mo, the unscaled relative
        index).
    """
    anom = _ersstv6_monthly_anomalies(nc, weighted=weighted)
    n34_anom = anom.n34_anom.rolling(time=3, center=True).mean()
    trop_anom = anom.trop_anom.rolling(time=3, center=True).mean()

    df = pd.DataFrame({
        "year": n34_anom.time.dt.year.values,
        "month": n34_anom.time.dt.month.values,
        "n34_3mo": n34_anom.values,
        "trop_3mo": trop_anom.values,
    }).dropna()
    df["diff_3mo"] = df.n34_3mo - df.trop_3mo
    df["SEAS"] = df.month.map(_season_label)
    # CPC labels each season by its center month's own calendar year (e.g.
    # NDJ 1997 = Nov97,Dec97,Jan98, center month Dec97 -> "1997"; DJF 1998 =
    # Dec97,Jan98,Feb98, center month Jan98 -> "1998") — verified against
    # RONI.ascii.txt's 1997-98 El Nino rows.
    df["YR"] = df.year
    return df


def monthly_diff_lookup(nc=None, weighted=True):
    """{(year, month): diff} for every ERSSTv6 calendar month, where
    diff = n34_anom - trop_anom (monthly, no running mean) — the per-month
    building block for backout_monthly_factors().
    """
    anom = _ersstv6_monthly_anomalies(nc, weighted=weighted)
    diff = anom.n34_anom - anom.trop_anom
    years = diff.time.dt.year.values
    months = diff.time.dt.month.values
    values = diff.values
    return {(int(y), int(m)): v for y, m, v in zip(years, months, values)}


def load_roni_table(path=None):
    path = config.RONI_TXT if path is None else path
    df = pd.read_csv(path, sep=r"\s+")
    print(f"load_roni_table: {path}, {len(df)} rows, "
          f"{df.SEAS.iloc[0]} {df.YR.iloc[0]} .. {df.SEAS.iloc[-1]} {df.YR.iloc[-1]}")
    return df


def load_rnino34_table(path=None):
    """CPC's published *monthly* (not 3-month-running-mean) relative
    Nino-3.4 series — Rnino34.ascii.txt, a different CPC product from the
    seasonal RONI.ascii.txt (config.RONI_TXT). Lets each calendar month's
    factor be fit directly against its own monthly target, rather than only
    indirectly through the seasonal aggregates (see
    backout_monthly_factors_direct()).
    """
    path = config.RNINO34_TXT if path is None else path
    df = pd.read_csv(path, sep=r"\s+")
    df = df.rename(columns={"MTH": "month", "YR": "year"})
    print(f"load_rnino34_table: {path}, {len(df)} rows, "
          f"{df.year.iloc[0]}-{df.month.iloc[0]:02d} .. {df.year.iloc[-1]}-{df.month.iloc[-1]:02d}")
    return df


def check_monthly_seasonal_consistency():
    """Self-consistency check between CPC's two published products, with no
    involvement of ERSSTv6 or any fitted factor: does a 3-month centered
    running mean of the published *monthly* series (Rnino34.ascii.txt)
    reproduce the published *seasonal* series (RONI.ascii.txt)?

    This is the first and most basic check, logically prior to any factor
    fitting: it establishes that CPC's scaling factor is applied at the
    monthly level (with the seasonal series nothing more than a 3-month
    average of already-scaled monthly values), which is the premise the
    rest of this module's monthly-factor fitting relies on. See
    docs/CPC_RONI_scaling.md Section 2.

    Returns
    -------
    pandas.DataFrame
        Matched rows: SEAS, YR, ANOM_seas (published seasonal), 3mo_avg
        (rolling mean of published monthly), diff.
    """
    rnino34 = load_rnino34_table().sort_values(["year", "month"]).reset_index(drop=True)
    rnino34["t"] = pd.to_datetime(dict(year=rnino34.year, month=rnino34.month, day=1))
    rnino34 = rnino34.set_index("t")
    n_gaps = len(pd.date_range(rnino34.index.min(), rnino34.index.max(), freq="MS")) - len(rnino34)
    assert n_gaps == 0, f"Rnino34.ascii.txt has {n_gaps} gap(s) in monthly coverage"

    rnino34["3mo_avg"] = rnino34.ANOM.rolling(3, center=True).mean()
    rnino34 = rnino34.dropna(subset=["3mo_avg"])
    rnino34["SEAS"] = rnino34.month.map(_season_label)
    rnino34["YR"] = rnino34.year

    roni = load_roni_table()
    merged = rnino34.reset_index().merge(
        roni, on=["SEAS", "YR"], how="inner", suffixes=("_mo", "_seas"))
    merged = merged.rename(columns={"ANOM_seas": "ANOM_seas", "ANOM_mo": "ANOM_monthly"})
    merged["diff"] = merged.ANOM_seas - merged["3mo_avg"]

    print(f"check_monthly_seasonal_consistency: matched {len(merged)}/{len(roni)}, "
          f"max|diff|={merged['diff'].abs().max():.4f}, "
          f"rms={np.sqrt((merged['diff']**2).mean()):.4f}, "
          f"mean={merged['diff'].mean():.6f}")
    return merged[["SEAS", "YR", "ANOM_seas", "3mo_avg", "diff"]]


def load_merged(weighted=True):
    """ERSSTv6-derived diff_3mo joined to published RONI, on (SEAS, YR)."""
    ersst = load_ersstv6_3mo_anomalies(weighted=weighted)
    roni = load_roni_table()
    merged = ersst.merge(roni, on=["SEAS", "YR"], how="inner", validate="one_to_one")
    print(f"load_merged: matched {len(merged)}/{len(roni)} published RONI rows "
          f"to ERSSTv6-derived seasons")
    return merged


def backout_factors(merged=None):
    """Back out each season's RONI factor from a single (year, season) row,
    then validate it against every *other* year of that same season.

    RONI = factor(season) * diff_3mo is a per-season scalar multiply, not a
    statistical fit — so factor(season) should be recoverable exactly (up to
    RONI's 2-decimal-place rounding) from one row, and that single-row
    factor should reproduce every other year's published value. This is a
    direct test of "each year has the same factor per season", rather than
    a regression that would silently paper over a factor that actually
    drifts across years.

    The reference row per season is the year with the largest |diff_3mo| —
    rounding noise in a 2-decimal published value is a fixed absolute error,
    so dividing by the largest-amplitude diff minimizes its relative effect
    on the recovered factor.
    """
    if merged is None:
        merged = load_merged()

    rows = []
    season_order = [_season_label(m) for m in range(1, 13)]
    for seas in season_order:
        g = merged[merged.SEAS == seas].sort_values("YR").reset_index(drop=True)
        x, y = g.diff_3mo.values, g.ANOM.values

        per_year_factor = y / x  # one ratio per (season, year) row

        ref = np.argmax(np.abs(x))
        factor_single_row = y[ref] / x[ref]
        ref_year = int(g.YR.values[ref])

        other = np.arange(len(g)) != ref
        resid_other = y[other] - factor_single_row * x[other]

        rows.append(dict(
            SEAS=seas, n=len(g), ref_year=ref_year,
            factor_single_row=factor_single_row,
            factor_mean_all_years=per_year_factor.mean(),
            factor_std_all_years=per_year_factor.std(ddof=0),
            max_abs_err_other_years=np.max(np.abs(resid_other)),
            rms_err_other_years=np.sqrt(np.mean(resid_other ** 2)),
        ))

    return pd.DataFrame(rows)


def backout_factors_recent(n_recent=20, merged=None):
    """Estimate each season's factor (and its accuracy) from the n_recent
    most recent years of that season, via ordinary least-squares regression
    through the origin: factor = sum(x*y) / sum(x**2), x = diff_3mo,
    y = published RONI.

    Regression through the origin naturally weights each row by x_i**2, so
    near-neutral years (diff_3mo ~ 0) — which make the simple ratio
    RONI/diff blow up, see backout_factors()'s factor_mean_all_years
    diagnostic — contribute almost nothing to the fit; strong-ENSO years
    dominate. Restricting to recent years rather than the full 1950-2026
    record checks whether the factor has drifted/been revised, and keeps
    the accuracy estimate representative of the record CPC is minting new
    RONI values against today.

    Accuracy: standard error of the slope, SE(factor) = sqrt(sigma2 /
    sum(x**2)) with sigma2 = sum(residuals**2) / (n - 1) (one estimated
    parameter, slope only), plus the resulting 95% CI (factor +/- 1.96*SE).
    """
    if merged is None:
        merged = load_merged()

    rows = []
    season_order = [_season_label(m) for m in range(1, 13)]
    for seas in season_order:
        g = merged[merged.SEAS == seas].sort_values("YR").tail(n_recent)
        x, y = g.diff_3mo.values, g.ANOM.values
        n = len(x)

        factor = np.sum(x * y) / np.sum(x * x)
        resid = y - factor * x
        sigma2 = np.sum(resid ** 2) / (n - 1)
        se = np.sqrt(sigma2 / np.sum(x * x))

        rows.append(dict(
            SEAS=seas, n=n, years=f"{int(g.YR.min())}-{int(g.YR.max())}",
            factor=factor, se=se,
            ci95_lo=factor - 1.96 * se, ci95_hi=factor + 1.96 * se,
            rms_resid=np.sqrt(np.mean(resid ** 2)),
        ))

    return pd.DataFrame(rows)


def _neighbor_month(year, month, offset):
    m = month + offset
    if m < 1:
        return year - 1, m + 12
    if m > 12:
        return year + 1, m - 12
    return year, m


def backout_monthly_factors(n_recent=20, merged=None, weighted=True):
    """Solve for the 12 calendar-month RONI factors jointly (per the user's
    2026-09-10 correction): CPC computes the *monthly* relative Niño-3.4
    value, factor(month) * diff_monthly(month), first, and only then
    averages 3 consecutive months into a season — not the reverse order
    backout_factors()/backout_factors_recent() assumed (average diff over
    3 months, then apply one per-season factor). The two are equivalent
    only if the factor is constant across a season's 3 months, which the
    fitted per-season values are not (e.g. FMA 1.3222 vs MAM 1.3598, from
    backout_factors_recent()) — so a joint monthly fit is needed to recover
    the true per-month factors.

    Still a linear least-squares problem, just with a different design
    matrix: instead of 12 independent single-column fits (one per season,
    using that season's 3-month-averaged diff as the sole regressor), this
    is one combined fit over all 12*n_recent rows at once, where row i
    (published RONI for season s, year y, center calendar month m) has
    three nonzero entries — at columns m-1, m, m+1 (mod 12, with year
    rollover) — each holding (1/3) * diff_monthly at that neighboring
    (year, month). Solved via numpy.linalg.lstsq; standard errors from the
    OLS covariance sigma2 * (X'X)^-1, sigma2 = SSR / (n_obs - 12).
    """
    if merged is None:
        merged = load_merged(weighted=weighted)
    diffs = monthly_diff_lookup(weighted=weighted)

    cutoff_year = merged.YR.max() - n_recent + 1
    recent = merged[merged.YR >= cutoff_year].sort_values(["SEAS", "YR"]).reset_index(drop=True)

    n_obs = len(recent)
    X = np.zeros((n_obs, 12))
    y = recent.ANOM.values.copy()
    for i, row in enumerate(recent.itertuples()):
        for off in (-1, 0, 1):
            ny, nm = _neighbor_month(row.year, row.month, off)
            X[i, nm - 1] += diffs[(ny, nm)] / 3.0

    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n_obs - 12
    sigma2 = np.sum(resid ** 2) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    out = pd.DataFrame({
        "month": month_names,
        "factor": beta,
        "se": se,
        "ci95_lo": beta - 1.96 * se,
        "ci95_hi": beta + 1.96 * se,
    })
    rms_resid = np.sqrt(np.mean(resid ** 2))
    return out, rms_resid, n_obs


MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def monthly_fit_frame(weighted=True, start_year=1950):
    """CPC's published *monthly* relative Niño-3.4 series joined to the
    ERSSTv6 monthly diff on (year, month) — the (x, y) table every
    monthly-resolution fit and check below is built from.
    """
    diffs = monthly_diff_lookup(weighted=weighted)
    rnino34 = load_rnino34_table().copy()
    rnino34["diff"] = [diffs.get((int(y), int(m))) for y, m in
                       zip(rnino34.year, rnino34.month)]
    rnino34 = rnino34.dropna(subset=["diff"])
    return rnino34[rnino34.year >= start_year].reset_index(drop=True)


def _recent_years(frame, n_recent):
    """The n_recent most recent years of each calendar month, stacked."""
    return pd.concat([frame[frame.month == m].sort_values("year").tail(n_recent)
                      for m in range(1, 13)])


def backout_monthly_factors_direct(n_recent=20, weighted=True):
    """Fit each calendar month's factor directly against CPC's published
    *monthly* relative Nino-3.4 series (Rnino34.ascii.txt), rather than only
    indirectly through the seasonal RONI aggregates (backout_monthly_
    factors()). Each calendar month's factor is now an independent
    OLS-through-origin fit — no shared design matrix needed, since every
    observation maps to exactly one column (no season-driven 3-month
    coupling) — over the n_recent most recent years of that calendar month.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        Per-month factor table (month, factor, se, ci95_lo, ci95_hi,
        rms_resid, n), and the merged (year, month, diff, ANOM) table used
        to fit it.
    """
    rnino34 = monthly_fit_frame(weighted=weighted)

    month_names = MONTH_NAMES
    rows = []
    for m in range(1, 13):
        g = rnino34[rnino34.month == m].sort_values("year").tail(n_recent)
        x, y = g["diff"].values, g.ANOM.values
        n = len(x)

        factor = np.sum(x * y) / np.sum(x * x)
        resid = y - factor * x
        sigma2 = np.sum(resid ** 2) / (n - 1)
        se = np.sqrt(sigma2 / np.sum(x * x))

        rows.append(dict(
            month=month_names[m - 1], n=n,
            years=f"{int(g.year.min())}-{int(g.year.max())}",
            factor=factor, se=se,
            ci95_lo=factor - 1.96 * se, ci95_hi=factor + 1.96 * se,
            rms_resid=np.sqrt(np.mean(resid ** 2)),
        ))

    return pd.DataFrame(rows), rnino34


def validate_direct_monthly_against_seasonal(monthly_direct, merged, n_recent=20,
                                             weighted=True):
    """Apply the directly-fit monthly factors (backout_monthly_factors_
    direct()) to monthly diff, 3-month-average the *scaled* values, and
    compare against the published seasonal RONI (merged.ANOM) — the full
    round-trip check of "fixed monthly factors applied to monthly data,
    then 3-month averaged" against an independent CPC product.
    """
    beta = dict(zip(range(1, 13), monthly_direct.factor.values))
    diffs = monthly_diff_lookup(weighted=weighted)

    cutoff_year = merged.YR.max() - n_recent + 1
    recent = merged[merged.YR >= cutoff_year].sort_values(["SEAS", "YR"]).reset_index(drop=True)

    pred = []
    for row in recent.itertuples():
        p = sum(beta[_neighbor_month(row.year, row.month, off)[1]]
                * diffs[_neighbor_month(row.year, row.month, off)] / 3.0
                for off in (-1, 0, 1))
        pred.append(p)
    resid = recent.ANOM.values - np.array(pred)

    rows = []
    season_order = [_season_label(m) for m in range(1, 13)]
    for seas in season_order:
        mask = recent.SEAS.values == seas
        rows.append(dict(
            SEAS=seas, n=int(mask.sum()),
            rms_resid=np.sqrt(np.mean(resid[mask] ** 2)),
        ))
    out = pd.DataFrame(rows)
    overall_rms = np.sqrt(np.mean(resid ** 2))
    return out, overall_rms


def _fit_slopes_with_offset(g):
    """OLS of the published monthly value on a per-calendar-month slope plus
    one intercept shared by all 12 months: y = f(month) * diff + c.

    The through-origin fits above force c = 0, so any constant offset
    between CPC's published values and factor * diff is absorbed into the
    slopes. Letting c float separates the two.

    Returns (factors dict {1..12}, offset, residuals).
    """
    n = len(g)
    X = np.zeros((n, 13))
    X[np.arange(n), g.month.values.astype(int) - 1] = g["diff"].values
    X[:, 12] = 1.0
    beta, _, _, _ = np.linalg.lstsq(X, g.ANOM.values, rcond=None)
    resid = g.ANOM.values - X @ beta
    return {m: beta[m - 1] for m in range(1, 13)}, beta[12], resid


def _quantize_match(pred, published):
    """How often each 2-decimal publication convention applied to `pred`
    reproduces `published` exactly, and the RMS it leaves.

    Round-to-nearest and floor (round *down*) are distinguishable here:
    they differ by half a step on average, which is the size of the bias
    left by the through-origin fits.
    """
    eps = 1e-9  # published values are 2-dp decimals read as floats
    out = {}
    for name, q in (("floor", np.floor(pred * 100 + eps) / 100),
                    ("round", np.round(pred, 2)),
                    ("trunc", np.trunc(pred * 100 + eps * np.sign(pred)) / 100)):
        out[f"{name}_exact"] = float(np.mean(np.isclose(q, published, atol=1e-9)))
        out[f"rms_{name}"] = float(np.sqrt(np.mean((q - published) ** 2)))
        out[f"max_{name}"] = float(np.max(np.abs(q - published)))
    return out


def residual_by_decade(factors, weighted=True):
    """Mean and RMS residual per decade, for a fixed set of factors.

    Under a pure publication convention the mean residual would sit at the
    same value (half a step) in every decade. It does not quite: it drifts,
    which is the part of the residual not explained in Section 5 of
    docs/CPC_RONI_scaling.md.
    """
    frame = monthly_fit_frame(weighted=weighted)
    resid = frame.ANOM.values - frame["diff"].values * frame.month.map(factors).values
    frame = frame.assign(resid=resid, decade=(frame.year // 10) * 10)
    return (frame.groupby("decade")
            .agg(n=("resid", "size"), mean_resid=("resid", "mean"),
                 rms_resid=("resid", lambda r: np.sqrt(np.mean(r ** 2))))
            .reset_index())


def rounding_convention_check(n_recent=20, weighted=True):
    """Is the residual left by the direct monthly fit a rounding convention?

    The through-origin fit (backout_monthly_factors_direct()) leaves a small
    *systematic* negative bias — the fit sits ~0.004 above the published
    values on average. Rounding to nearest cannot do that: it is symmetric.
    Flooring can, and leaves exactly half a step (0.005).

    Two diagnostics, each on the 20 most recent years and on the full
    1950-present record:

    - refit with a free intercept (_fit_slopes_with_offset): if the residual
      is a publication convention, the intercept lands near -0.005 (half of
      the 0.01 step), the slopes barely move, and the residual scatter drops
      to the quantization floor (0.01/sqrt(12) = 0.0029).
    - apply each candidate convention to the prediction and count exact
      reproductions of the published 2-decimal value.

    Caveat, stated in docs/CPC_RONI_scaling.md: flooring with no offset and
    round-to-nearest with a genuine -0.005 offset in the anomaly difference
    are observationally identical. The offset landing at exactly half a
    publication step is what favours flooring.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        One row per (sample, factor model) with the fitted offset, raw
        residual statistics and per-convention exact-match rates; and the
        refit (offset-model) factor table.
    """
    frame = monthly_fit_frame(weighted=weighted)
    recent = _recent_years(frame, n_recent)

    fac_origin = {m: np.sum(g["diff"] * g.ANOM) / np.sum(g["diff"] ** 2)
                  for m, g in recent.groupby("month")}
    fac_offset, offset_recent, _ = _fit_slopes_with_offset(recent)

    rows = []
    for sample, g in (("last %d yr" % n_recent, recent), ("1950-present", frame)):
        _, offset_here, resid_off = _fit_slopes_with_offset(g)
        for model, fac in (("through-origin", fac_origin), ("offset refit", fac_offset)):
            pred = g["diff"].values * g.month.map(fac).values
            resid = g.ANOM.values - pred
            row = dict(sample=sample, model=model, n=len(g),
                       rms_raw=np.sqrt(np.mean(resid ** 2)),
                       bias_raw=resid.mean(),
                       offset_fit=offset_here,
                       rms_after_offset=np.sqrt(np.mean(resid_off ** 2)))
            row.update(_quantize_match(pred, g.ANOM.values))
            rows.append(row)

    refit = pd.DataFrame({
        "month": MONTH_NAMES,
        "factor_through_origin": [fac_origin[m] for m in range(1, 13)],
        "factor_offset_refit": [fac_offset[m] for m in range(1, 13)],
    })
    refit["delta"] = refit.factor_offset_refit - refit.factor_through_origin
    return pd.DataFrame(rows), refit


def reconstruct_from_factors(factors, weighted=True):
    """The forward recipe: scale each ERSSTv6 monthly diff by its calendar
    month's factor, then take a 3-month centered running mean of the
    *scaled* values.

    Parameters
    ----------
    factors : dict {1..12: float}

    Returns
    -------
    pandas.DataFrame
        year, month, diff, recon_monthly, recon_3mo, SEAS, YR — the last
        two labelling each 3-month mean by its center month, CPC-style.
    """
    diffs = monthly_diff_lookup(weighted=weighted)
    df = pd.DataFrame([(y, m, v) for (y, m), v in diffs.items()],
                      columns=["year", "month", "diff"]).sort_values(["year", "month"])
    df = df.reset_index(drop=True)
    df["recon_monthly"] = df["diff"] * df.month.map(factors)
    df["recon_3mo"] = df.recon_monthly.rolling(3, center=True).mean()
    df["SEAS"] = df.month.map(_season_label)
    df["YR"] = df.year
    return df


def seasonal_convention_check(factors, weighted=True):
    """The seasonal counterpart of rounding_convention_check(): does the
    published *seasonal* RONI equal the floored, or the rounded, 3-month
    mean of the unrounded scaled monthly values?
    """
    recon = reconstruct_from_factors(factors, weighted=weighted).dropna(subset=["recon_3mo"])
    roni = load_roni_table()
    g = roni.merge(recon, on=["SEAS", "YR"], how="inner", validate="one_to_one")
    resid = g.ANOM.values - g.recon_3mo.values
    row = dict(n=len(g), rms_raw=np.sqrt(np.mean(resid ** 2)), bias_raw=resid.mean())
    row.update(_quantize_match(g.recon_3mo.values, g.ANOM.values))
    return pd.DataFrame([row])


def reconstruction_tables(factors, start_year=2020, weighted=True):
    """Published vs. reconstructed seasonal RONI, in CPC's own year x season
    layout, for the window shown on CPC's RONI product page.

    Returns
    -------
    (published, reconstructed, stats, mismatches_1dp)
        Two pivot tables; a one-row stats frame (n, max|err|, RMS, count
        agreeing at the 1-decimal precision the product page displays); and
        the rows that disagree at 1 decimal.
    """
    recon = reconstruct_from_factors(factors, weighted=weighted).dropna(subset=["recon_3mo"])
    roni = load_roni_table()
    g = roni.merge(recon, on=["SEAS", "YR"], how="inner", validate="one_to_one")
    g = g[g.YR >= start_year].copy()
    g["err"] = g.recon_3mo - g.ANOM
    g["agree_1dp"] = np.round(g.recon_3mo, 1) == np.round(g.ANOM, 1)

    season_order = [_season_label(m) for m in range(1, 13)]
    published = g.pivot(index="YR", columns="SEAS", values="ANOM")[season_order]
    reconstructed = (g.pivot(index="YR", columns="SEAS", values="recon_3mo")[season_order]
                     .round(2))
    for t in (published, reconstructed):
        t.index.name = "Year"

    stats = pd.DataFrame([dict(
        n=len(g), max_abs_err=g.err.abs().max(),
        rms=np.sqrt(np.mean(g.err ** 2)),
        n_agree_1dp=int(g.agree_1dp.sum()),
        n_recon_above_published=int((g.err > 0.005).sum()),
        n_recon_below_published=int((g.err < -0.005).sum()),
    )])
    mismatches = g.loc[~g.agree_1dp, ["SEAS", "YR", "ANOM", "recon_3mo"]]
    return published, reconstructed, stats, mismatches


def weighting_sensitivity(n_recent=20):
    """Cosine-latitude weighted vs. plain grid-cell box means, side by side.

    Reported because the reference implementation this recipe was checked
    against uses unweighted box means: the choice barely changes how well
    the fit reproduces the published series, but it does move the fitted
    factor itself by more than that factor's standard error.
    """
    w_tab, _ = backout_monthly_factors_direct(n_recent=n_recent, weighted=True)
    u_tab, _ = backout_monthly_factors_direct(n_recent=n_recent, weighted=False)
    return pd.DataFrame({
        "month": w_tab.month,
        "factor_weighted": w_tab.factor,
        "factor_unweighted": u_tab.factor,
        "delta": u_tab.factor - w_tab.factor,
        "se_weighted": w_tab.se,
        "rms_weighted": w_tab.rms_resid,
        "rms_unweighted": u_tab.rms_resid,
    })


def factor_table_by_year(merged=None):
    """Per-(year, season) backed-out factor = published RONI / ERSSTv6
    diff_3mo, pivoted into the same year x season layout as CPC's own
    RONI.ascii.txt table / roni web page — for visually inspecting where
    the factor is stable versus where a near-zero diff_3mo year makes the
    single-row ratio blow up (see backout_factors()'s factor_mean_all_years
    diagnostic for the same instability, aggregated instead of per-year).
    """
    if merged is None:
        merged = load_merged()
    merged = merged.copy()
    merged["factor"] = merged.ANOM / merged.diff_3mo
    season_order = [_season_label(m) for m in range(1, 13)]
    table = merged.pivot(index="YR", columns="SEAS", values="factor")[season_order]
    table.index.name = "Year"
    return table


def _emit(text, log):
    print(text)
    log.append(text)


if __name__ == "__main__":
    pd.set_option("display.float_format", lambda v: f"{v:8.4f}")
    pd.set_option("display.width", 120)

    out_dir = config.PLOTS_DIR_RONI_FACTOR_BACKOUT
    out_dir.mkdir(parents=True, exist_ok=True)
    log = []

    # Step 1 (docs/CPC_RONI_scaling.md Section 2): do CPC's own published
    # monthly and seasonal series agree with each other, with no ERSSTv6
    # involved at all? Establishes that the scaling factor is applied
    # monthly, before any temporal averaging.
    consistency = check_monthly_seasonal_consistency()
    _emit("\n-- Step 1: published seasonal RONI vs. 3-month avg of published monthly RONI --", log)
    _emit(f"n={len(consistency)}, max|diff|={consistency['diff'].abs().max():.4f}, "
          f"rms={np.sqrt((consistency['diff']**2).mean()):.4f}, "
          f"mean={consistency['diff'].mean():.6f}", log)
    consistency.to_csv(out_dir / "monthly_seasonal_consistency.csv", index=False)

    merged = load_merged()  # weighted=True (cos-lat area weighting) — the
    # unweighted variant tested 2026-09-10 gave no improvement (RMS 0.00488
    # vs 0.00482) and is not the default.

    table = backout_factors(merged=merged)
    _emit("\n-- single-row backout, validated against all other years --", log)
    _emit(table.to_string(index=False), log)
    table.to_csv(out_dir / "seasonal_factors_singlerow.csv", index=False)

    recent = backout_factors_recent(n_recent=20, merged=merged)
    _emit("\n-- OLS-through-origin on the 20 most recent years, with accuracy --", log)
    _emit(recent.to_string(index=False), log)
    recent.to_csv(out_dir / "seasonal_factors_recent.csv", index=False)

    monthly, monthly_rms, monthly_n = backout_monthly_factors(n_recent=20, merged=merged)
    _emit(f"\n-- joint monthly-factor fit, 20 most recent years, n_obs={monthly_n} --", log)
    _emit(monthly.to_string(index=False), log)
    _emit(f"overall RMS residual: {monthly_rms:.4f}  "
          f"(cf. per-season RMS 0.0041-0.0128 from the season-block fit above)", log)
    monthly.to_csv(out_dir / "monthly_factors.csv", index=False)

    # Cross-check: does treating the factor as constant within each season
    # (the old block-per-season fit) actually cost accuracy relative to the
    # joint monthly fit? Compare season-by-season RMS under both models.
    beta_by_month = dict(zip(range(1, 13), monthly.factor.values))
    diffs = monthly_diff_lookup()
    rows = []
    season_order = [_season_label(m) for m in range(1, 13)]
    cutoff_year = merged.YR.max() - 20 + 1
    for seas in season_order:
        g = merged[(merged.SEAS == seas) & (merged.YR >= cutoff_year)]
        pred = []
        for row in g.itertuples():
            p = sum(beta_by_month[_neighbor_month(row.year, row.month, off)[1]]
                     * diffs[_neighbor_month(row.year, row.month, off)] / 3.0
                     for off in (-1, 0, 1))
            pred.append(p)
        resid_monthly = g.ANOM.values - np.array(pred)
        rows.append(dict(
            SEAS=seas,
            rms_seasonal_fit=recent[recent.SEAS == seas].rms_resid.values[0],
            rms_monthly_fit=np.sqrt(np.mean(resid_monthly ** 2)),
        ))
    compare = pd.DataFrame(rows)
    _emit("\n-- per-season RMS: seasonal-block fit vs. joint monthly fit --", log)
    _emit(compare.to_string(index=False), log)
    compare.to_csv(out_dir / "seasonal_vs_monthly_rms.csv", index=False)

    # Third check (2026-09-10): CPC also publishes the *monthly* relative
    # Nino-3.4 series directly (Rnino34.ascii.txt) — fit each month's factor
    # against that directly (no seasonal-aggregate coupling needed), then
    # verify the full round trip: scale monthly, 3-month-average, compare
    # to the independently-published seasonal RONI.
    monthly_direct, rnino34_matched = backout_monthly_factors_direct(n_recent=20)
    _emit(f"\n-- monthly factor fit DIRECT against Rnino34.ascii.txt, "
          f"20 most recent years --", log)
    _emit(monthly_direct.to_string(index=False), log)
    monthly_direct.to_csv(out_dir / "monthly_factors_direct.csv", index=False)

    # Pooled accuracy of that fit, over every fitted (month, year) point —
    # the summary quoted in docs/CPC_RONI_scaling.md §3. The bias is
    # systematic, not sampling noise; §5 there identifies its source.
    factors_direct = dict(zip(range(1, 13), monthly_direct.factor.values))
    fitted = _recent_years(monthly_fit_frame(), 20)
    resid_pooled = fitted.ANOM.values - fitted["diff"].values * fitted.month.map(factors_direct).values
    _emit(f"pooled over n={len(resid_pooled)} fitted points: "
          f"rms={np.sqrt(np.mean(resid_pooled ** 2)):.4f}, "
          f"bias={resid_pooled.mean():+.4f} "
          f"(published minus fit; the fit sits above the published values)", log)

    seasonal_check, seasonal_check_rms = validate_direct_monthly_against_seasonal(
        monthly_direct, merged, n_recent=20)
    _emit(f"\n-- round-trip check: direct monthly factors, 3-month averaged, "
          f"vs. published seasonal RONI --", log)
    _emit(seasonal_check.to_string(index=False), log)
    _emit(f"overall RMS: {seasonal_check_rms:.4f}  "
          f"(cf. {monthly_rms:.4f} for the joint-fit-from-seasonal-only approach)", log)
    seasonal_check.to_csv(out_dir / "direct_monthly_vs_seasonal_rms.csv", index=False)

    # Cos-lat weighting vs. plain grid-cell means (docs §1): immaterial for
    # fit quality, but not for the factor's own value.
    weighting = weighting_sensitivity(n_recent=20)
    _emit("\n-- box-mean weighting sensitivity (cos-lat vs. plain grid-cell) --", log)
    _emit(weighting.to_string(index=False), log)
    _emit(f"max |delta factor| = {weighting.delta.abs().max():.4f}, "
          f"max |delta rms| = {(weighting.rms_unweighted - weighting.rms_weighted).abs().max():.5f}", log)
    weighting.to_csv(out_dir / "weighting_sensitivity.csv", index=False)

    # Publication convention (docs §5): the residual bias above is CPC
    # flooring its published values to 2 decimals, not rounding them.
    conv, refit = rounding_convention_check(n_recent=20)
    _emit("\n-- publication convention: monthly series --", log)
    _emit(conv.to_string(index=False), log)
    _emit("\n-- factors refit with a free intercept vs. through the origin --", log)
    _emit(refit.to_string(index=False), log)
    conv.to_csv(out_dir / "rounding_convention.csv", index=False)
    refit.to_csv(out_dir / "monthly_factors_offset_refit.csv", index=False)

    factors_refit = dict(zip(range(1, 13), refit.factor_offset_refit.values))
    seas_conv = pd.concat([
        seasonal_convention_check(factors_direct).assign(model="through-origin"),
        seasonal_convention_check(factors_refit).assign(model="offset refit"),
    ])
    _emit("\n-- publication convention: seasonal series --", log)
    _emit(seas_conv.to_string(index=False), log)
    seas_conv.to_csv(out_dir / "rounding_convention_seasonal.csv", index=False)

    decades = residual_by_decade(factors_refit)
    _emit("\n-- residual by decade, refit factors (docs §7: the drift is unexplained) --", log)
    _emit(decades.to_string(index=False), log)
    decades.to_csv(out_dir / "residual_by_decade.csv", index=False)

    # Round trip, in CPC's own table layout (docs §3, Tables 2-3).
    published_tab, recon_tab, recon_stats, recon_miss = reconstruction_tables(
        factors_direct, start_year=2020)
    two_dp = lambda v: f"{v:6.2f}"
    _emit("\n-- published seasonal RONI, 2020-present --", log)
    _emit(published_tab.to_string(float_format=two_dp), log)
    _emit("\n-- reconstruction from ERSSTv6 + the direct monthly factors --", log)
    _emit(recon_tab.to_string(float_format=two_dp), log)
    _emit(recon_stats.to_string(index=False), log)
    _emit("rows disagreeing at the 1-decimal precision of CPC's product page:", log)
    _emit(recon_miss.to_string(index=False), log)
    published_tab.to_csv(out_dir / "seasonal_published_recent.csv")
    recon_tab.to_csv(out_dir / "seasonal_reconstructed_recent.csv")
    recon_stats.to_csv(out_dir / "seasonal_reconstruction_stats.csv", index=False)

    by_year = factor_table_by_year(merged=merged)
    by_year.to_csv(out_dir / "factor_by_year.csv")
    pd.set_option("display.max_rows", None)
    pd.set_option("display.float_format", lambda v: f"{v:6.2f}")
    _emit("\n-- per-(year, season) backed-out factor (RONI / ERSSTv6 diff_3mo) --", log)
    _emit(by_year.to_string(), log)

    report_path = out_dir / "report.txt"
    report_path.write_text("\n".join(log) + "\n")
    print(f"\nwrote {report_path} and CSVs to {out_dir}/")
