#!/usr/bin/env python3
"""
BSISO Bimodal Index from NOAA OLR Anomalies
===========================================

Author
------
Sandro W. Lubis
Pacific Northwest National Laboratory (PNNL)

Purpose
-------
This script calculates the boreal summer intraseasonal oscillation (BSISO)
bimodal index from daily NOAA outgoing longwave radiation (OLR) anomalies.

The workflow:
    1. Reads daily NOAA OLR anomalies.
    2. Selects the tropical domain 30°S-30°N.
    3. Applies a 25-90 day Lanczos band-pass filter with 141 weights.
    4. Constructs extended EOFs (EEOFs) using lags -10, -5, and 0 days.
    5. Calculates and standardizes the leading two principal components.
    6. Derives BSISO amplitude and phases 1-8.
    7. Saves the EEOF and BSISO index to NetCDF.
    8. Produces an EEOF verification plot centered at 120°E.


Citation
--------
If you use this script or the BSISO index generated from it, please cite:

Lubis, S. W., Chen, Z., Lu, J., Hagos, S., Chang, C.-C., & Leung, L. R.
(2024). Enhanced Pacific Northwest heat extremes and wildfire risks induced
by the boreal summer intraseasonal oscillation. npj Climate and Atmospheric
Science, 7, 232. https://doi.org/10.1038/s41612-024-00766-3
"""

from pathlib import Path
from datetime import datetime
import os

import numpy as np
import xarray as xr
from scipy.ndimage import convolve1d
from scipy.sparse.linalg import svds

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point


# =============================================================================
# USER SETTINGS
# =============================================================================

# ONLY input file used by this script
OLR_FILE = "/global/homes/s/slubis/link2/OLR/olr.anom.1979_2021.nc"
OLR_VAR = "olr"

# Analysis period.
# Filtering uses the full 1979-2021 file as temporal padding, but the BSISO
# index itself is produced for 1980-2020.
YEAR_START = 1980
YEAR_END = 2020

# Tropical EEOF domain
LAT_MIN = -30.0
LAT_MAX = 30.0

# EEOF training season June-October (JJASO).
TRAIN_MONTHS = (6, 7, 8, 9, 10)

# EEOF lags
LAGS = np.array([-10, -5, 0], dtype=np.int32)
NEOF = 2

# Lanczos band-pass filter
FILTER_PERIOD_LOW = 25.0    # days
FILTER_PERIOD_HIGH = 90.0   # days
FILTER_NWEIGHTS = 141       # must be odd
FILTER_NSIGMA = 1.0

# BSISO active threshold.
# Original supplied NCL uses amplitude > 1 (<=1 is inactive).
AMP_THRESHOLD = 1.0

# Sign convention: the regional-average EEOF pattern to choose a consistent sign orientation.
SIGN_METHOD = "mean"

# Optional: save the filtered 1980-2020 OLR field.
# This file can be relatively large.
SAVE_FILTERED_OLR = False

# Outputs
FILTERED_FILE = "./olr.anom.25_90d.1980_2020.nc"
EEOF_FILE = "./cal_eeof_olr_inReanalysis.nc"
PHASE_FILE = "./bsiso_phase.nc"
CHECK_FIG = "./bsiso_eeof_check_120E.png"

# Plot settings
CENTER_LON = 120.0
PLOT_LEVELS = [-10, -8, -6, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 6, 8, 10]
PLOT_TICKS = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
PLOT_CMAP = "RdBu_r"

# Use float32 for large gridded arrays
DTYPE = np.float32


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def subset_lat(da, lat_min, lat_max):
    """Select latitude range regardless of coordinate order."""
    return da.where(
        (da["lat"] >= lat_min) & (da["lat"] <= lat_max),
        drop=True
    )


def check_no_missing(a, name):
    """Require complete fields for the EEOF calculation."""
    if not np.all(np.isfinite(a)):
        nbad = np.size(a) - np.count_nonzero(np.isfinite(a))
        raise ValueError(
            f"{name} contains {nbad} NaN/Inf/missing values. "
            "The EEOF calculation requires complete data."
        )


def check_daily_time(time):
    """Make sure the NOAA input is daily and continuous."""
    vals = np.asarray(time.values)

    if vals.size < 2:
        raise ValueError("Time coordinate is too short.")

    # NOAA standard calendar normally decodes to numpy datetime64.
    if np.issubdtype(vals.dtype, np.datetime64):
        dt_hours = np.diff(vals).astype("timedelta64[h]").astype(np.int64)
        if not np.all(dt_hours == 24):
            bad = np.where(dt_hours != 24)[0]
            raise ValueError(
                f"Input time is not continuous daily data. "
                f"Found {bad.size} non-24-hour interval(s)."
            )


# =============================================================================
# 25-90 DAY LANCZOS BAND-PASS FILTER
# =============================================================================

def lanczos_bandpass_weights(
    nweights=141,
    period_low=25.0,
    period_high=90.0,
    nsigma=1.0,
):
    """
    Construct a standard Lanczos band-pass FIR filter.

    For daily data:
        low-frequency cutoff  = 1 / period_high
        high-frequency cutoff = 1 / period_low

    This follows the standard ideal band-pass kernel multiplied by the
    Lanczos sigma factor. nsigma=1 corresponds to the commonly used
    Lanczos window.

    Parameters
    ----------
    nweights : int
        Odd number of weights.
    period_low : float
        Short-period cutoff in days (25).
    period_high : float
        Long-period cutoff in days (90).
    nsigma : float
        Power applied to the Lanczos sigma factor.

    Returns
    -------
    weights : ndarray
        Symmetric FIR weights.
    """
    if nweights < 3 or nweights % 2 == 0:
        raise ValueError("FILTER_NWEIGHTS must be an odd integer >= 3.")

    if not (0 < period_low < period_high):
        raise ValueError("Need 0 < period_low < period_high.")

    m = (nweights - 1) // 2
    k = np.arange(-m, m + 1, dtype=np.float64)

    # Cycles per day
    fca = 1.0 / period_high   # low-frequency cutoff
    fcb = 1.0 / period_low    # high-frequency cutoff

    # Ideal band-pass = low-pass(fcb) - low-pass(fca)
    # np.sinc(x) = sin(pi*x)/(pi*x)
    ideal = (
        2.0 * fcb * np.sinc(2.0 * fcb * k)
        - 2.0 * fca * np.sinc(2.0 * fca * k)
    )

    # Lanczos sigma factor. End weights go to zero.
    sigma = np.sinc(k / m)
    sigma = sigma ** nsigma

    weights = ideal * sigma

    return weights.astype(DTYPE)


def bandpass_filter_olr(olr):
    """
    Apply centered 25-90-day Lanczos filtering along time.

    The full NOAA 1979-2021 record is filtered first so that 1980-2020
    has adequate temporal padding. End points affected by the half-window
    are set to NaN and are never used in the 1980-2020 analysis.
    """
    print("\nApplying 25-90 day Lanczos band-pass filter...")
    print(f"  Number of weights = {FILTER_NWEIGHTS}")
    print(f"  Half window       = {(FILTER_NWEIGHTS - 1)//2} days")

    weights = lanczos_bandpass_weights(
        nweights=FILTER_NWEIGHTS,
        period_low=FILTER_PERIOD_LOW,
        period_high=FILTER_PERIOD_HIGH,
        nsigma=FILTER_NSIGMA,
    )

    data = np.asarray(olr.values, dtype=DTYPE)
    check_no_missing(data, "Raw NOAA OLR anomaly")

    # Centered symmetric convolution along time.
    filtered = convolve1d(
        data,
        weights=weights,
        axis=0,
        mode="constant",
        cval=0.0,
    ).astype(DTYPE, copy=False)

    # Match the behavior of a running Lanczos filter: no valid output where
    # the complete window is unavailable.
    half = (FILTER_NWEIGHTS - 1) // 2
    filtered[:half, :, :] = np.nan
    filtered[-half:, :, :] = np.nan

    da = xr.DataArray(
        filtered,
        dims=("time", "lat", "lon"),
        coords={
            "time": olr["time"],
            "lat": olr["lat"],
            "lon": olr["lon"],
        },
        name="olr_25_90",
        attrs={
            "long_name": "25-90 day Lanczos band-pass filtered OLR anomaly",
            "units": olr.attrs.get("units", "W/m^2"),
            "filter": "Lanczos band-pass",
            "period": "25-90 days",
            "nweights": FILTER_NWEIGHTS,
            "nsigma": FILTER_NSIGMA,
        },
    )

    return da, weights


# =============================================================================
# EEOF CONSTRUCTION
# =============================================================================

def prepare_training_data(filtered_olr):
    """
    Select training months and reshape to:
        year x seasonal_day x lat x lon
    """
    year = filtered_olr["time"].dt.year
    month = filtered_olr["time"].dt.month

    mask = (
        year.isin(np.arange(YEAR_START, YEAR_END + 1))
        & month.isin(TRAIN_MONTHS)
    )

    train = filtered_olr.where(mask, drop=True)

    years = np.asarray(train["time"].dt.year.values)
    unique_years, counts = np.unique(years, return_counts=True)
    expected_years = np.arange(YEAR_START, YEAR_END + 1)

    if not np.array_equal(unique_years, expected_years):
        raise ValueError(
            "Training period contains missing years. "
            f"Expected {YEAR_START}-{YEAR_END}."
        )

    if not np.all(counts == counts[0]):
        bad = dict(zip(unique_years, counts))
        raise ValueError(
            "Training seasons do not all contain the same number of days: "
            f"{bad}"
        )

    nseason = int(counts[0])
    nyears = unique_years.size
    nlat = train.sizes["lat"]
    nlon = train.sizes["lon"]

    arr = np.asarray(train.values, dtype=DTYPE)
    check_no_missing(arr, "Filtered OLR training data")

    arr = arr.reshape(nyears, nseason, nlat, nlon)

    print("\nEEOF training data:")
    print(f"  years          = {YEAR_START}-{YEAR_END} ({nyears})")
    print(f"  months         = {TRAIN_MONTHS}")
    print(f"  days/season    = {nseason}")
    print(f"  latitude       = {nlat}")
    print(f"  longitude      = {nlon}")

    return arr, nseason


def compute_eeof(train_4d, lags=LAGS, neof=NEOF):
    """
    Compute covariance EEOFs.

    Instead of constructing the NCL 3*nlat x 3*nlon sparse block array,
    concatenate the three populated lag blocks directly:

       [lag -10 spatial field | lag -5 field | lag 0 field]

    This is mathematically equivalent but uses much less memory.
    """
    nyears, ndays, nlat, nlon = train_4d.shape

    min_lag = int(np.min(lags))
    nsample_days = ndays + min_lag

    if nsample_days <= 0:
        raise ValueError("Season is too short for requested lags.")

    offsets = lags - min_lag
    nspace = nlat * nlon

    blocks = []

    for offset in offsets:
        start = int(offset)
        stop = start + nsample_days

        block = train_4d[:, start:stop, :, :].reshape(
            nyears * nsample_days,
            nspace
        )

        blocks.append(block)

    X = np.concatenate(blocks, axis=1).astype(DTYPE, copy=False)
    check_no_missing(X, "EEOF matrix")

    print("\nEEOF matrix:")
    print(f"  samples  = {X.shape[0]}")
    print(f"  features = {X.shape[1]}")

    # Covariance EOF: remove the temporal mean at each EEOF feature.
    X -= X.mean(axis=0, dtype=np.float64).astype(DTYPE)

    # Total covariance variance for percentage variance.
    total_variance = float(
        np.var(X, axis=0, ddof=1, dtype=np.float64).sum()
    )

    print("\nCalculating leading EEOFs...")

    # X = U S V^T. Rows of V^T are unit-normalized EOF vectors.
    _, singular_values, vt = svds(
        X,
        k=neof,
        which="LM",
    )

    # svds returns singular values in ascending order.
    order = np.argsort(singular_values)[::-1]
    singular_values = singular_values[order]
    vt = vt[order, :]

    eigenvalues = (
        singular_values.astype(np.float64) ** 2
        / (X.shape[0] - 1)
    )

    pcvar = 100.0 * eigenvalues / total_variance

    eof = vt.reshape(
        neof,
        len(lags),
        nlat,
        nlon
    ).astype(DTYPE)

    print("  eigenvalues =", eigenvalues)
    print("  variance %  =", pcvar)

    return eof, eigenvalues, pcvar


# =============================================================================
# EEOF SIGN CONVENTION
# =============================================================================

def region_values(field2d, lat, lon, lat1, lat2, lon1, lon2):
    lat_mask = (
        (lat >= min(lat1, lat2))
        & (lat <= max(lat1, lat2))
    )

    # NOAA longitude is expected to be 0-360.
    lon360 = np.mod(lon, 360.0)
    lon1 = lon1 % 360.0
    lon2 = lon2 % 360.0

    if lon1 <= lon2:
        lon_mask = (lon360 >= lon1) & (lon360 <= lon2)
    else:
        lon_mask = (lon360 >= lon1) | (lon360 <= lon2)

    if not np.any(lat_mask):
        raise ValueError(
            f"No latitude points found in region {lat1}:{lat2}."
        )

    if not np.any(lon_mask):
        raise ValueError(
            f"No longitude points found in region {lon1}:{lon2}."
        )

    return field2d[np.ix_(lat_mask, lon_mask)]


def enforce_sign_convention(eof, lat, lon):
    """
    Orient EEOF1/EEOF2 consistently with the supplied NCL script.

    The supplied NCL script used ANY grid point to decide whether to reverse
    the sign, which can be unstable. By default this optimized version uses
    the REGIONAL MEAN sign instead.

    Set SIGN_METHOD = "legacy_any" above for a literal reproduction.
    """
    eof = eof.copy()

    # Lag index 0 = day -10
    e1_a = region_values(eof[0, 0], lat, lon, 12, 20, 72, 84)
    e1_b = region_values(eof[0, 0], lat, lon, -2.5, 0.5, 75, 85)

    e2_a = region_values(eof[1, 0], lat, lon, -2.5, 5, 80, 90)
    e2_b = region_values(eof[1, 0], lat, lon, 12.5, 17.5, 128, 132)

    if SIGN_METHOD == "legacy_any":
        reverse_eof1 = np.any(e1_a < 0) or np.any(e1_b > 0)
        reverse_eof2 = np.any(e2_a < 0) or np.any(e2_b > 0)

    elif SIGN_METHOD == "mean":
        reverse_eof1 = (np.nanmean(e1_a) < 0) or (np.nanmean(e1_b) > 0)
        reverse_eof2 = (np.nanmean(e2_a) < 0) or (np.nanmean(e2_b) > 0)

    else:
        raise ValueError(
            "SIGN_METHOD must be 'mean' or 'legacy_any'."
        )

    if reverse_eof1:
        eof[0] *= -1.0
        print("  Reversed EEOF1 sign.")

    if reverse_eof2:
        eof[1] *= -1.0
        print("  Reversed EEOF2 sign.")

    return eof


# =============================================================================
# PROJECT FULL 1980-2020 RECORD ONTO EEOFS
# =============================================================================

def project_full_record(filtered_olr, eof, lags=LAGS, chunk_size=512):
    """
    Project the complete filtered 1980-2020 record onto the EEOFs.

    With lags [-10,-5,0], output date t uses:
        OLR(t-10), OLR(t-5), OLR(t)

    Therefore an input record of 1980-01-01 ... 2020-12-31 produces
    PCs from 1980-01-11 ... 2020-12-31 = 14966 days.
    """
    ntime = filtered_olr.sizes["time"]
    min_lag = int(np.min(lags))
    nout = ntime + min_lag

    if nout <= 0:
        raise ValueError("Analysis record too short for requested lags.")

    offsets = lags - min_lag
    neof = eof.shape[0]

    pcs = np.empty((neof, nout), dtype=np.float64)

    print("\nProjecting full 1980-2020 record onto EEOFs...")

    for start in range(0, nout, chunk_size):
        end = min(start + chunk_size, nout)
        score = np.zeros((neof, end - start), dtype=np.float64)

        for ilag, offset in enumerate(offsets):
            field = np.asarray(
                filtered_olr.isel(
                    time=slice(
                        start + int(offset),
                        end + int(offset)
                    )
                ).values,
                dtype=DTYPE,
            )

            check_no_missing(
                field,
                f"projection field at lag {lags[ilag]}"
            )

            score += np.einsum(
                "tij,mij->mt",
                field,
                eof[:, ilag, :, :],
                optimize=True,
            )

        pcs[:, start:end] = score

    # NCL eofunc_ts_n subtracts the mean of each PC time series.
    pcs -= pcs.mean(axis=1, keepdims=True)

    raw_pc_std = pcs.std(axis=1, ddof=0)

    if np.any(raw_pc_std == 0):
        raise ValueError("A PC has zero standard deviation.")

    # NCL dim_standardize(..., opt=1): population standard deviation.
    pcs /= raw_pc_std[:, None]

    pc_time = filtered_olr["time"].isel(
        time=slice(-min_lag, None)
    )

    print(f"  PC start = {pc_time.values[0]}")
    print(f"  PC end   = {pc_time.values[-1]}")
    print(f"  Ntime    = {pc_time.size}")

    return pcs.astype(DTYPE), pc_time, raw_pc_std


# =============================================================================
# BSISO PHASE
# =============================================================================

def calculate_bsiso_phase(pcs):
    """Calculate standardized-PC amplitude, angle, and phases 1-8."""
    pc1 = pcs[0]
    pc2 = pcs[1]

    amplitude = np.sqrt(pc1**2 + pc2**2)

    # Exact quadrant-aware equivalent of my_atan(-PC1, -PC2).
    theta = np.degrees(
        np.arctan2(-pc1, -pc2)
    )

    # Guarantee [-180, 180)
    theta = (theta + 180.0) % 360.0 - 180.0

    # Eight 45-degree sectors:
    # [-180,-135) -> 1, ..., [135,180) -> 8
    phase = (
        np.floor((theta + 180.0) / 45.0).astype(np.int32) + 1
    )

    phase = np.clip(phase, 1, 8)

    # Match the supplied NCL: amplitude <= 1 is inactive.
    phase_active = phase.copy()
    phase_active[amplitude <= AMP_THRESHOLD] = 0

    return (
        amplitude.astype(DTYPE),
        theta.astype(DTYPE),
        phase.astype(np.int32),
        phase_active.astype(np.int32),
    )


# =============================================================================
# PHYSICAL EEOF PATTERNS
# =============================================================================

def physical_eeof_patterns(eof, eigenvalues):
    """
    Convert unit-normalized EOF vectors to physical loading patterns.

    covariance EOF loading = EOF * sqrt(eigenvalue)

    Units are therefore the same as OLR: W m-2.
    """
    scale = np.sqrt(eigenvalues).astype(DTYPE)

    return (
        eof
        * scale[:, None, None, None]
    ).astype(DTYPE)


# =============================================================================
# CHECK PLOT
# =============================================================================

def make_check_plot(eof_physical, pcvar, lat, lon):
    """
    Create a polished 3x2 EEOF check figure.

    Layout
    ------
    Left column : EEOF2
    Right column: EEOF1

    Rows:
      day -10
      day  -5
      day   0

    Improvements
    ------------
    * Centered at 120E.
    * Longitude is normalized to 0-360 and sorted before plotting.
    * A cyclic longitude point is added to remove vertical seam artifacts.
    * Wider +/-10 W m-2 plotting range reduces color saturation.
    * Shared colorbar is placed well below the bottom row.
    * More compact figure aspect reduces unnecessary vertical whitespace.
    """
    print(f"\nMaking check plot: {CHECK_FIG}")

    projection = ccrs.PlateCarree(
        central_longitude=CENTER_LON
    )
    data_crs = ccrs.PlateCarree()

    # ------------------------------------------------------------------
    # Normalize longitude before plotting.
    #
    # This is important if the input happens to use -180..180 instead of
    # 0..360.  Sorting prevents an interior longitude jump, and the cyclic
    # point then closes the final 360/0 boundary cleanly.
    # ------------------------------------------------------------------
    lon360 = np.mod(np.asarray(lon, dtype=np.float64), 360.0)
    lon_order = np.argsort(lon360)
    lon_plot = lon360[lon_order]

    # Six wide tropical panels; a shorter figure avoids large empty gaps.
    fig, axes = plt.subplots(
        3,
        2,
        figsize=(15.5, 7.0),
        subplot_kw={"projection": projection},
        constrained_layout=False,
    )

    # Match the reference layout:
    # left column = EEOF2, right column = EEOF1
    mode_order = [1, 0]

    panel_left = ["(a)", "(b)", "(c)"]
    panel_right = ["(d)", "(e)", "(f)"]

    # A 120E-centered global map has its seam at 60W.
    xticks = [300, 0, 60, 120, 180, 240]
    yticks = [-30, -15, 0, 15, 30]

    cf = None

    for row, lag in enumerate(LAGS):
        for col, mode in enumerate(mode_order):
            ax = axes[row, col]

            # Reorder longitude into monotonically increasing 0..360.
            field = np.asarray(
                eof_physical[mode, row, :, :],
                dtype=np.float64,
            )[:, lon_order]

            # ----------------------------------------------------------
            # Add cyclic longitude.
            #
            # This closes the 360 -> 0 boundary and prevents the narrow
            # white vertical seam that contourf can otherwise leave.
            # ----------------------------------------------------------
            field_cyclic, lon_cyclic = add_cyclic_point(
                field,
                coord=lon_plot,
                axis=-1,
            )

            cf = ax.contourf(
                lon_cyclic,
                lat,
                field_cyclic,
                levels=PLOT_LEVELS,
                cmap=PLOT_CMAP,
                extend="both",
                transform=data_crs,
                antialiased=False,
            )

            ax.coastlines(
                resolution="110m",
                linewidth=0.75,
            )

            # Tropical global domain; 120E is physically at panel center.
            ax.set_global()
            ax.set_ylim(LAT_MIN, LAT_MAX)

            ax.set_xticks(
                xticks,
                crs=data_crs,
            )
            ax.set_yticks(
                yticks,
                crs=data_crs,
            )

            ax.xaxis.set_major_formatter(
                LongitudeFormatter(
                    degree_symbol="°",
                    zero_direction_label=False,
                    dateline_direction_label=False,
                )
            )

            ax.yaxis.set_major_formatter(
                LatitudeFormatter(
                    degree_symbol="°"
                )
            )

            ax.tick_params(
                axis="both",
                labelsize=10,
                length=4,
                width=0.8,
                pad=3,
            )

            # Slightly stronger map frame.
            try:
                ax.spines["geo"].set_linewidth(1.0)
            except Exception:
                pass

            panel = (
                panel_left[row]
                if col == 0
                else panel_right[row]
            )

            # Panel/mode title
            ax.text(
                0.00,
                1.035,
                f"{panel} BSISO mode: EEOF{mode + 1} "
                f"({pcvar[mode]:.2f}%)",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=13,
            )

            # Lag label
            ax.text(
                1.00,
                1.035,
                f"day {int(lag)}",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=13,
            )

    # ------------------------------------------------------------------
    # Reserve explicit space for a colorbar below all six panels.
    # ------------------------------------------------------------------
    fig.subplots_adjust(
        left=0.055,
        right=0.985,
        top=0.955,
        bottom=0.185,
        wspace=0.10,
        hspace=0.43,
    )

    # Dedicated colorbar axis:
    # [left, bottom, width, height] in figure coordinates.
    cbar_ax = fig.add_axes(
        [0.16, 0.065, 0.68, 0.028]
    )

    cbar = fig.colorbar(
        cf,
        cax=cbar_ax,
        orientation="horizontal",
        extend="both",
    )

    cbar.set_ticks(PLOT_TICKS)
    cbar.ax.tick_params(
        labelsize=10,
        length=3,
        pad=3,
    )

    cbar.set_label(
        r"EEOF OLR loading [W m$^{-2}$]",
        fontsize=12,
        labelpad=5,
    )

    fig.savefig(
        CHECK_FIG,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.06,
    )

    plt.close(fig)


# =============================================================================
# NETCDF OUTPUT
# =============================================================================

def write_outputs(
    filtered_analysis,
    eof,
    eof_physical,
    eigenvalues,
    pcvar,
    pcs,
    pc_time,
    amplitude,
    theta,
    phase,
    phase_active,
):
    """Write EEOF and BSISO output NetCDF files."""
    history = (
        f"{datetime.now().isoformat()} "
        f"{os.getcwd()}/cal_bsiso_from_noaa_olr.py"
    )

    time_encoding = {
        "units": "hours since 1800-01-01 00:00:0.0",
        "calendar": "standard",
        "dtype": "float64",
    }

    # -------------------------------------------------------------------------
    # Optional filtered OLR
    # -------------------------------------------------------------------------
    if SAVE_FILTERED_OLR:
        print(f"\nWriting filtered OLR: {FILTERED_FILE}")

        ds_filtered = xr.Dataset(
            {
                "olr_25_90": filtered_analysis
            },
            attrs={
                "history": history,
                "source": OLR_FILE,
            },
        )

        Path(FILTERED_FILE).unlink(missing_ok=True)

        ds_filtered.to_netcdf(
            FILTERED_FILE,
            encoding={
                "olr_25_90": {
                    "dtype": "float32",
                    "zlib": True,
                    "complevel": 2,
                },
                "time": time_encoding,
            },
        )

    # -------------------------------------------------------------------------
    # EEOF file
    # -------------------------------------------------------------------------
    print(f"Writing EEOF file: {EEOF_FILE}")

    ds_eeof = xr.Dataset(
        data_vars={
            "eof": (
                ("evn", "lag_lead", "lat", "lon"),
                eof,
                {
                    "long_name": "unit-normalized covariance EEOF",
                    "note": "sum of squares of each EEOF is approximately one",
                },
            ),
            "eof_physical": (
                ("evn", "lag_lead", "lat", "lon"),
                eof_physical,
                {
                    "long_name": "physical EEOF loading",
                    "units": filtered_analysis.attrs.get("units", "W/m^2"),
                    "definition": "eof * sqrt(eigenvalue)",
                },
            ),
            "ev_ts": (
                ("evn", "time"),
                pcs,
                {
                    "long_name": "standardized EEOF principal component",
                    "units": "1",
                },
            ),
            "eigenvalue": (
                ("evn",),
                eigenvalues.astype(np.float64),
            ),
            "pcvar": (
                ("evn",),
                pcvar.astype(np.float64),
                {
                    "long_name": "percent variance explained",
                    "units": "%",
                },
            ),
        },
        coords={
            "evn": np.arange(1, NEOF + 1, dtype=np.int32),
            "lag_lead": LAGS.astype(np.int32),
            "lat": filtered_analysis["lat"],
            "lon": filtered_analysis["lon"],
            "time": pc_time,
        },
        attrs={
            "history": history,
            "source_olr": OLR_FILE,
            "olr_variable": OLR_VAR,
            "filter": (
                f"{FILTER_PERIOD_LOW:g}-{FILTER_PERIOD_HIGH:g} day "
                f"Lanczos band-pass, {FILTER_NWEIGHTS} weights"
            ),
            "training_years": f"{YEAR_START}-{YEAR_END}",
            "training_months": ",".join(map(str, TRAIN_MONTHS)),
            "eeof_lags_days": ",".join(map(str, LAGS.tolist())),
            "sign_method": SIGN_METHOD,
        },
    )

    Path(EEOF_FILE).unlink(missing_ok=True)

    ds_eeof.to_netcdf(
        EEOF_FILE,
        encoding={
            "eof": {
                "dtype": "float32",
                "zlib": True,
                "complevel": 2,
            },
            "eof_physical": {
                "dtype": "float32",
                "zlib": True,
                "complevel": 2,
            },
            "ev_ts": {
                "dtype": "float32",
                "zlib": True,
                "complevel": 2,
            },
            "time": time_encoding,
        },
    )

    # -------------------------------------------------------------------------
    # BSISO phase file
    # -------------------------------------------------------------------------
    print(f"Writing phase file: {PHASE_FILE}")

    ds_phase = xr.Dataset(
        data_vars={
            "bsiso_phase": (
                ("time",),
                phase,
                {
                    "long_name": "BSISO phase from standardized PC1-PC2 angle"
                },
            ),
            "bsiso_phase_nAmpGe1": (
                ("time",),
                phase_active,
                {
                    "long_name": (
                        "BSISO phase; 0 when standardized-PC amplitude <= 1"
                    )
                },
            ),
            "bsiso_amplitude": (
                ("time",),
                amplitude,
                {
                    "long_name": "BSISO amplitude sqrt(PC1^2 + PC2^2)",
                    "units": "1",
                },
            ),
            "bsiso_angle": (
                ("time",),
                theta,
                {
                    "long_name": "BSISO phase-space angle",
                    "units": "degree",
                },
            ),
            "pc1": (
                ("time",),
                pcs[0],
                {
                    "long_name": "standardized EEOF1 principal component",
                    "units": "1",
                },
            ),
            "pc2": (
                ("time",),
                pcs[1],
                {
                    "long_name": "standardized EEOF2 principal component",
                    "units": "1",
                },
            ),
        },
        coords={
            "time": pc_time,
        },
        attrs={
            "history": history,
            "source_olr": OLR_FILE,
            "filter": (
                f"{FILTER_PERIOD_LOW:g}-{FILTER_PERIOD_HIGH:g} day "
                f"Lanczos band-pass, {FILTER_NWEIGHTS} weights"
            ),
            "active_definition": f"amplitude > {AMP_THRESHOLD:g}",
        },
    )

    Path(PHASE_FILE).unlink(missing_ok=True)

    phase_fill = -2147483647

    ds_phase.to_netcdf(
        PHASE_FILE,
        encoding={
            "bsiso_phase": {
                "dtype": "int32",
                "_FillValue": phase_fill,
            },
            "bsiso_phase_nAmpGe1": {
                "dtype": "int32",
                "_FillValue": phase_fill,
            },
            "bsiso_amplitude": {
                "dtype": "float32",
            },
            "bsiso_angle": {
                "dtype": "float32",
            },
            "pc1": {
                "dtype": "float32",
            },
            "pc2": {
                "dtype": "float32",
            },
            "time": time_encoding,
        },
    )


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("BSISO EEOF CALCULATION FROM NOAA OLR ANOMALY ONLY")
    print("=" * 78)
    print(f"Input: {OLR_FILE}")

    # -------------------------------------------------------------------------
    # 1. Read NOAA OLR anomaly
    # -------------------------------------------------------------------------
    ds = xr.open_dataset(
        OLR_FILE,
        decode_times=True,
    )

    if OLR_VAR not in ds:
        raise KeyError(
            f"Variable {OLR_VAR!r} not found in {OLR_FILE}. "
            f"Available variables: {list(ds.data_vars)}"
        )

    olr = subset_lat(
        ds[OLR_VAR],
        LAT_MIN,
        LAT_MAX,
    ).transpose(
        "time",
        "lat",
        "lon"
    )

    check_daily_time(olr["time"])

    print("\nNOAA OLR anomaly:")
    print(f"  shape      = {olr.shape}")
    print(f"  start      = {olr.time.values[0]}")
    print(f"  end        = {olr.time.values[-1]}")
    print(f"  lat range  = {float(olr.lat.min()):g} to {float(olr.lat.max()):g}")
    print(f"  lon range  = {float(olr.lon.min()):g} to {float(olr.lon.max()):g}")

    # -------------------------------------------------------------------------
    # 2. Filter full NOAA record first
    # -------------------------------------------------------------------------
    filtered_full, weights = bandpass_filter_olr(olr)

    # -------------------------------------------------------------------------
    # 3. Restrict analysis record to 1980-2020
    # -------------------------------------------------------------------------
    filtered_analysis = filtered_full.sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31",
        )
    )

    filtered_array = np.asarray(
        filtered_analysis.values,
        dtype=DTYPE,
    )

    check_no_missing(
        filtered_array,
        f"filtered {YEAR_START}-{YEAR_END} OLR"
    )

    del filtered_array

    print("\nFiltered analysis period:")
    print(f"  start = {filtered_analysis.time.values[0]}")
    print(f"  end   = {filtered_analysis.time.values[-1]}")
    print(f"  days  = {filtered_analysis.time.size}")

    # For 1980-2020 this should be 14976 daily values.
    expected_days = (
        np.datetime64(f"{YEAR_END + 1}-01-01")
        - np.datetime64(f"{YEAR_START}-01-01")
    ).astype("timedelta64[D]").astype(int)

    if filtered_analysis.time.size != expected_days:
        raise ValueError(
            f"Expected {expected_days} daily values in "
            f"{YEAR_START}-{YEAR_END}, got {filtered_analysis.time.size}."
        )

    # -------------------------------------------------------------------------
    # 4. EEOF training data
    # -------------------------------------------------------------------------
    train_4d, nseason = prepare_training_data(
        filtered_analysis
    )

    # -------------------------------------------------------------------------
    # 5. Calculate covariance EEOFs
    # -------------------------------------------------------------------------
    eof, eigenvalues, pcvar = compute_eeof(
        train_4d,
        lags=LAGS,
        neof=NEOF,
    )

    # -------------------------------------------------------------------------
    # 6. Sign orientation
    # -------------------------------------------------------------------------
    lat = np.asarray(
        filtered_analysis["lat"].values,
        dtype=np.float64,
    )
    lon = np.asarray(
        filtered_analysis["lon"].values,
        dtype=np.float64,
    )

    print(f"\nApplying EEOF sign convention ({SIGN_METHOD})...")
    eof = enforce_sign_convention(
        eof,
        lat,
        lon,
    )

    # -------------------------------------------------------------------------
    # 7. Project entire 1980-2020 filtered record
    # -------------------------------------------------------------------------
    pcs, pc_time, raw_pc_std = project_full_record(
        filtered_analysis,
        eof,
        lags=LAGS,
    )

    # With 14976 input days and lag -10, result should be 14966.
    expected_pc_days = filtered_analysis.time.size + int(np.min(LAGS))

    if pc_time.size != expected_pc_days:
        raise RuntimeError(
            f"Expected {expected_pc_days} PC days, got {pc_time.size}."
        )

    print("\nPC diagnostics:")
    print(f"  raw PC std = {raw_pc_std}")
    print(
        "  corr(PC1,PC2) = "
        f"{np.corrcoef(pcs[0], pcs[1])[0,1]:.6f}"
    )

    # -------------------------------------------------------------------------
    # 8. BSISO phase and amplitude
    # -------------------------------------------------------------------------
    print("\nCalculating BSISO amplitude and phase...")

    amplitude, theta, phase, phase_active = calculate_bsiso_phase(
        pcs
    )

    print(
        f"  amplitude min/max = "
        f"{float(amplitude.min()):.3f} / "
        f"{float(amplitude.max()):.3f}"
    )

    unique, counts = np.unique(
        phase_active,
        return_counts=True,
    )

    print("\nBSISO active-phase counts:")
    for p, n in zip(unique, counts):
        label = "inactive" if p == 0 else f"phase {p}"
        print(f"  {label:8s}: {n}")

    # -------------------------------------------------------------------------
    # 9. Physical EEOF loading patterns for check plot
    # -------------------------------------------------------------------------
    eof_physical = physical_eeof_patterns(
        eof,
        eigenvalues,
    )

    # -------------------------------------------------------------------------
    # 10. Write NetCDF outputs
    # -------------------------------------------------------------------------
    write_outputs(
        filtered_analysis=filtered_analysis,
        eof=eof,
        eof_physical=eof_physical,
        eigenvalues=eigenvalues,
        pcvar=pcvar,
        pcs=pcs,
        pc_time=pc_time,
        amplitude=amplitude,
        theta=theta,
        phase=phase,
        phase_active=phase_active,
    )

    # -------------------------------------------------------------------------
    # 11. Check plot
    # -------------------------------------------------------------------------
    make_check_plot(
        eof_physical=eof_physical,
        pcvar=pcvar,
        lat=lat,
        lon=lon,
    )

    # -------------------------------------------------------------------------
    # Done
    # -------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("DONE")
    print("=" * 78)
    print(f"EEOF file : {EEOF_FILE}")
    print(f"Phase file: {PHASE_FILE}")
    print(f"Check plot: {CHECK_FIG}")

    if SAVE_FILTERED_OLR:
        print(f"Filtered  : {FILTERED_FILE}")

    print(
        f"\nExpected bsiso_phase.nc time dimension: "
        f"{expected_pc_days}"
    )

    ds.close()


if __name__ == "__main__":
    main()
