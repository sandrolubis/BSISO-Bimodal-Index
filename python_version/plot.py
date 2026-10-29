import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ==========================================================
# USER SETTINGS
# ==========================================================

YEAR1 = 2011
YEAR2 = 2012

file_my  = "../bsiso_phase.nc"
file_ref = "../BSISO_25-90bpfil.rt_pc.txt"

fout = f"compare_BSISO_JJA_{YEAR1}_{YEAR2}.png"


# ==========================================================
# READ YOUR BSISO INDEX
# ==========================================================

ds = xr.open_dataset(file_my)

dates_my = pd.DatetimeIndex(
    pd.to_datetime(ds.time.values)
).normalize()

my_index = pd.DataFrame(
    {
        "PC1":   ds.PC1.values,
        "PC2":   ds.PC2.values,
        "amp":   ds.amp.values,
        "phase": ds.bsiso_phase.values,
    },
    index=dates_my,
)

ds.close()


# ==========================================================
# READ REFERENCE BSISO INDEX
# ==========================================================

columns = [
    "year", "month", "day",
    "PCx", "PCy", "phase",
    "amp_norm", "amp_raw"
]

ref = pd.read_csv(
    file_ref,
    sep=r"\s+",
    skiprows=1,
    header=None,
    names=columns,
)

ref["date"] = pd.to_datetime(
    ref[["year", "month", "day"]]
)

ref = ref.set_index("date").sort_index()


# ==========================================================
# VARIABLES TO COMPARE
# ==========================================================

variables = [
    ("PC1",   "PCx",      "PC1"),
    ("PC2",   "PCy",      "PC2"),
    ("amp",   "amp_norm", "Amplitude"),
    ("phase", "phase",    "BSISO Phase"),
]


# ==========================================================
# CREATE FOUR-PANEL FIGURE
# ==========================================================

fig, axes = plt.subplots(
    4, 1,
    figsize=(13, 10),
    sharex=True
)

years = [YEAR1, YEAR2]

for iyear, year in enumerate(years):

    # JJA dates
    dates = pd.date_range(
        f"{year}-06-01",
        f"{year}-08-31",
        freq="D"
    )

    # 92 days per summer, plus a gap between years
    x = np.arange(92) + iyear * 93

    # Match data by calendar date
    mine  = my_index.reindex(dates)
    other = ref.reindex(dates)

    print(
        f"JJA {year}: "
        f"My index = {mine['PC1'].notna().sum()} days; "
        f"Reference = {other['PCx'].notna().sum()} days"
    )

    # Plot four variables
    for ax, (my_var, ref_var, ylabel) in zip(axes, variables):

        # Your BSISO: red
        ax.plot(
            x,
            mine[my_var],
            color="red",
            linewidth=1.7,
            label="My BSISO" if iyear == 0 else None
        )

        # Reference BSISO: black
        ax.plot(
            x,
            other[ref_var],
            color="black",
            linewidth=1.7,
            label="Reference BSISO" if iyear == 0 else None
        )

        ax.set_ylabel(ylabel, fontsize=11)


# ==========================================================
# FIGURE FORMATTING
# ==========================================================

ticks = [0, 30, 61, 93, 123, 154]

tick_labels = [
    f"Jun {YEAR1}", "Jul", "Aug",
    f"Jun {YEAR2}", "Jul", "Aug"
]

for ax in axes:

    # Separator between summers
    ax.axvline(
        92,
        color="gray",
        linestyle="--",
        linewidth=0.9
    )

    ax.grid(
        True,
        color="0.88",
        linewidth=0.6
    )

    ax.set_xlim(-1, 185)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# Legend
axes[0].legend(
    loc="upper right",
    frameon=False,
    fontsize=10,
    ncol=2
)

# Amplitude
axes[2].set_ylim(bottom=0)

# Phase
axes[3].set_ylim(0.5, 8.5)
axes[3].set_yticks(np.arange(1, 9))

# Time labels
axes[3].set_xticks(ticks)
axes[3].set_xticklabels(tick_labels)
axes[3].set_xlabel("Month")

# Main title
fig.suptitle(
    f"BSISO Comparison: JJA {YEAR1} and {YEAR2}",
    fontsize=15
)

fig.tight_layout(rect=[0, 0, 1, 0.97])

# Save
plt.savefig(
    fout,
    dpi=250,
    bbox_inches="tight"
)

plt.show()

print(f"Saved: {fout}")
