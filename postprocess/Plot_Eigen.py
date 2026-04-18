# ====================================================
# Plot growth rate and phase speed of the eigenmodes
# ====================================================

# ====================================================
# Import packages
# ====================================================
import sys
import h5py
import numpy as np

from pathlib import Path
from typing import cast

from matplotlib import pyplot as plt

style_path = Path(__file__).resolve().parent / "SingleLine.mplstyle"
plt.style.use(["seaborn-v0_8-colorblind", str(style_path)])

# ====================================================
# Main function
# ====================================================
def main(output_path_str: str, fig_path_str: str) -> None:

    # Convert to Path objects
    output_path: Path = Path(output_path_str)
    fig_path: Path = Path(fig_path_str)

    # ----------------------------------------------------
    # Load data
    # ----------------------------------------------------
    with h5py.File(output_path / "EigenAnalysis.h5", "r") as f:

        # load data
        k: np.ndarray = cast(h5py.Dataset, f["k"])[...]
        growth : np.ndarray = cast(h5py.Dataset, f["GrowthRates"])[...]
        phase_speed : np.ndarray = cast(h5py.Dataset, f["PhaseSpeeds"])[...]

    # ----------------------------------------------------
    # Transform wavenumber to non-dimensional form
    # ----------------------------------------------------

    k_nondim = k * (40000.0 / (2*np.pi*4320.0))

    # ----------------------------------------------------
    # Visualization
    # ----------------------------------------------------
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))

    ax.plot(
        k_nondim,
        growth[:, 0],
        marker="o",
        markevery=1,
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=1,
        color="black",
        linestyle="-",
        linewidth=2,
    )
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 0.14)
    ax.set_xlabel("Non-dimensional Wavenumber")
    ax.set_ylabel("Growth Rate (1/day)")
    ax.tick_params(direction="in", length=6, width=1, top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y")
    plt.savefig(fig_path / "GrowthRate.png", dpi=300, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    ax.plot(
        k_nondim,
        phase_speed[:, 0],
        marker="o",
        markevery=1,
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=1,
        color="black",
        linestyle="-",
        linewidth=2,
    )
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 60)
    ax.set_xlabel("Non-dimensional Wavenumber")
    ax.set_ylabel("Phase Speed (m/s)")
    ax.tick_params(direction="in", length=6, width=1, top=True, right=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y")
    plt.savefig(fig_path / "PhaseSpeed.png", dpi=300, bbox_inches="tight")
    plt.close()

# ====================================================
# Run main function
# ====================================================
if __name__ == "__main__":
    output_path = str(sys.argv[1])
    fig_path = str(sys.argv[2])

    main(output_path, fig_path)