"""
Paired GTEx/TCGA performance plots for journal figures.

Input: one row per cancer and feature source, with mean and standard deviation.
Error bars denote standard deviations, not confidence intervals.
Requires matplotlib, pandas and numpy.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

################################################################################
def plotPairedPerformance(
        dataFrame: pd.DataFrame,
        title: str,
        meanCol: str = "aucMean",
        stdCol: str = "aucStd",
        cancerCol: str = "cancer",
        sourceCol: str = "source",
        cancerOrder: tuple[str, ...] = (
            "Colorectal", "Esophagus", "Liver", "Lung", "Stomach"
        ),
        ax=None,
        xLabel: str = "Mean AUC",
        xLimits: tuple[float, float] = (0.5, 1.0),
        panelLabel: str | None = None,
        showLegend: bool = True,
        yOffset: float = 0.0,
    ) -> tuple[plt.figure, plt.axes]:
    """
    Plot GTEx (blue circles) and TCGA (red squares) as mean ± SD.

    Input dataFrame columns default to: cancer, source, aucMean, aucStd.
    The source values must be 'GTEx' or 'TCGA'. Missing source pairs are allowed.
    A supplied ax lets the same function make panels A and B on one figure.
    
    Returns (fig, ax).
    """
    requiredCols = [cancerCol, sourceCol, meanCol, stdCol]
    missingCols = [col for col in requiredCols if col not in dataFrame.columns]
    if missingCols:
        raise ValueError(f"Missing DataFrame columns: {missingCols}")
    if dataFrame.duplicated([cancerCol, sourceCol]).any():
        raise ValueError("Expected at most one row for each cancer/source pair")

    if ax is None:
        fig, ax = plt.subplots(figsize=(90 / 25.4, 78 / 25.4))
        fig.subplots_adjust(left=0.28, right=0.97, top=0.84, bottom=0.17)
    else:
        fig = ax.figure

    colors = {"GTEx": "#0864b5", "TCGA": "#d62828"}
    markers = {"GTEx": "o", "TCGA": "s"}
    offsets = {"GTEx": -yOffset, "TCGA": yOffset}
    plotData = dataFrame.set_index([cancerCol, sourceCol])

    for rowNumber, cancer in enumerate(cancerOrder):
        for source in ("GTEx", "TCGA"):
            dataKey = (cancer, source)
            if dataKey not in plotData.index:
                continue
            meanValue = plotData.loc[dataKey, meanCol]
            stdValue = plotData.loc[dataKey, stdCol]
            if pd.isna(meanValue):
                continue
            meanValue = float(meanValue)
            if not 0 <= meanValue <= 1:
                raise ValueError(f"Metric must be in [0, 1]: {dataKey}")
            if pd.notna(stdValue) and float(stdValue) < 0:
                raise ValueError(f"Standard deviation cannot be negative: {dataKey}")

            errorValues = None
            if pd.notna(stdValue):
                stdValue = float(stdValue)
                errorValues = [
                    [min(stdValue, meanValue)],
                    [min(stdValue, 1.0 - meanValue)],
                ]

            ax.errorbar(
                meanValue,
                rowNumber + offsets[source],
                xerr=errorValues,
                fmt=markers[source],
                color=colors[source],
                markersize=4.5,
                markeredgewidth=0,
                elinewidth=1.0,
                capsize=2.4,
                capthick=1.0,
                linestyle="none",
                zorder=3 if source == "TCGA" else 2,
                clip_on=True,
            )

    ax.set_xlim(*xLimits)
    ax.set_ylim(len(cancerOrder) - 0.5, -0.5)
    ax.set_yticks(np.arange(len(cancerOrder)), labels=cancerOrder)
    ax.set_xticks(np.linspace(xLimits[0], xLimits[1], 6))
    ax.set_xlabel(xLabel, fontsize=7)
    ax.set_title(title, fontsize=7, fontweight="bold", pad=7)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.6, linestyle="--")
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=6.5, width=0.8, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    if panelLabel:
        ax.text(-0.13, 1.09, panelLabel, transform=ax.transAxes,
                fontsize=8, fontweight="bold", va="top", ha="left")

    if showLegend:
        legendHandles = [
            Line2D([], [], marker=markers[source], color=colors[source],
                   markersize=4.5, linewidth=1.0, label=f"{source}-derived features")
            for source in ("GTEx", "TCGA")
        ]
        ax.legend(handles=legendHandles, loc="lower left", fontsize=5.5,
                  frameon=True, facecolor="white", edgecolor="#cccccc",
                  borderpad=0.45, labelspacing=0.35, handlelength=1.5)

    return fig, ax


