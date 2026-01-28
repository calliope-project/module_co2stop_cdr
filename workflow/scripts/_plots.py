"""Plotting helpers."""

import geopandas as gpd
import pandas as pd
from cmap import Colormap
from matplotlib import pyplot as plt


def plot_aggregate(
    shapes: gpd.GeoDataFrame, aggregated: pd.DataFrame, cmap="cmasher:sepia_r"
):
    """Plot the aggregation result."""
    fig, ax = plt.subplots(layout="compressed")
    combined = shapes.merge(aggregated, how="inner", on="shape_id")

    shapes.boundary.plot(lw=0.5, color="grey", ax=ax)
    combined.plot(
        "max_sequestered_mtco2",
        cmap=Colormap(cmap).to_mpl(),
        ax=ax,
        legend=True,
        legend_kwds={"label": "$MtCO_2$"},
    )
    ax.set_axis_off()
    return fig, ax
