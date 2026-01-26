"""Simple generic utility functions for the module."""

import geopandas as gpd


def get_padded_bounds(
    gdf: gpd.GeoDataFrame, *, pad_frac: float = 0.01
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Get padded boundaries of a spatial dataset.

    Args:
        gdf (gpd.GeoDataFrame): dataset to get boundaries from.
        pad_frac (float, optional): ratio to pad around the box. Defaults to 0.01.

    Returns:
        tuple[tuple[float,float], tuple[float,float]]: min, max for X, Y.
    """
    xmin, ymin, xmax, ymax = map(float, gdf.total_bounds)
    dx, dy = (xmax - xmin), (ymax - ymin)

    # Avoid degenerate cases
    if dx == 0:
        dx = 1.0
    if dy == 0:
        dy = 1.0

    pad_x = dx * pad_frac
    pad_y = dy * pad_frac

    return (xmin - pad_x, xmax + pad_x), (ymin - pad_y, ymax + pad_y)
