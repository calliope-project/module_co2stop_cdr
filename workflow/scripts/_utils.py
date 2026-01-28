"""Simple generic utility functions for the module."""

from dataclasses import dataclass

import geopandas as gpd


@dataclass
class CDRGroup:
    """Configuration for a given type of storage case."""

    primary: dict[str, str]
    """Main data points (min, mean, max)."""
    fallback: dict[str, str]
    """Fallback data points (min, mean, max)."""
    methods: list[str]
    """Columns specifying calculation method."""


AQUIFER = CDRGroup(
    primary={
        "low_mtco2": "EST_STORECAP_MIN",
        "medium_mtco2": "EST_STORECAP_MEAN",
        "high_mtco2": "EST_STORECAP_MAX",
    },
    fallback={
        "low_mtco2": "STORE_CAP_MIN",
        "medium_mtco2": "STORE_CAP_MEAN",
        "high_mtco2": "STORE_CAP_MAX",
    },
    methods=["CAP_EST_METHOD", "CAP_CAL_METHOD"],
)
GAS = CDRGroup(
    primary={
        "low_mtco2": "MIN_EST_STORE_CAP_GAS",
        "medium_mtco2": "MEAN_EST_STORE_CAP_GAS",
        "high_mtco2": "MAX_EST_STORE_CAP_GAS",
    },
    fallback={
        "low_mtco2": "MIN_CALC_STORE_CAP_GAS",
        "medium_mtco2": "MEAN_CALC_STORE_CAP_GAS",
        "high_mtco2": "MAX_CALC_STORE_CAP_GAS",
    },
    methods=["EST_METHOD_GAS", "CALC_METHOD_GAS"],
)
OIL = CDRGroup(
    primary={
        "low_mtco2": "MIN_EST_STORE_CAP_OIL",
        "medium_mtco2": "MEAN_EST_STORE_CAP_OIL",
        "high_mtco2": "MAX_EST_STORE_CAP_OIL",
    },
    fallback={
        "low_mtco2": "MIN_CALC_STORE_CAP_OIL",
        "medium_mtco2": "MEAN_CALC_STORE_CAP_OIL",
        "high_mtco2": "MAX_CALC_STORE_CAP_OIL",
    },
    methods=["EST_METHOD_OIL", "CALC_METHOD_OIL"],
)
# Handle snakemake wildcards
CDR_GROUP: dict[str, CDRGroup] = {"aquifer": AQUIFER, "gas": GAS, "oil": OIL}


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
