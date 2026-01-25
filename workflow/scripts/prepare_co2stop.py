"""Prepare the CO2Stop data so it fits our schemas.

Dataset-wide imputations happen here.

This code was adapted from PyPSA-Eur (https://github.com/pypsa/pypsa-eur)
Copyright (c) 2017-2024 The PyPSA-Eur Authors
Licensed under the MIT License
Commit: 630d37dd061fda6fff93c3a2c458dcdfdc9dcedd
File: scripts/build_co2_sequestration_potentials.py
"""

import sys
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry as sg
from shapely.ops import unary_union

if TYPE_CHECKING:
    snakemake: Any

AQUIFERS = {
    "primary": {
        "conservative_mtco2": "EST_STORECAP_MIN",
        "neutral_mtco2": "EST_STORECAP_MEAN",
        "optimistic_mtco2": "EST_STORECAP_MAX",
    },
    "fallback": {
        "conservative_mtco2": "STORE_CAP_MIN",
        "neutral_mtco2": "STORE_CAP_MEAN",
        "optimistic_mtco2": "STORE_CAP_MAX",
    },
}
GAS = {
    "primary": {
        "conservative_mtco2": "MIN_EST_STORE_CAP_GAS",
        "neutral_mtco2": "MEAN_EST_STORE_CAP_GAS",
        "optimistic_mtco2": "MAX_EST_STORE_CAP_GAS",
    },
    "fallback": {
        "conservative_mtco2": "MIN_CALC_STORE_CAP_GAS",
        "neutral_mtco2": "MEAN_CALC_STORE_CAP_GAS",
        "optimistic_mtco2": "MAX_CALC_STORE_CAP_GAS",
    },
}
OIL = {
    "primary": {
        "conservative_mtco2": "MIN_EST_STORE_CAP_OIL",
        "neutral_mtco2": "MEAN_EST_STORE_CAP_OIL",
        "optimistic_mtco2": "MAX_EST_STORE_CAP_OIL",
    },
    "fallback": {
        "conservative_mtco2": "MIN_CALC_STORE_CAP_OIL",
        "neutral_mtco2": "MEAN_CALC_STORE_CAP_OIL",
        "optimistic_mtco2": "MAX_CALC_STORE_CAP_OIL",
    },
}


def _mask_surface_issues(df: pd.DataFrame) -> pd.Series:
    """Only empty or 'None' cases are considered safe."""
    sufrace_issues = df["SURF_ISSUES"].replace("None", np.nan)
    return sufrace_issues.isna()


def _mask_subsurface_interferance(df: pd.DataFrame) -> pd.Series:
    """Only empty or 'No' cases are considered safe."""
    subsurface_interf = df["SUBSURF_INTERF"].replace("No", np.nan)
    return subsurface_interf.isna()


def _mask_unassessed(df: pd.DataFrame, cols: str | list[str]):
    """Detect unassessed cases.

    Based on the methods described in Tumara et al 2024 (p.12).
    https://dx.doi.org/10.2760/582433
    """
    if isinstance(cols, str):
        cols = [cols]

    empty_masks = []
    for c in cols:
        s = df[c]
        if pd.api.types.is_string_dtype(s) or s.dtype == object:
            empty = s.isna() | (s.astype(str).str.strip() == "")
        else:
            empty = s.isna()
        empty_masks.append(empty)

    all_empty = pd.concat(empty_masks, axis=1).all(axis=1)
    return ~all_empty


def estimate_storage_scenarios(
    df: pd.DataFrame,
    primary_cols: dict[str, str],
    fallback_cols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Get minimum, mean and maximum CO2 capacity per storage unit.

    - Data for minimum, mean, and maximum CO2 capacity is taken from the primary column.
    - If no primary data is present, the fallback column will be used.
    - If at least one value is present (i.e., only mean is given), other categories will
    be filled with it.
    Smaller values are given priority over larger ones.
    E.g.: if min and max are present, mean will be filled with min first.
    - Corrections are applied to ensure monotonic behaviour per row
    (i.e., `conservative <= neutral <= optimistic`).
    """
    if len(primary_cols) != 3:
        raise ValueError(
            "`primary_cols` must have length 3, ordered as min < mean < max."
        )

    out = pd.DataFrame(index=df.index)

    for name, col in primary_cols.items():
        s = df[col].replace(0, np.nan)
        if fallback_cols is not None:
            s = s.fillna(df[fallback_cols[name]].replace(0, np.nan))
        out[name] = s

    # Bidirectional propagation within each row
    out = out.ffill(axis=1).bfill(axis=1)

    lo, mid, hi = list(primary_cols.keys())

    # Enforce lo <= mid <= hi (preserving NaNs)
    m = out[lo].notna() & out[mid].notna() & (out[lo] > out[mid])
    out[lo] = out[lo].where(~m, out[mid])

    m = out[hi].notna() & out[mid].notna() & (out[hi] < out[mid])
    out[hi] = out[hi].where(~m, out[mid])

    return out


def main() -> None:
    """Main snakemake process."""
    pass


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    main()
