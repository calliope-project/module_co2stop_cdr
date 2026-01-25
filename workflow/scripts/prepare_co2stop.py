"""Prepare the CO2Stop data so it fits our schemas.

Dataset-wide imputations happen here.
"""

import re
import sys
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd

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
    "methods": ["CAP_EST_METHOD", "CAP_CAL_METHOD"],
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
    "methods": ["EST_METHOD_GAS", "CALC_METHOD_GAS"],
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
    "methods": ["EST_METHOD_OIL", "CALC_METHOD_OIL"],
}


def get_surface_issues(df: pd.DataFrame) -> pd.Series:
    """Detect surface issues per row.

    Columns used:
    - SURF_ISSUES: Only empty or 'None' cases are considered safe.
    - REMARKS: additional remarks by CO2Stop authors (e.g., land ownership issues).

    Args:
        df (pd.DataFrame): dataframe with CO2Stop data.

    Returns:
        pd.Series: True if issue is present. False otherwise.
    """
    issues = (
        df["SURF_ISSUES"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"none": np.nan, "": np.nan})
    )
    unsafe = issues.notna()

    problems = ["surface issue ="]
    pattern = "|".join(re.escape(i) for i in problems)

    flagged_in_remarks = (
        df["REMARKS_DATA"].fillna("").astype(str).str.lower().str.contains(pattern, regex=True)
    )
    unsafe = unsafe | flagged_in_remarks

    return unsafe


def get_subsurface_interference(df: pd.DataFrame) -> pd.Series:
    """Detect subsurface issues per row.

    Columns used:
    - SUBSURF_INTERF: Only empty or 'No' cases are considered safe.
    - REMARKS: additional remarks by CO2Stop authors (e.g., groundwater source).

    Args:
        df (pd.DataFrame): dataframe with CO2Stop data.

    Returns:
        pd.Series: True if issue is present. False otherwise.
    """
    subsurface_interf = (
        df["SUBSURF_INTERF"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"no": np.nan, "": np.nan})
    )
    unsafe = subsurface_interf.notna()

    problems = ["subsurface issue =", "geothermal", "groundwater", "potable water"]
    pattern = "|".join(re.escape(i) for i in problems)
    flagged_in_remarks = (
        df["REMARKS_DATA"].fillna("").astype(str).str.lower().str.contains(pattern, regex=True)
    )
    unsafe |= flagged_in_remarks

    return unsafe



def get_fake_polygons(df: pd.DataFrame) -> pd.Series:
    """Detect cases where the polygon is artificial.

    Uses:
    - REMARKS (from polygons)
    - REMARKS_DATA (from CSV merge, if present)
    """
    checks = [
        ("REMARKS", [
            "polygon does not represent",
            "polygon in no way represents",
            "arbitrary storage unit polygon",
        ]),
        ("REMARKS_DATA", [
            "fictive saline aquifer",
            "polygon not available",
            "aproximated polygon",
            "polygon aproximated",
        ]),
    ]

    fake = pd.Series(False, index=df.index)

    for col, problems in checks:
        pattern = "|".join(re.escape(p) for p in problems)
        flagged = (
            df[col]
            .fillna("")
            .astype(str)
            .str.contains(pattern, case=False, regex=True)
        )
        fake |= flagged

    return fake


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

    all_empty = pd.concat(empty_masks, axis="columns").all(axis="columns")
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
    out = out.ffill(axis="columns").bfill(axis="columns")

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
