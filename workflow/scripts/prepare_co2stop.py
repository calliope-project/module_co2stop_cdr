"""Prepare the CO2Stop data so it fits our schemas.

Dataset-wide imputations happen here.
"""

import re
import sys
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

import _schemas
import geopandas as gpd
import numpy as np
import pandas as pd
from _utils import CDR_GROUP, StorageGroup, get_padded_bounds
from cmap import Colormap
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def _surface_issues(df: pd.DataFrame) -> pd.Series:
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
        df["REMARKS_DATA"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(pattern, regex=True)
    )
    unsafe = unsafe | flagged_in_remarks

    return unsafe


def _subsurface_interference_issues(df: pd.DataFrame) -> pd.Series:
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
        df["REMARKS_DATA"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(pattern, regex=True)
    )
    unsafe |= flagged_in_remarks

    return unsafe


def _artificial_polygon_issues(df: pd.DataFrame) -> pd.Series:
    """Detect cases where the polygon is artificial.

    Uses:
    - REMARKS (from polygons)
    - REMARKS_DATA (from CSV merge, if present)
    """
    checks = [
        (
            "REMARKS",
            [
                "polygon does not represent",
                "polygon in no way represents",
                "arbitrary storage unit polygon",
            ],
        ),
        (
            "REMARKS_DATA",
            [
                "fictive saline aquifer",
                "polygon not available",
                "aproximated polygon",
                "polygon aproximated",
            ],
        ),
    ]

    fake = pd.Series(False, index=df.index)

    for col, problems in checks:
        pattern = "|".join(re.escape(p) for p in problems)
        flagged = (
            df[col].fillna("").astype(str).str.contains(pattern, case=False, regex=True)
        )
        fake |= flagged

    return fake


def _ambiguous_duplicate_issues(df: pd.DataFrame, id_col: str) -> pd.Series:
    """Detect ambiguous cases with repeated IDs.

    All duplicates are eliminated because attribution is uncertain.
    """
    return df[id_col].duplicated(keep=False)


def get_assessed(df: pd.DataFrame, cols: Iterable[str]):
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

    all_empty = pd.concat(empty_masks, axis="columns").any(axis="columns")
    return ~all_empty


def estimate_storage_scenarios(
    df: pd.DataFrame, storage_group: StorageGroup
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
    if len(storage_group.primary) != 3:
        raise ValueError(
            "`primary_cols` must have length 3, ordered as min < mean < max."
        )
    out = pd.DataFrame(index=df.index)

    for name, col in storage_group.primary.items():
        s = df[col].replace(0, np.nan)
        s = s.fillna(df[storage_group.fallback[name]].replace(0, np.nan))
        out[name] = s

    # Bidirectional propagation within each row
    out = out.ffill(axis="columns").bfill(axis="columns")

    # Enforce lo <= mid <= hi (preserving NaNs)
    lo, mid, hi = list(storage_group.primary.keys())
    m = out[lo].notna() & out[mid].notna() & (out[lo] > out[mid])
    out[lo] = out[lo].where(~m, out[mid])

    m = out[hi].notna() & out[mid].notna() & (out[hi] < out[mid])
    out[hi] = out[hi].where(~m, out[mid])

    return out


def harmonise_stopco2_dataset(
    map_file: str, data_file: str, id_col: str, crs: str | int, suffix: str = "_DATA"
) -> gpd.GeoDataFrame:
    """Open and combine paired CO2Stop datasets."""
    gdf = gpd.read_file(map_file).rename({"id": id_col}, axis="columns").to_crs(crs)
    gdf.geometry = gdf.geometry.force_2d().make_valid()
    return gdf.merge(
        pd.read_csv(data_file), how="inner", on=id_col, suffixes=("", suffix)
    )


def plot_polygon_issues(
    countries: gpd.GeoDataFrame,
    data: gpd.GeoDataFrame,
    *,
    issues_col: str = "issues",
    cmap: str = "tol:high_contrast_alt",
) -> tuple[Figure, Axes]:
    """Show a combination of all dropped cases."""
    fig, ax = plt.subplots(layout="constrained")
    countries.plot(color="grey", alpha=0.5, ax=ax)
    countries.boundary.plot(color="black", lw=0.5, ax=ax)
    data.plot(issues_col, legend=True, ax=ax, cmap=Colormap(cmap).to_mpl())

    x_lim, y_lim = get_padded_bounds(data, pad_frac=0.02)
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_axis_off()
    return fig, ax


def plot_scenarios(data: pd.DataFrame) -> tuple[Figure, list[Axes]]:
    """Show a quick comparison between each scenario."""
    axes: list[Axes]
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), layout="constrained")

    scen_names = data.columns.str.split("_", n=1).str[0].tolist()
    axes[0].bar(x=scen_names, height=data.sum().values)
    axes[0].set_title("Aggregate")
    tmp = data.T.copy()
    tmp.index = scen_names
    tmp.plot(ax=axes[1], legend=False, color="grey", alpha=0.5, marker="o")
    axes[1].set_title("Per polygon")
    for ax in axes:
        ax.set_ylabel("$MtCO_2$")
        ax.tick_params(axis="x", which="both", length=0)
    return fig, axes


def main() -> None:
    """Main snakemake process."""
    geo_crs = snakemake.params.geo_crs

    if not CRS.from_user_input(geo_crs).is_geographic:
        raise ValueError(f"Expected geographic CRS, got {geo_crs!r}.")

    dataset_name = snakemake.wildcards.dataset
    storage_group = snakemake.wildcards.cdr_group

    match dataset_name:
        case "storage_units":
            data_id = "STORAGE_UNIT_ID"
            id_columns = {data_id: "storage_unit_id"}
            validation_method = _schemas.StorageUnitsSchema.validate
        case "traps":
            data_id = "TRAP_ID"
            id_columns = {data_id: "trap_id", "STORAGE_UNIT_ID": "storage_unit_id"}
            validation_method = _schemas.TrapsSchema.validate
        case _:
            raise ValueError(f"Invalid dataset requested: {dataset_name!r}.")

    dataset = harmonise_stopco2_dataset(
        snakemake.input.polygons, snakemake.input.table, data_id, geo_crs
    )

    # Try to catch problematic cases
    dataset["issues"] = _surface_issues(dataset)
    dataset["issues"] |= _subsurface_interference_issues(dataset)
    dataset["issues"] |= _artificial_polygon_issues(dataset)
    dataset["issues"] |= _ambiguous_duplicate_issues(dataset, data_id)
    # Plot cases with identified problems.
    countries = gpd.read_file(snakemake.input.countries).to_crs(geo_crs)
    fig, ax = plot_polygon_issues(countries, dataset)
    ax.set_title(f"'{dataset_name}:{storage_group}': polygons with identified issues.")
    fig.savefig(snakemake.output.plot_issues, dpi=300)
    # Remove bad apples
    dataset = dataset[~dataset["issues"]]

    # Estimate storage capacity
    capacity_scenarios = estimate_storage_scenarios(dataset, CDR_GROUP[storage_group])
    # Plot scenarios
    fig, _ = plot_scenarios(capacity_scenarios)
    fig.suptitle(f"'{dataset_name}:{storage_group}': scenario comparison")
    fig.savefig(snakemake.output.plot_scenarios, dpi=300)

    # Construct the final dataset
    result = dataset[list(id_columns.keys()) + ["geometry"]].join(capacity_scenarios)
    result = result.dropna(
        subset=capacity_scenarios.columns, how="all", ignore_index=True
    ).rename(id_columns, axis="columns")
    result["dataset"] = dataset_name
    result["storage_group"] = storage_group
    result = validation_method(result)
    result.to_parquet(snakemake.output.mtco2)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    main()
