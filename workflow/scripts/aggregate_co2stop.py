"""Aggregate CO2Stop to provided shapes."""

import sys
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

import _schemas
import geopandas as gpd
import pandas as pd
from cmap import Colormap
from matplotlib import pyplot as plt
from pyproj import CRS

if TYPE_CHECKING:
    snakemake: Any


def build_scenario_gdf(
    storage_unit_file: str,
    trap_files: list[str],
    *,
    datasets: Iterable[str],
    groups: Iterable[str],
    scenario: str,
) -> gpd.GeoDataFrame:
    """Load and combine requested datasets into a scenario.

    Args:
        storage_unit_file (str): path to storage unit dataset.
        trap_files (list[str]): paths to all trap datasets.
        datasets (Iterable[str]): list of datasets to combine.
        groups (Iterable[str]): list of CO2 sink groups to use.
        scenario (str): scenario name. One of: low, medium, high.

    Raises:
        ValueError: the dataset / group combination lead to an empty scenario.

    Returns:
        gpd.GeoDataFrame: resulting scenario combination.
    """
    cols = [scenario, "dataset", "storage_group", "geometry"]

    traps = pd.concat([gpd.read_parquet(p) for p in trap_files], ignore_index=True)
    storage_units = gpd.read_parquet(storage_unit_file)

    # Filter to requested case
    traps = traps.loc[
        traps["storage_group"].isin(groups) & traps["dataset"].isin(datasets)
    ]
    storage_units = storage_units.loc[
        storage_units["storage_group"].isin(groups)
        & storage_units["dataset"].isin(datasets)
    ]

    # Remove traps already represented by storage_units
    traps = traps.loc[~traps["storage_unit_id"].isin(storage_units["storage_unit_id"])]

    scenario_gdf = gpd.GeoDataFrame(
        pd.concat([storage_units[cols], traps[cols]], ignore_index=True),
        geometry="geometry",
        crs=storage_units.crs,
    )
    if scenario_gdf.empty:
        raise ValueError(
            f"Request '{datasets}:{groups}:{scenario}' produced empty scenario."
        )

    scenario_gdf["scenario_id"] = scenario_gdf.index
    return scenario_gdf


def aggregate_scenario_into_shapes(
    shapes: gpd.GeoDataFrame, scenario_gdf: gpd.GeoDataFrame, *, scenario_col: str
) -> pd.DataFrame:
    """Overlay scenario polygons with target shapes and MtCO2 per area."""
    if (not shapes.crs.is_projected) or (not shapes.crs.equals(scenario_gdf.crs)):
        raise ValueError("Provided files must share a projected CRS.")

    overlay = gpd.overlay(
        shapes[["shape_id", "geometry"]],
        scenario_gdf,
        how="intersection",
        keep_geom_type=False,
    )

    overlay["piece_area"] = overlay.area
    overlay = overlay.loc[overlay["piece_area"] > 0].copy()

    scenario_area = scenario_gdf.area
    overlay["source_area"] = overlay["scenario_id"].map(scenario_area)
    overlay = overlay.loc[overlay["source_area"] > 0].copy()

    overlay["max_sequestered_mtco2"] = (
        overlay[scenario_col] * overlay["piece_area"] / overlay["source_area"]
    )

    return overlay.groupby("shape_id", as_index=False)["max_sequestered_mtco2"].sum()


def plot(shapes: gpd.GeoDataFrame, aggregated: pd.DataFrame, cmap="cmasher:amber_r"):
    """Plot the aggregation result."""
    fig, ax = plt.subplots(layout="constrained")
    combined = shapes.merge(aggregated, how="inner", on="shape_id")

    shapes.boundary.plot(lw=0.5, color="black", ax=ax)
    combined.plot(
        "max_sequestered_mtco2", legend=True, cmap=Colormap(cmap).to_mpl(), ax=ax
    )
    ax.set_title("Aggregated $MtCO_2$")
    return fig, ax


def main() -> None:
    """Main snakemake process."""
    proj_crs = snakemake.params.proj_crs
    if not CRS.from_user_input(proj_crs).is_projected:
        raise ValueError(f"Expected projected CRS, got {proj_crs!r}.")

    datasets = snakemake.params.included_datasets
    groups = snakemake.params.included_groups
    scenario = f"{snakemake.wildcards.scenario}_mtco2"

    shapes = _schemas.ShapeSchema.validate(gpd.read_parquet(snakemake.input.shapes))
    shapes = shapes.to_crs(proj_crs)

    scenario_gdf = build_scenario_gdf(
        storage_unit_file=snakemake.input.storage_units,
        trap_files=snakemake.input.all_traps,
        datasets=datasets,
        groups=groups,
        scenario=scenario,
    )

    aggregated = aggregate_scenario_into_shapes(
        shapes=shapes, scenario_gdf=scenario_gdf.to_crs(proj_crs), scenario_col=scenario
    )
    aggregated = _schemas.AggregatedSchema.validate(aggregated)
    aggregated.to_parquet(snakemake.output.aggregated)

    fig, _ = plot(shapes, aggregated)
    fig.savefig(snakemake.output.plot, dpi=300)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    main()
