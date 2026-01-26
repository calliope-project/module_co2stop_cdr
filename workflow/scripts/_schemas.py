"""Schemas for key files."""

from _utils import CDR_GROUP
from pandera.pandas import DataFrameModel, Field, check
from pandera.typing.geopandas import GeoSeries
from pandera.typing.pandas import Index, Series


class ShapeSchema(DataFrameModel):
    class Config:
        coerce = True
        strict = "filter"

    shape_id: Series[str] = Field(unique=True)
    "Unique ID for this shape."
    geometry: GeoSeries
    "Shape polygon."

    @check("geometry", element_wise=True)
    def geom_not_empty(cls, geom):
        return (geom is not None) and (not geom.is_empty) and geom.is_valid


class CO2StopSchema(DataFrameModel):
    """Generic CO2Stop values."""

    class Config:
        coerce = True
        strict = True

    index: Index[int] = Field(unique=True)

    low_mtco2: Series[float]
    "Low storage scenario."
    medium_mtco2: Series[float]
    "Medium storage scenario."
    high_mtco2: Series[float]
    "High storage scenario."
    geometry: GeoSeries
    "Shape polygon."

    @check("geometry", element_wise=True)
    def geom_not_empty(cls, geom):
        return (geom is not None) and (not geom.is_empty) and geom.is_valid


class StorageUnitsSchema(CO2StopSchema):
    """Storage unit specifics."""

    storage_unit_id: Series[str] = Field(unique=True)
    "Unique storage unit identifier."
    dataset: Series[str] = Field(eq="storage_units")
    "Parent dataset."
    storage_group: Series[str] = Field(eq="aquifer")
    "Storage group."


class TrapsSchema(CO2StopSchema):
    """Trap specifics."""

    trap_id: Series[str] | None = Field(unique=True)
    "Unique trap identifier."
    storage_unit_id: Series[str]
    "Storage unit identifier (may repeat)."
    dataset: Series[str] = Field(eq="traps")
    "Parent dataset."
    storage_group: Series[str] = Field(isin=list(CDR_GROUP.keys()))
    "Storage group."


class AggregatedSchema(DataFrameModel):
    """Validates the final output dataframes."""

    class Config:
        coerce = True
        strict = True

    shape_id: Series[str] = Field(unique=True)
    "Unique ID for this shape."
    max_sequestered_mtco2: Series[float]
