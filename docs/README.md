# Module CO2Stop - Carbon Dioxide Removal

Welcome to the documentation of the `module_co2stop_cdr` data module!
This module performs geospatial aggregation of $CO_2$ sequestration potentials for Europe using the [CO2Stop dataset](https://energy.ec.europa.eu/publications/assessment-co2-storage-potential-europe-co2stop_en).

Please consult the [specification guidelines](./specification.md) and the [`clio` documentation](https://clio.readthedocs.io/) for more information.

## Overview

![rulegraph](./rulegraph.png)

The analysis of the module is structured as follows:

1. The CO2Stop dataset is downloaded and cleaned up following methods described in the [JRC - CO2 transport report](https://publications.jrc.ec.europa.eu/repository/handle/JRC136709).
Storage units and traps are removed according to the following.
    - Cases that were not assessed or for which data was undisclosed.
    - Ambiguous duplicates (these are two small traps located in the north sea with repeated IDs and capacities).
    - Optionally, details from the dataset are used to remove the following, if configured:
        - Qualitative cases marked as having surface or subsurface issues, and artificially created polygons.
        - Quantitative values (e.g., porosity, depth, ...).
1. To avoid double-counting, traps within the remaining storage units are removed as their capacity is already included in the storage unit total.
Please consult the [CO2Stop Final report](https://energy.ec.europa.eu/publications/assessment-co2-storage-potential-europe-co2stop_en) (section 2.3.1) for details.
1. Three scenarios (`low`, `medium`, `high`) are created for each sequestration type (`aquifer`, `gas`, `oil`) for the remaining CO2Stop data.
User-configured lower and upper bounds are applied per-polygon at this stage.
See `bounds_mtco2: co2stop_polygons` in the configuration schema for more information.
![scenarios](./aquifer_scenarios.png)
1. The resulting sequestration potential is aggregated per scenario into user provided shapes.
These shapes should follow the schema provided by the [geo-boundaries module](https://github.com/calliope-project/module_geo_boundaries/tree/v0.1.6).

>[!WARNING]
>Estimates from the CO2Stop dataset are biased by disclosure (or lack thereof), and the filtering settings used.
>Some countries are affected more than others, with Germany having particularly poor disclosure.
>
>Similarly, CO2Stop suffers from poor data handling practices that make unavailable data and 'true' zero values indistinguishable from each other, amplifying the uneven assignation of sequestration. For example: setting `porosity_ratio: 0.1` will completely remove France in most cases.
>
>We provide automated figures and logging (in `logs/storage_units/` and `logs/traps/`) so users can evaluate how their settings affect polygon selection.
>Below is an example for storage unit aquifers where only undisclosed and artificial polygons have been removed. This can be seen as a _MINIMUM_ amount of removals.
>![filters](./aquifer_kept.png)

## Configuration

See [the configuration README](./../config/README.md).

## Outputs

See the [interface file](./../INTERFACE.yaml).

## Data sources

- Poulsen, N., Holloway, S., Neele, F., Smith, N.A., Kirk, K., 2012. CO2StoP Executive Summary (No. ENER/C1/154-2011-SI2.611598). GEOLOGICAL SURVEY OF DENMARK AND GREENLAND. <https://energy.ec.europa.eu/publications/assessment-co2-storage-potential-europe-co2stop_en>.

## References

This module relies on code and methods from the following sources.

- **Shape schema definition:**
Ruiz Manuel, I. clio - module_geo_boundaries [Computer software]. <https://github.com/calliope-project/module_geo_boundaries/>.
- **Filtering minimum defaults:**
Van Den Broek, M., Brederode, E., Ramírez, A., Kramers, L., Van Der Kuip, M., Wildenborg, T., Turkenburg, W., Faaij, A., 2010. Designing a cost-effective CO2 storage infrastructure using a GIS based linear optimization energy model. Environmental Modelling & Software 25, 1754–1768. <https://doi.org/10.1016/j.envsoft.2010.06.015>.
