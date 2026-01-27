"""Rules for database preparation."""

CDR_GROUP = ["aquifer", "gas", "oil"]


rule prepare_co2stop_storage_units:
    message:
        "Harmonising CO2Stop {wildcards.dataset}:{wildcards.cdr_group}."
    params:
        bounds=lambda wc: config["imputation"]["bounds_mtco2"]["co2stop_polygons"],
        filters=config["imputation"]["filter"],
        geo_crs=config["crs"]["geographic"],
    input:
        table="resources/automatic/co2stop/storage_table.csv",
        polygons="resources/automatic/co2stop/storage_map.kml",
        countries="resources/automatic/co2stop/countries.kml",
    output:
        mtco2="resources/automatic/co2stop/{dataset}/{cdr_group}.parquet",
        plot_issues=report(
            "resources/automatic/co2stop/{dataset}/{cdr_group}_kept.png",
            caption="../report/prepare_co2stop_kept.rst",
            category="CO2Stop module",
        ),
        plot_scenarios=report(
            "resources/automatic/co2stop/{dataset}/{cdr_group}_scenarios.png",
            caption="../report/prepare_co2stop_scenarios.rst",
            category="CO2Stop module",
        ),
    log:
        "logs/{dataset}/{cdr_group}/prepare_co2stop.log",
    wildcard_constraints:
        dataset="storage_units",
        cdr_group="aquifer",
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/prepare_co2stop.py"


rule prepare_co2stop_traps:
    message:
        "Harmonising CO2Stop {wildcards.dataset}:{wildcards.cdr_group}."
    params:
        bounds=config["imputation"]["bounds_mtco2"]["co2stop_polygons"],
        filters=config["imputation"]["filter"],
        geo_crs=config["crs"]["geographic"],
    input:
        table=rules.unzip_co2stop.output.traps_data,
        polygons=rules.unzip_co2stop.output.traps_map,
        countries=rules.unzip_co2stop.output.country_map,
    output:
        mtco2="resources/automatic/co2stop/{dataset}/{cdr_group}.parquet",
        plot_issues="resources/automatic/co2stop/{dataset}/{cdr_group}_issues.png",
        plot_scenarios="resources/automatic/co2stop/{dataset}/{cdr_group}_scenarios.png",
    log:
        "logs/{dataset}/{cdr_group}/prepare_co2stop.log",
    wildcard_constraints:
        dataset="traps",
        cdr_group="|".join(CDR_GROUP),
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/prepare_co2stop.py"
