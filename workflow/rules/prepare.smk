"""Rules for database preparation."""

rule prepare_co2stop_storage_units:
    message:
        "Harmonising CO2Stop {wildcards.dataset}:{wildcards.cdr_group}."
    params:
        geo_crs=config["crs"]["geographic"],
    input:
        table="resources/automatic/co2stop/storage_table.csv",
        polygons="resources/automatic/co2stop/storage_map.kml",
        countries="resources/automatic/co2stop/countries.kml",
    output:
        mtco2="resources/automatic/co2stop/{dataset}/{cdr_group}.parquet",
        plot_issues="resources/automatic/co2stop/{dataset}/{cdr_group}_issues.png",
    log:
        "logs/{dataset}/{cdr_group}/prepare_co2stop.log",
    wildcard_constraints:
        dataset="storage_units",
        cdr_group="aquifer"
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/prepare_co2stop.py"
