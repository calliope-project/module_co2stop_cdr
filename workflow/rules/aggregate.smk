"""Rules for shape aggregation."""

rule aggregate_co2stop:
    message:
        "Aggregating configured inputation to '{wildcards.shapes}-{wildcards.scenario}'."
    params:
        proj_crs=config["crs"]["projected"],
        included_datasets=config["imputation"]["datasets"],
        included_groups=config["imputation"]["groups"]
    input:
        shapes="resources/user/{shapes}/shapes.parquet",
        storage_units="resources/automatic/co2stop/storage_units/aquifer.parquet",
        all_traps=expand("resources/automatic/co2stop/traps/{cdr_group}.parquet", cdr_group=CDR_GROUP),
    output:
        aggregated="results/{shapes}/{scenario}/aggregated.parquet",
        plot="results/{shapes}/{scenario}/aggregated.png"
    log:
        "logs/{shapes}/{scenario}/aggregate_co2stop.log",
    wildcard_constraints:
        scenario="|".join(["low", "medium", "high"])
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/aggregate_co2stop.py"
