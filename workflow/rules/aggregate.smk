"""Rules for shape aggregation."""


rule aggregate_co2stop:
    message:
        "Aggregating '{wildcards.shapes}-{wildcards.scenario}-{wildcards.cdr_group}'."
    params:
        bounds_mtco2=config["imputation"]["bounds_mtco2"]["aggregated_polygons"],
        proj_crs=config["crs"]["projected"],
    input:
        shapes="resources/user/{shapes}/shapes.parquet",
        storage_units="resources/automatic/co2stop/storage_units/aquifer.parquet",
        traps="resources/automatic/co2stop/traps/{cdr_group}.parquet",
    output:
        aggregated="results/{shapes}/{scenario}/{cdr_group}.parquet",
        plot=report(
            "results/{shapes}/{scenario}/{cdr_group}.png",
            caption="../report/aggregate_co2stop.rst",
            category="CO2Stop module",
            subcategory="aggregated {cdr_group}",
        ),
    log:
        "logs/{shapes}/{scenario}/{cdr_group}/aggregate_co2stop.log",
    wildcard_constraints:
        scenario="|".join(["low", "medium", "high"]),
        cdr_group="|".join(CDR_GROUP),
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/aggregate_co2stop.py"


rule aggregate_totals:
    message:
        "Aggregating totals for '{wildcards.shapes}-{wildcards.scenario}'."
    params:
        proj_crs=config["crs"]["projected"],
    input:
        shapes="resources/user/{shapes}/shapes.parquet",
        aggregates=expand(
            "results/{{shapes}}/{{scenario}}/{cdr_group}.parquet", cdr_group=CDR_GROUP
        ),
    output:
        totals="results/{shapes}/{scenario}/totals.parquet",
        plot=report(
            "results/{shapes}/{scenario}/totals.png",
            caption="../report/aggregate_co2stop.rst",
            category="CO2Stop module",
            subcategory="aggregated totals",
        ),
    log:
        "logs/{shapes}/{scenario}/totals/aggregate_co2stop.log",
    wildcard_constraints:
        scenario="|".join(["low", "medium", "high"]),
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/aggregate_totals.py"
