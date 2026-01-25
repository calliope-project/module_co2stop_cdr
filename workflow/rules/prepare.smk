"""Rules for database preparation."""

rule prepare_co2stop:
    message:
        "Harmonising CO2Stop dataset."
    params:
        geo_crs=config["crs"]["geographic"],
    input:
        storage_data=rules.unzip_co2stop.output.storage_data,
        storage_map=rules.unzip_co2stop.output.storage_map,
        traps_data=rules.unzip_co2stop.output.traps_data,
        traps_map=rules.unzip_co2stop.output.traps_map,
    output:
        cdr_potential="resources/automatic/co2stop.csv"
    log:
        "logs/prepare_co2stop.log",
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/prepare_co2stop.py"
