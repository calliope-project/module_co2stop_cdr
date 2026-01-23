"""Rules to used to download automatic resource files."""


rule download_co2stop:
    message:
        "Download the open CO2Stop dataset."
    params:
        url=internal["resources"]["automatic"]["co2stop"],
    output:
        zipfile="resources/automatic/co2stop.zip",
    log:
        "logs/download_co2stop.log",
    conda:
        "../envs/shell.yaml"
    shell:
        'curl -sSLo {output.zipfile:q} {params.url:q}'
