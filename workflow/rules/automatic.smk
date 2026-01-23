"""Rules to used to download automatic resource files."""


rule download_co2stop:
    message:
        "Downloading the open CO2Stop dataset."
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


rule unzip_co2stop:
    message:
        "Unzipping necessary data from the CO2Stop file."
    params:
        storage_table="CO2Stop_DataInterrogationSystem/Hydrocarbon_Storage_Units.csv",
        storage_map="CO2Stop_Polygons Data/StorageUnits_March13.kml",
        traps_table1="CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps.csv",
        traps_table2="CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps_Temp.csv",
        traps_table3="CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps1.csv",
        traps_map="CO2Stop_Polygons Data/DaughterUnits_March13.kml",
    input:
        zipfile=rules.download_co2stop.output.zipfile,
    output:
        storage_table="resources/automatic/co2stop/storage_table.csv",
        storage_map="resources/automatic/co2stop/storage_map.kml",
        traps_table1="resources/automatic/co2stop/traps_table1.csv",
        traps_table2="resources/automatic/co2stop/traps_table2.csv",
        traps_table3="resources/automatic/co2stop/traps_table3.csv",
        traps_map="resources/automatic/co2stop/traps_map.kml",
    log:
        "logs/automatic/unzip_co2stop.log",
    conda:
        "../envs/co2stop.yaml"
    script:
        "../scripts/unzip_co2stop.py"
