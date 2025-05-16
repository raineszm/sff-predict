# Load configuration
configfile: "config.yaml"

rule all:
    input:
        config["paths"]["filtered_data"]

rule get_raw_data:
    input:
        storage.http(config["remote"]["bookdata"])
    output:
        config["paths"]["raw_data"]
    message:
        "Fetching raw data"
    shell:
        "gunzip -c {input} > {output}"

rule filter_scifi_fantasy:
    input:
        config["paths"]["raw_data"]
    output:
        config["paths"]["filtered_data"]
    message:
        "Filtering for top scifi/fantasy books"
    script:
        "scripts/filter_books.py" 