import os.path

# Load configuration
configfile: "config.yaml"

def remote_data_url(basename):
    """Get the full URL for a remote data file.
    
NOTE: this function takes the BASENAME not the key because
it is called with a matched wildcard"""
    return os.path.join(config["remote"]["goodreads_data_root"], basename)

def raw_data_path(key):
    """Get the full path for a raw data file."""
    basename = config["paths"]["basenames"][key]
    return os.path.join(config["paths"]["roots"]["raw_data"], basename)

def local_data_path(key):
    """Get the full path for a processed data file."""
    basename = config["paths"]["basenames"][key]
    return os.path.join(config["paths"]["roots"]["data"], basename)

rule all:
    input:
        local_data_path("augmented_works_data")

rule clean_all:
    input:
        local_data_path("selected_works_data"),
        local_data_path("augmented_works_data"),
        local_data_path("cleaned_awards_data")

rule download_data:
    input:
        raw_data_path("book_data"),
        raw_data_path("author_data"),
        raw_data_path("works_data"),
        raw_data_path("awards_data")

rule download_raw_awards_data:
    input:
        sparql='scripts/wikidata_awards.sparql'
    output:
        awards=raw_data_path("awards_data")
    message:
        "Downloading awards data from Wikidata"
    script:
        "scripts/download_wikidata_awards.py"
rule clean_awards_data:
    input:
        awards=raw_data_path("awards_data"),
        works=raw_data_path("works_data")
    output:
        local_data_path("cleaned_awards_data")
    script:
        "scripts/clean_awards.py"

rule download_raw_goodreads_data:
    output:
        "data/raw/{basename}"
    params:
        url=lambda wildcards: remote_data_url(wildcards.basename)
    message:
        "Fetching raw goodreads data: {wildcards.basename}"
    script:
        "scripts/download_data.py"

rule combine_data:
    input:
        books=raw_data_path("book_data"),
        works=raw_data_path("works_data"),
        authors=raw_data_path("author_data")
    output:
        selected_works=local_data_path("selected_works_data"),
        augmented_works=local_data_path("augmented_works_data")
    message:
        "Filtering and combining input datasets"
    params:
        # Minimum number of ratings for a work to be included
        ratings_threshold=config["filter"]["ratings_threshold"],
        # Minimum number of times a work must be tagged as scifi/fantasy
        # to be included
        tag_threshold=config["filter"]["tag_threshold"]
    shell:
        """
        export RATINGS_THRESHOLD={params.ratings_threshold}
        export TAG_THRESHOLD="{params.tag_threshold}"
        export INPUT_WORKS="{input.works}"
        export INPUT_BOOKS="{input.books}"
        export INPUT_AUTHORS="{input.authors}"
        export OUTPUT_SELECTED_WORKS="{output.selected_works}"
        export OUTPUT_AUGMENTED_WORKS="{output.augmented_works}"

        envsubst < scripts/combine_data.sql | duckdb
        """
