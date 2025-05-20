import os.path
from tqdm.auto import tqdm

# Load configuration
configfile: "config.yaml"

rule all:
    input:
        config["paths"]["augmented_works_data"]

rule clean_all:
    input:
        config["paths"]["selected_works_data"],
        config["paths"]["augmented_works_data"],
        config["paths"]["cleaned_awards_data"]

rule download_data:
    input:
        config["paths"]["raw_book_data"],
        config["paths"]["raw_author_data"],
        config["paths"]["raw_awards_data"],
        config["paths"]["raw_works_data"]

rule scrape_raw_awards_data:
    output:
        awards=config["paths"]["raw_awards_data"]
    params:
        rel_path=lambda wildcards, input, output: os.path.relpath(output.awards, "data/bookdata")
    shell:
        "cd data/bookdata && pixi run scrapy crawl sfadb -o {params.rel_path}"

rule clean_awards_data:
    input:
        config["paths"]["raw_awards_data"]
    output:
        config["paths"]["cleaned_awards_data"]
    script:
        "scripts/clean_awards.py"

rule download_raw_goodreads_data:
    output:
        "data/raw/{data_name}"
    params:
        url=lambda wildcards: os.path.join(config["remote"]["goodreads_data"], f"{wildcards.data_name}")
    message:
        "Fetching raw goodreads data: {wildcards.data_name}"
    script:
        "scripts/download_data.py"

rule combine_data:
    input:
        books=config["paths"]["raw_book_data"],
        works=config["paths"]["raw_works_data"]
    output:
        selected_works=config["paths"]["selected_works_data"],
        augmented_works=config["paths"]["augmented_works_data"]
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
        export OUTPUT_SELECTED_WORKS="{output.selected_works}"
        export OUTPUT_AUGMENTED_WORKS="{output.augmented_works}"

        envsubst < scripts/combine_data.sql | duckdb
        """
