import os.path
from tqdm.auto import tqdm

# Load configuration
configfile: "config.yaml"

rule all:
    input:
        config["paths"]["filtered_book_data"]

rule clean_all:
    input:
        config["paths"]["filtered_book_data"],
        config["paths"]["filtered_review_data"],
        config["paths"]["cleaned_awards_data"]

rule download_data:
    input:
        config["paths"]["raw_book_data"],
        config["paths"]["raw_review_data"],
        config["paths"]["raw_author_data"],
        config["paths"]["raw_awards_data"]

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

rule filter_scifi_fantasy:
    input:
        config["paths"]["raw_book_data"]
    output:
        config["paths"]["filtered_book_data"]
    message:
        "Filtering for top scifi/fantasy books"
    script:
        "scripts/filter_books.py" 

rule filter_reviews:
    input:
        config["paths"]["raw_review_data"],
        books=config["paths"]["filtered_book_data"]
    output:
        config["paths"]["filtered_review_data"]
    script:
        "scripts/filter_reviews.py"