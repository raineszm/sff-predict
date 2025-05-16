import os.path
from tqdm.auto import tqdm

# Load configuration
configfile: "config.yaml"

rule all:
    input:
        config["paths"]["filtered_book_data"]

rule download_data:
    input:
        config["paths"]["raw_book_data"],
        config["paths"]["raw_review_data"]

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