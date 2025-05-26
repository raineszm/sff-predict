import os.path

# Configuration
# ------------

GOODREADS_DATA_ROOT = "https://mcauleylab.ucsd.edu/public_datasets/gdrive/goodreads"

DATA_ROOT = "data"
RAW_DATA_ROOT = os.path.join(DATA_ROOT, "raw")

GOODREADS_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'books': 'goodreads_books.json.gz',
        'authors': 'goodreads_book_authors.json.gz',
        'works': 'goodreads_book_works.json.gz',
    }.items()
}

WIKIDATA_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'award_novels': 'award_novels.json',
        'awards_all': 'awards_all.json',
        'nominee_biographical': 'nominee_biographical.json',
    }.items()
}

RAW_DATA = {
    **GOODREADS_DATA,
    **WIKIDATA_DATA,
}

PROCESSED_DATA = {
    k: os.path.join(DATA_ROOT, v) for k, v in
    {
        'selected_works': 'sff_works.parquet',
        'selected_books': 'sff_books.parquet',
        'augmented_works': 'sff_works_augmented.parquet',
        'identifiers': 'identifiers.parquet',
        'cleaned_awards': 'awards.csv',
    }.items()
}

# Phony rules
# ------------

rule all:
    input:
        expand(PROCESSED_DATA.values())

rule download_all_data:
    input:
        expand(RAW_DATA.values())

rule download_goodreads_data:
    input:
        expand(GOODREADS_DATA.values())

rule download_raw_awards_data:
    input:
        expand(WIKIDATA_DATA.values())

# Data download rules
# ------------------

rule download_raw_goodreads_data:
    output:
        RAW_DATA_ROOT + "/{basename}"
    wildcard_constraints:
        basename=r"goodreads_\w+\.json\.gz"
    params:
        url=lambda wildcards: os.path.join(GOODREADS_DATA_ROOT, wildcards.basename)
    message:
        "Fetching raw goodreads data: {wildcards.basename}"
    script:
        "scripts/download_data.py"

# Data processing rules
# --------------------

rule combine_data:
    input:
        **RAW_DATA,
    output:
        selected_works=PROCESSED_DATA['selected_works'],
        augmented_works=PROCESSED_DATA['augmented_works'],
        identifiers=PROCESSED_DATA['identifiers'],
        selected_books=PROCESSED_DATA['selected_books']
    message:
        "Filtering and combining input datasets"
    params:
        # Minimum number of ratings for a work to be included
        ratings_threshold=100,
        # Minimum number of times a work must be tagged as scifi/fantasy
        # to be included
        tag_threshold=2
    shell:
        """
        export RATINGS_THRESHOLD={params.ratings_threshold}
        export TAG_THRESHOLD="{params.tag_threshold}"
        export INPUT_WORKS="{input.works}"
        export INPUT_BOOKS="{input.books}"
        export INPUT_AUTHORS="{input.authors}"
        export OUTPUT_SELECTED_WORKS="{output.selected_works}"
        export OUTPUT_AUGMENTED_WORKS="{output.augmented_works}"
        export OUTPUT_IDENTIFIERS="{output.identifiers}"
        export OUTPUT_SELECTED_BOOKS="{output.selected_books}"

        envsubst < scripts/queries/combine_data.sql | duckdb
        """

# Wikidata data download rules
# ---------------------------

for sparql_query, data_name in [
    ('wikidata_awards.sparql', "award_novels"),
    ('wikidata_winners.sparql', "awards_all"),
    ('wikidata_authors.sparql', "nominee_biographical")
]:
    rule:
        name: f'download_{data_name}'
        input:
            sparql=f'scripts/queries/{sparql_query}'
        output:
            json=WIKIDATA_DATA[data_name]
        message:
            f"Downloading {data_name} from Wikidata"
        script:
            "scripts/download_sparql_query.py"

# Awards processing rules
# ----------------------

rule clean_awards_data:
    input:
        awards=WIKIDATA_DATA['award_novels'],
        identifiers=PROCESSED_DATA['identifiers']
    output:
        PROCESSED_DATA['cleaned_awards']
    script:
        "scripts/clean_awards.py"