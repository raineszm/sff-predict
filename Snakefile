import os.path

# Configuration
# ------------

DATA_ROOT = "data"
RAW_DATA_ROOT = os.path.join(DATA_ROOT, "raw")

# Model configurations
# ------------
EMBEDDING_MODEL = "all-MiniLM-L6-v2"
SENTIMENT_MODEL = {
    "model": "distilbert-base-uncased-finetuned-sst-2-english",
    "revision": "714eb0f",
    "truncation": True,
    "max_length": 512
}

WIKIDATA_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'award_novels': 'award_novels.json',
        'awards_all': 'awards_all.json',
        'nominee_biographical': 'nominee_biographical.json',
    }.items()
}

KAGGLE_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'openlibrary_works': 'ol_works.parquet',
    }.items()
}

NYT_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'headlines': 'headlines.json',
    }.items()
}

HEADLINE_DIR = os.path.join(RAW_DATA_ROOT, "headlines")

RAW_DATA = {
    **WIKIDATA_DATA,
    **KAGGLE_DATA,
    **NYT_DATA,
}

TRAIN_DATA = {
    k: os.path.join(DATA_ROOT, v) for k, v in
    {
        'train_desc': 'train_desc_labeled.csv',
    }.items()
}

PROCESSED_DATA = {
    k: os.path.join(DATA_ROOT, v) for k, v in
    {
        'nominated_novels': 'nominated_novels.csv',
        'openlibrary_ids': 'openlibrary_ids.csv',
        'isbns': 'isbns.csv',
        'wikipedia': 'wikipedia.csv',
        'cumulative_awards': 'cumulative_awards.csv',
        'descriptions': 'descriptions.csv',
        'headline_embeddings': 'headline_embeddings.parquet',
        'description_embeddings': 'description_embeddings.parquet',
        'descriptions_debiased': 'descriptions_debiased.parquet',
    }.items()
}

# General targets
# ------------

rule all:
    input:
        expand(PROCESSED_DATA.values())

rule download_all_data:
    input:
        expand(RAW_DATA.values())

rule download_raw_awards_data:
    input:
        expand(WIKIDATA_DATA.values())

# ---------------------------
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


# Kaggle data download rules
# -------------------------

rule download_openlibrary_works:
    output:
        protected(KAGGLE_DATA['openlibrary_works'])
    script:
        "scripts/download_openlibrary.py"


# Awards processing rules
# ----------------------

rule clean_nominated_novels:
    input:
        novels=WIKIDATA_DATA['award_novels'],
        wins_as_of=PROCESSED_DATA['cumulative_awards'],
        authors=WIKIDATA_DATA['nominee_biographical']
    output:
        cleaned_novels=PROCESSED_DATA['nominated_novels'],
        openlibrary_ids=PROCESSED_DATA['openlibrary_ids'],
        isbns=PROCESSED_DATA['isbns'],
        wikipedia=PROCESSED_DATA['wikipedia']
    script:
        "scripts/clean_awards.py"

rule tally_cumulative_awards:
    input:
        all_awards=WIKIDATA_DATA['awards_all'],
        sql=f'scripts/queries/awards_as_of.sql'
    output:
        cumulative_awards=PROCESSED_DATA['cumulative_awards']
    shell:
        """
        export INPUT_AWARDS={input.all_awards}
        export OUTPUT_CUMULATIVE={output.cumulative_awards}
        envsubst < {input.sql} | duckdb
        """

rule collect_descriptions:
    input:
        openlibrary_ids=PROCESSED_DATA['openlibrary_ids'],
        ol_works=KAGGLE_DATA['openlibrary_works'],
        nominated_novels=PROCESSED_DATA['nominated_novels'],
        wikipedia=PROCESSED_DATA['wikipedia'],
        isbns=PROCESSED_DATA['isbns']
    output:
        descriptions=PROCESSED_DATA['descriptions']
    script:
        "scripts/collect_descriptions.py"

rule embed_descriptions:
    input:
        descriptions=PROCESSED_DATA['descriptions']
    output:
        description_embeddings=PROCESSED_DATA['description_embeddings']
    params:
        embedding_model=EMBEDDING_MODEL
    script:
        "scripts/embed_descriptions.py"

rule debias_descriptions:
    input:
        description_embeddings=PROCESSED_DATA['description_embeddings'],
        train_desc=TRAIN_DATA['train_desc']
    output:
        descriptions_debiased=PROCESSED_DATA['descriptions_debiased']
    script:
        "scripts/debias_descriptions.py"

# World State download rules
# --------------------------

rule fetch_month_headlines:
    output:
        month_headlines=os.path.join(HEADLINE_DIR, "{year}-{month}.json")
    resources:
        nyt_api=1
    params:
        year="{year}",
        month="{month}"
    script:
        "scripts/fetch_month_headlines.py"

rule download_headlines:
    input:
        months=expand(
            "{HEADLINE_DIR}/{year}-{month}.json",
            HEADLINE_DIR=HEADLINE_DIR,
            year=range(1954, 2025),
            month=range(1, 13)
        )
    log:
        "logs/download_headlines.log"
    output:
        headlines=NYT_DATA['headlines']
    run:
        with open(output.headlines, "w") as outfile:
            for month_headlines in input.months:
                with open(month_headlines, "r") as infile:
                    for line in infile:
                        outfile.write(line)


rule embed_headlines:
    input:
        headlines=NYT_DATA['headlines']
    output:
        headline_embeddings=protected(PROCESSED_DATA['headline_embeddings'])
    params:
        embedding_model=EMBEDDING_MODEL,
        sentiment_model=SENTIMENT_MODEL
    script:
        "scripts/embed_headlines.py"
