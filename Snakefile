import os.path

# Configuration
# ------------

DATA_ROOT = "data"
RAW_DATA_ROOT = os.path.join(DATA_ROOT, "raw")

WIKIDATA_DATA = {
    k: os.path.join(RAW_DATA_ROOT, v) for k, v in
    {
        'award_novels': 'award_novels.json',
        'awards_all': 'awards_all.json',
        'nominee_biographical': 'nominee_biographical.json',
    }.items()
}

RAW_DATA = {
    **WIKIDATA_DATA,
}

PROCESSED_DATA = {
    k: os.path.join(DATA_ROOT, v) for k, v in
    {
        'nominated_novels': 'nominated_novels.csv',
        'cumulative_awards': 'cumulative_awards.csv',
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

# Awards processing rules
# ----------------------

rule clean_nominated_novels:
    input:
        novels=WIKIDATA_DATA['award_novels']
    output:
        PROCESSED_DATA['nominated_novels']
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