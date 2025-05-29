from snakemake.script import snakemake
import duckdb
from utils.metadata import WikipediaDescriptionProvider
import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import requests_cache

duckdb.sql(
    f"CREATE OR REPLACE TABLE openlibrary_ids AS FROM '{snakemake.input['openlibrary_ids']}'"
)

openlibrary_descriptions = duckdb.sql(
    """
SELECT DISTINCT openlibrary_ids.work_qid, openlibrary_works.description
FROM openlibrary_ids
LEFT JOIN '{ol_works}' as openlibrary_works
ON '/works/' || openlibrary_ids.openlibrary_id = openlibrary_works.key
WHERE openlibrary_works.description IS NOT NULL
""".format(
        ol_works=snakemake.input["ol_works"],
    )
).df()

print(
    "Fetched {} descriptions from openlibrary data".format(
        openlibrary_descriptions.shape[0]
    )
)

wikipedia_ids = duckdb.sql(
    """
    SELECT DISTINCT work_qid, wikipedia_url
    FROM '{wikipedia}'
    ANTI JOIN openlibrary_ids USING (work_qid)
""".format(wikipedia=snakemake.input["wikipedia"])
).df()

wikipedia_ids["description"] = np.empty_like(wikipedia_ids.wikipedia_url)

with requests_cache.enabled(
    backend="filesystem",
    cache_control=False,
    cache_name=".cache/wikipedia",
    expire_after=requests_cache.NEVER_EXPIRE,
):
    wiki = WikipediaDescriptionProvider()
    for i, row in tqdm(
        wikipedia_ids.iterrows(),
        desc="Getting descriptions from wikipedia",
        total=len(wikipedia_ids),
    ):
        wikipedia_ids.loc[i, "description"] = wiki.get_description(row["wikipedia_url"])

wikipedia_ids = wikipedia_ids.drop(columns=["wikipedia_url"])

print("Fetched {} descriptions from wikipedia".format(wikipedia_ids.shape[0]))


descriptions = (
    pd.concat([openlibrary_descriptions, wikipedia_ids])
    .dropna(subset=["description"])
    .drop_duplicates(subset=["work_qid"])
).set_index("work_qid")

n_works = duckdb.sql(
    "SELECT DISTINCT COUNT(work_qid) FROM '{nominated_novels}';".format(
        nominated_novels=snakemake.input["nominated_novels"]
    )
).fetchone()[0]
print("Descriptions found for {}/{} works".format(len(descriptions), n_works))

descriptions.to_csv(snakemake.output["descriptions"])
