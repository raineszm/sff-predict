from snakemake.script import snakemake
import duckdb
from utils.metadata import WikipediaDescriptionProvider, GoogleBooksDescriptionProvider
import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import dotenv
import os

dotenv.load_dotenv()

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
duckdb.register("openlibrary_descriptions", openlibrary_descriptions)

print(
    "Fetched {} descriptions from openlibrary data".format(
        openlibrary_descriptions.shape[0]
    )
)

wikipedia_ids = duckdb.sql(
    """
    SELECT DISTINCT work_qid, wikipedia_url
    FROM '{wikipedia}'
    ANTI JOIN openlibrary_descriptions USING (work_qid)
""".format(wikipedia=snakemake.input["wikipedia"])
).df()

wikipedia_ids["description"] = np.empty_like(wikipedia_ids.wikipedia_url)

with WikipediaDescriptionProvider() as wiki:
    for i, row in tqdm(
        wikipedia_ids.iterrows(),
        desc="Getting descriptions from wikipedia",
        total=len(wikipedia_ids),
    ):
        wikipedia_ids.loc[i, "description"] = wiki.get_description(row["wikipedia_url"])

wikipedia_descriptions = wikipedia_ids.dropna().drop(columns=["wikipedia_url"])
duckdb.register("wikipedia_descriptions", wikipedia_descriptions)

print("Fetched {} descriptions from wikipedia".format(wikipedia_descriptions.shape[0]))

isbns = duckdb.sql(
    """
    SELECT DISTINCT work_qid, isbn
    FROM '{isbns}'
    ANTI JOIN openlibrary_descriptions USING (work_qid)
    ANTI JOIN wikipedia_descriptions USING (work_qid)
""".format(isbns=snakemake.input["isbns"])
).df()

with GoogleBooksDescriptionProvider(os.getenv("GOOGLE_BOOKS_API_KEY")) as google_books:
    work_qids = isbns.work_qid.unique()
    isbn_descriptions = pd.DataFrame(
        {
            "work_qid": work_qids,
            "description": np.empty_like(work_qids, dtype="object"),
        }
    ).set_index("work_qid")
    isbn_groups = isbns.groupby("work_qid")
    for work_qid, group in tqdm(
        isbn_groups,
        desc="Getting descriptions from google books",
        total=len(isbn_groups),
    ):
        for isbn in group["isbn"]:
            if description := google_books.get_description(isbn):
                isbn_descriptions.loc[work_qid, "description"] = description
                break


isbn_descriptions = isbn_descriptions.dropna().reset_index()

print("Fetched {} descriptions from google books".format(isbn_descriptions.shape[0]))


descriptions = (
    pd.concat(
        [openlibrary_descriptions, wikipedia_descriptions, isbn_descriptions]
    ).drop_duplicates(subset=["work_qid"])
).set_index("work_qid")

n_works = duckdb.sql(
    "SELECT COUNT(DISTINCT work_qid) FROM '{nominated_novels}';".format(
        nominated_novels=snakemake.input["nominated_novels"]
    )
).fetchone()[0]
print("Descriptions found for {}/{} works".format(descriptions.index.size, n_works))

descriptions.to_csv(snakemake.output["descriptions"])
