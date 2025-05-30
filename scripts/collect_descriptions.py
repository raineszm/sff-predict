from snakemake.script import snakemake
import duckdb
from utils.metadata import WikipediaDescriptionProvider, GoogleBooksDescriptionProvider
import numpy as np
from tqdm.auto import tqdm
import pandas as pd
import dotenv
import os
from typing import Callable, Optional


def descriptions_from_ids(
    ids: pd.DataFrame, key_col: str, provider_fn: Callable[[str], Optional[str]]
):
    work_qids = ids.work_qid.unique()
    work_descriptions = pd.DataFrame(
        {
            "work_qid": work_qids,
            "description": np.empty_like(work_qids, dtype="object"),
        }
    ).set_index("work_qid")
    work_groups = ids.groupby("work_qid")
    for work_qid, group in tqdm(
        work_groups,
        desc="Getting descriptions from wikipedia",
        total=len(work_groups),
    ):
        for key in group[key_col]:
            if description := provider_fn(key):
                work_descriptions.loc[work_qid, "description"] = description
                break

    return work_descriptions.dropna().reset_index()


dotenv.load_dotenv()

duckdb.sql(
    f"CREATE OR REPLACE TABLE openlibrary_ids AS FROM '{snakemake.input['openlibrary_ids']}'"
)

openlibrary_descriptions = duckdb.sql(
    # The second part of the where clause
    #  is a hack to remove a duplicate work
    # For some reason there is a seperate work for the Italian translation of
    # 'The Stone Sky' by N.K. Jemisin
    """
SELECT DISTINCT openlibrary_ids.work_qid, openlibrary_works.description
FROM openlibrary_ids
INNER JOIN '{ol_works}' as openlibrary_works
ON '/works/' || openlibrary_ids.openlibrary_id = openlibrary_works.key
WHERE openlibrary_works.description IS NOT NULL
AND openlibrary_ids.openlibrary_id != 'OL28171187W'
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
    FROM read_csv('{wikipedia}', header := true)
    ANTI JOIN openlibrary_descriptions USING (work_qid)
""".format(wikipedia=snakemake.input["wikipedia"])
).df()


with WikipediaDescriptionProvider() as wiki:
    wikipedia_descriptions = descriptions_from_ids(
        wikipedia_ids, "wikipedia_url", wiki.get_description
    )

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
    isbn_descriptions = descriptions_from_ids(
        isbns, "isbn", google_books.get_description_from_isbn
    )

duckdb.register("isbn_descriptions", isbn_descriptions)
print(
    "Fetched {} descriptions from google books by isbn".format(
        isbn_descriptions.shape[0]
    )
)

print("Falling back to searching google books by title and author for remaining works")

remainder = duckdb.sql(
    """
    SELECT DISTINCT work_qid, title || ';' || FIRST(authorLabel) as title_author
    FROM '{nominated_novels}'
    ANTI JOIN openlibrary_descriptions USING (work_qid)
    ANTI JOIN wikipedia_descriptions USING (work_qid)
    ANTI JOIN isbn_descriptions USING (work_qid)
    GROUP BY work_qid, title
""".format(nominated_novels=snakemake.input["nominated_novels"])
).df()

with GoogleBooksDescriptionProvider(os.getenv("GOOGLE_BOOKS_API_KEY")) as google_books:
    remainder_descriptions = descriptions_from_ids(
        remainder,
        "title_author",
        google_books.get_description_from_title_author_string,
    )

print(
    "Fetched {} descriptions from google books by title and author".format(
        remainder_descriptions.shape[0]
    )
)

descriptions = (
    pd.concat(
        [
            openlibrary_descriptions,
            wikipedia_descriptions,
            isbn_descriptions,
            remainder_descriptions,
        ]
    ).drop_duplicates(subset=["work_qid"])
).set_index("work_qid")

n_works = duckdb.sql(
    "SELECT COUNT(DISTINCT work_qid) FROM read_csv('{nominated_novels}', header := true);".format(
        nominated_novels=snakemake.input["nominated_novels"]
    )
).fetchone()[0]
print("Descriptions found for {}/{} works".format(descriptions.index.size, n_works))

descriptions.to_csv(snakemake.output["descriptions"])
