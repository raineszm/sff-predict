from typing import Optional
from snakemake.script import snakemake
import pandas as pd
import httpx
from tqdm.auto import tqdm
from httpx_retries import RetryTransport, Retry
import duckdb

# load the identifiers data
identifiers = pd.read_parquet(snakemake.input["identifiers"])

WORK_IDS = set(identifiers.work_id)

# Be a good citizen
retry = Retry(
    total=3,
    backoff_factor=0.5,
    status_forcelist=[429],
)

transport = RetryTransport(
    retry=retry,
)

client = httpx.Client(
    http2=True,
    headers={
        "User-Agent": "scifi-fantasy/0.1 (dev@zmraines.com)",
    },
    transport=transport,
)


def query_openlibrary(query: str) -> dict:
    response = client.get(
        "https://openlibrary.org/search.json",
        params={
            "q": query,
            "fields": "key,id_goodreads",
            "limit": 1,
        },
    )
    response.raise_for_status()
    json = response.json()
    if json["numFound"] == 0 or "id_goodreads" not in json["docs"][0]:
        return None

    for work_id in json["docs"][0]["id_goodreads"]:
        # check if this is a valid work id
        if int(work_id) in WORK_IDS:
            return int(work_id)

        # sometimes openlibrary returns a book id instead of a work id
        # so we check if this is a valid book id and if so, return the work id
        by_book_id = identifiers.loc[
            (identifiers.kind == "book_id") & (identifiers.value == work_id)
        ].work_id
        if not by_book_id.empty:
            return by_book_id.iloc[0]


def find_work_id_from_ol_ids(ol_ids: list[str]) -> Optional[str]:
    for ol_id in ol_ids:
        work_id = query_openlibrary(f"key:/works/{ol_id}")
        if work_id is not None:
            return work_id
    return None


def find_work_id_from_isbns(isbns: list[str]) -> Optional[str]:
    for isbn in isbns:
        isbn = isbn.replace("-", "")
        work_id = query_openlibrary(f"isbn:{isbn}")
        if work_id is not None:
            return work_id
    return None


def find_work_id(row: pd.Series) -> Optional[str]:
    if row.ol_ids != "":
        return find_work_id_from_ol_ids(row.ol_ids.split(";"))
    if row.isbns != "":
        return find_work_id_from_isbns(row.isbns.split(";"))
    return pd.NA


# Reshape the awards data so we have a single row per novel
# and keep track of a count of how many awards each novel has
# been nominated for and won in the year
#
# This could be done with pandas too, but duckdb is faster
# and the expression is easier to write this way for me
awards = duckdb.execute(
    f"""
        PIVOT 
        '{snakemake.input["awards"]}'
        ON status
        GROUP BY
            title,
            authors,
            year,
            work_ids,
            ol_ids,
            isbns,
            isdfb_ids;
        """
).df()

# clean up the column names
awards.rename(
    columns={
        "work_ids": "work_id",
        "nominated": "n_nom",
        "winner": "n_win",
    },
    inplace=True,
)


# If a novel has multiple work ids, we only keep the first one
has_multiple_work_ids = awards.work_id.str.contains(";")
awards.loc[has_multiple_work_ids, "work_id"] = awards.loc[
    has_multiple_work_ids, "work_id"
].map(lambda x: x.split(";")[0])


# convert the work_id column to an integer type
# and set any values that are not valid work ids to NA
awards.work_id = awards.work_id.replace("", pd.NA).astype("Int64")
awards.loc[awards.work_id.notna() & ~awards.work_id.isin(WORK_IDS), "work_id"] = pd.NA

print(f"{awards.work_id.isna().sum()} awards are missing work_ids")
# fill in the missing work_ids
for i, row in tqdm(
    awards[awards.work_id.isna()].iterrows(),
    total=awards.work_id.isna().sum(),
    desc="Filling in missing work_ids for awards",
):
    awards.loc[i, "work_id"] = find_work_id(row)

print(f"{awards.work_id.isna().sum()} awards are still missing work_ids")

awards.to_csv(snakemake.output[0], index=False)
