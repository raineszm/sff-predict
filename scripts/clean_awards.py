from typing import Optional
from snakemake.script import snakemake
import pandas as pd
import httpx
from tqdm.auto import tqdm
from httpx_retries import RetryTransport, Retry

# we'll need to test for what works exist
works = pd.read_json(snakemake.input["works"], lines=True)

WORK_IDS = set(works.work_id)

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


def find_work_id_from_ol_ids(ol_ids: list[str]) -> Optional[str]:
    for ol_id in ol_ids:
        response = client.get(
            "https://openlibrary.org/search.json",
            params={
                "q": f"key:/works/{ol_id}",
                "fields": "key,id_goodreads",
                "limit": 1,
            },
        )
        response.raise_for_status()
        json = response.json()
        if json["numFound"] == 0 or "id_goodreads" not in json["docs"][0]:
            continue

        for work_id in json["docs"][0]["id_goodreads"]:
            if int(work_id) in WORK_IDS:
                return int(work_id)


def find_work_id(row: pd.Series) -> Optional[str]:
    if row.ol_ids != "":
        return find_work_id_from_ol_ids(row.ol_ids.split(";"))
    return pd.NA


# Let's clean the awards data
# so we can match it with the books data
awards = pd.read_json(snakemake.input["awards"], lines=True)

awards.rename(columns={"work_ids": "work_id"}, inplace=True)

has_multiple_work_ids = awards.work_id.str.contains(";")
awards.loc[has_multiple_work_ids, "work_id"] = awards.loc[
    has_multiple_work_ids, "work_id"
].map(lambda x: x.split(";")[0])

has_work_id = awards.work_id != ""
awards.loc[has_work_id, "work_id"] = awards.loc[has_work_id, "work_id"].astype("Int64")

for i, row in tqdm(
    awards[~has_work_id].iterrows(),
    total=(~has_work_id).sum(),
    desc="Filling in missing work_ids for awards",
):
    awards.loc[i, "work_id"] = find_work_id(row)

awards.to_json(snakemake.output[0], orient="records", lines=True)
