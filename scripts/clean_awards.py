from pandas._libs.missing import NAType
from snakemake.script import snakemake
import pandas as pd
from tqdm.auto import tqdm
import duckdb
from utils.metadata import (
    AwardIdentifiers,
    LocalIdentifierMetadataProvider,
    OpenLibraryMetadataProvider,
)

# load the identifiers data
identifiers = pd.read_parquet(snakemake.input["identifiers"])

local_provider = LocalIdentifierMetadataProvider(identifiers)
openlibrary_provider = OpenLibraryMetadataProvider(local_provider)


def find_work_id(row: pd.Series) -> int | NAType:
    identifiers = AwardIdentifiers.from_award(row)
    if work_id := local_provider.lookup(identifiers):
        return work_id
    if work_id := openlibrary_provider.lookup(identifiers):
        return work_id
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
awards.loc[
    awards.work_id.notna() & ~awards.work_id.isin(set(identifiers.work_id)), "work_id"
] = pd.NA

# fill in the missing work_ids

print(f"{awards.work_id.isna().sum()} awards are missing work_ids")

for i, row in tqdm(
    awards[awards.work_id.isna()].iterrows(),
    total=awards.work_id.isna().sum(),
    desc="Filling in missing work_ids for awards",
):
    awards.loc[i, "work_id"] = find_work_id(row)

print(f"{awards.work_id.isna().sum()} awards are still missing work_ids")

awards.to_csv(snakemake.output[0], index=False)
