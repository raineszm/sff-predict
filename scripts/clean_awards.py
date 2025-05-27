from snakemake.script import snakemake
import duckdb
import pandas as pd

# Reshape the awards data so we have a single row per novel
# and keep track of a count of how many awards each novel has
# been nominated for and won in the year
#
# This could be done with pandas too, but duckdb is faster
# and the expression is easier to write this way for me
awards = duckdb.execute(
    f"""
        PIVOT 
        '{snakemake.input["novels"]}'
        ON status
        GROUP BY
            work_qid,
            title,
            author_qids,
            authors,
            year,
            pubDate,
            openlibrary_ids,
            isbns
        """
).df()

# clean up the column names
awards.rename(
    columns={
        "nominated": "n_nom",
        "winner": "n_win",
    },
    inplace=True,
)

# see how many works are missing an openlibrary id
n_missing = (
    awards.groupby("work_qid").openlibrary_ids.first().replace("", pd.NA).isna().sum()
)

print(f"{n_missing} works are missing an openlibrary id")

awards.to_csv(snakemake.output[0], index=False)
