from snakemake.script import snakemake
import duckdb
import pandas as pd

# Reshape the awards data so we have a single row per novel
# and keep track of a count of how many awards each novel has
# been nominated for and won in the year
#
# This could be done with pandas too, but duckdb is faster
# and the expression is easier to write this way for me
#
# Since we aren't grouping by the nominated and winner columns,
# they are counted by the pivot operation
awards = (
    duckdb.execute(
        f"""
        PIVOT '{snakemake.input["novels"]}'
        ON status
        GROUP BY work_qid, title, author_qids, authors, year, pubDate, openlibrary_ids, isbns
        """
    )
    .df()
    .rename(columns={"nominated": "n_nom", "winner": "n_win"})
)

# Some clean up of the awards data
awards["year"] = awards["year"].astype(int)
awards.isbns = awards.isbns.str.replace("-", "").replace("", pd.NA).str.split(";")
awards.openlibrary_ids = awards.openlibrary_ids.replace("", pd.NA).str.split(";")

# Read cumulative awards
cumulative_awards = pd.read_csv(snakemake.input["wins_as_of"])

# Explode author QIDs and join with cumulative awards
work_author_join = (
    awards.assign(author_qid=awards.author_qids.str.split(";"))
    .explode("author_qid")[["work_qid", "author_qid", "year"]]
    .rename(columns={"year": "year_of_award"})
)

merged = work_author_join.merge(cumulative_awards, on="author_qid")

# Find, for each work and author, the latest cumulative award year before this award
previous_awards = (
    merged[merged.year_of_award > merged.year]
    .sort_values("year")
    .groupby(["work_qid", "author_qid"], as_index=False)
    .last()
    .groupby("work_qid", as_index=False)["awards_as_of_year"]
    .sum()
)

# Join with awards DataFrame
awards = awards.merge(previous_awards, on="work_qid", how="left")
awards.awards_as_of_year = awards.awards_as_of_year.fillna(0).astype(int)

# Add author biography information
authors = pd.read_json(snakemake.input["authors"], lines=True)
authors.dob = pd.to_datetime(authors.dob, format="%Y-%m-%dT%H:%M:%SZ")

author_info = (
    work_author_join.merge(authors, on="author_qid")
    .assign(age=lambda x: x.year_of_award - x.dob.dt.year)
    .groupby("work_qid")
    .agg(
        {
            "genderLabel": list,
            "birthCountryLabel": lambda x: x.dropna().tolist(),
            "age": list,
        }
    )
    .reset_index()
    .rename(
        columns={
            "age": "ages",
            "birthCountryLabel": "nationalities",
            "genderLabel": "genders",
        }
    )
)

awards = awards.merge(author_info, on="work_qid", how="left")

# Count works missing an openlibrary ID
n_missing = (
    awards.groupby("work_qid")["openlibrary_ids"]
    .first()
    .replace("", pd.NA)
    .isna()
    .sum()
)
print(f"{n_missing} works are missing an openlibrary id")

awards.to_json(snakemake.output[0], orient="records", lines=True)
