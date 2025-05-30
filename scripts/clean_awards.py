from snakemake.script import snakemake
import pandas as pd


def get_identifier_table(df, id_col, index_col=None):
    table = df[["work_qid", id_col]].dropna().drop_duplicates(subset=[id_col])
    if index_col is not None:
        table = table.set_index(index_col)
    return table


awards = pd.read_json(snakemake.input["novels"], lines=True).rename(
    columns={"itemLabel": "title", "article": "wikipedia_url"}
)

# Split off identifying information into separate dataframes
openlibrary_ids = get_identifier_table(
    awards, "openlibrary_id", index_col="openlibrary_id"
)
isbns = get_identifier_table(awards, "isbn")
wikipedia = get_identifier_table(awards, "wikipedia_url", index_col="wikipedia_url")

openlibrary_ids.to_csv(snakemake.output["openlibrary_ids"])
isbns.to_csv(snakemake.output["isbns"], index=False)
wikipedia.to_csv(snakemake.output["wikipedia"])

# Drop identifying information from awards dataframe
awards = awards.drop(
    columns=["openlibrary_id", "isbn", "wikipedia_url"]
).drop_duplicates()

group_keys = awards.columns.drop(["awardLabel", "status", "year"]).tolist()

awards = (
    awards.groupby(group_keys, dropna=False)
    .agg(
        # handle the edge cases where a novel's award year is not
        # the calendar year
        year=("year", lambda x: x.mode().iat[0]),
        n_nom=("status", lambda s: (s == "nominated").sum()),
        n_win=("status", lambda s: (s == "winner").sum()),
    )
    .reset_index()
)

# Read cumulative awards
cumulative_awards = pd.read_csv(snakemake.input["wins_as_of"])

work_author_join = awards[["work_qid", "author_qid", "year"]].rename(
    columns={"year": "year_of_award"}
)

merged = work_author_join.merge(cumulative_awards, on="author_qid", how="left")


# Find, for each work and author, the latest cumulative award year before this award
previous_awards = (
    merged[merged.year_of_award > merged.year]
    .sort_values("year")
    .groupby(["work_qid", "author_qid"], as_index=False)
    .last()
    .reset_index()
)[["work_qid", "author_qid", "awards_as_of_year"]]

# Join with awards DataFrame
awards = awards.merge(previous_awards, on=["work_qid", "author_qid"], how="left")
awards.awards_as_of_year = awards.awards_as_of_year.fillna(0).astype(int)

# Add author biography information
authors = pd.read_json(snakemake.input["authors"], lines=True)
authors.dob = pd.to_datetime(authors.dob, format="%Y-%m-%dT%H:%M:%SZ")

author_info = (
    work_author_join.merge(authors, on="author_qid")
    .assign(age=lambda x: x.year_of_award - x.dob.dt.year)
    .rename(
        columns={
            "birthCountryLabel": "birth_country",
            "genderLabel": "gender",
        }
    )
)[["work_qid", "author_qid", "age", "birth_country", "gender"]]

awards = awards.merge(author_info, on=["work_qid", "author_qid"], how="left")

# Count number of works which have no identifying information
# i.e., no openlibrary id or isbn
n_missing = (
    awards.work_qid.nunique()
    - pd.concat(
        [openlibrary_ids.work_qid, isbns.work_qid, wikipedia.work_qid]
    ).nunique()
)
print(
    f"{n_missing}/{awards.work_qid.nunique()} works lack identifying information (openlibrary/isbn/etc.)"
)

awards.to_csv(snakemake.output["cleaned_novels"], index=False)
