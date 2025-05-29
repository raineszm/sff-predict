from snakemake.script import snakemake
import pandas as pd


def explode_id_column(df, id_column, single_id_col=None):
    """
    Takes a dataframe with an ID column containing semicolon-separated values and explodes it into separate rows.

    Args:
        df: DataFrame containing the ID column
        id_column: Name of column containing semicolon-separated IDs (should end in 's')

    Returns:
        DataFrame with exploded ID column set as index
    """
    if not id_column.endswith("s"):
        raise ValueError(f"ID column {id_column} should end in 's'")
    if single_id_col is None:
        single_id_col = id_column[:-1]  # Trim trailing 's'
    return (
        df.loc[df[id_column] != "", ["work_qid", id_column]]
        .assign(**{single_id_col: lambda x: x[id_column].str.split(";")})
        .explode(single_id_col)
        .drop(columns=[id_column])
        .drop_duplicates()
        .set_index(single_id_col)
    )


# Reshape the awards data so we have a single row per novel
# and keep track of a count of how many awards each novel has
# been nominated for and won in the year
#
awards = pd.read_json(snakemake.input["novels"], lines=True).drop(
    columns=["awardLabel"]
)
group_keys = awards.columns.drop(["status", "year"]).tolist()

awards = (
    awards.groupby(group_keys)
    .agg(
        # handle the edge cases where a novel's award year is not
        # the calendar year
        year=("year", lambda x: x.mode().iat[0]),
        n_nom=("status", lambda s: (s == "nominated").sum()),
        n_win=("status", lambda s: (s == "winner").sum()),
    )
    .reset_index()
)

# Some clean up of the awards data
awards["year"] = awards["year"].astype(int)
awards.isbns = awards.isbns.str.replace("-", "")

# Split off openlibrary and isbn ids into separate dataframes
openlibrary_ids = explode_id_column(awards, "openlibrary_ids")
isbns = explode_id_column(awards, "isbns")

openlibrary_ids.to_csv(snakemake.output["openlibrary_ids"])
isbns.to_csv(snakemake.output["isbns"])

# Drop the original openlibrary and isbn columns
# and explode the author qids
awards = (
    awards.drop(columns=["openlibrary_ids", "isbns"])
    .explode("author_qids")
    .rename(columns={"author_qids": "author_qid"})
)


# Read cumulative awards
cumulative_awards = pd.read_csv(snakemake.input["wins_as_of"])

# Explode author QIDs and join with cumulative awards
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
    - pd.concat([openlibrary_ids.work_qid, isbns.work_qid]).nunique()
)
print(
    f"{n_missing}/{awards.work_qid.nunique()} works lack identifying information (openlibrary/isbn/etc.)"
)

awards.to_csv(snakemake.output["cleaned_novels"], index=False)
