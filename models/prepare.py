import pandas as pd

AUTHOR_COLUMNS = ["age", "gender", "awards_as_of_year", "birth_country"]


def reduce_authors(df: pd.DataFrame) -> pd.DataFrame:
    author_data = (
        df.groupby("work_qid")[AUTHOR_COLUMNS]
        .agg(
            mean_age=("age", "mean"),
            gender=("gender", list),
            birth_country=("birth_country", list),
            awards_as_of_year=("awards_as_of_year", "sum"),
        )
        .reset_index()
    )

    non_author_data = df.drop(columns=AUTHOR_COLUMNS).drop_duplicates(
        subset=["work_qid"], keep="first"
    )

    return pd.merge(non_author_data, author_data, on="work_qid", how="left")


def prepare_bestseller_stats(df: pd.DataFrame) -> pd.DataFrame:
    return df.fillna({"max": 0, "count": 0, "median": 0}).rename(
        columns={
            "max": "max_bestseller_rank",
            "count": "months_on_bestseller",
            "median": "median_bestseller_rank",
        }
    )


def compute_cohort_stats(df: pd.DataFrame) -> pd.DataFrame:
    return df.merge(
        df.groupby("year")
        .agg(
            tot_cohort_nom=("n_nom", "sum"),
            tot_cohort_awards=("n_win", "sum"),
        )
        .reset_index(),
        on="year",
        how="outer",
    )


def estimate_target(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign().drop(columns=["n_win"])


def parse_dates(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(month=pd.to_datetime(df.pubDate).dt.month).drop(
        columns=["pubDate"]
    )


def prepare_df(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.pipe(reduce_authors)
        .pipe(prepare_bestseller_stats)
        .pipe(parse_dates)
        .pipe(compute_cohort_stats)
    )
