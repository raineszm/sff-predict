import pandas as pd

AUTHOR_COLUMNS = ["age", "gender", "awards_as_of_year", "birth_country"]


def one_hot_encoding(df: pd.DataFrame) -> pd.DataFrame:
    return (
        pd.get_dummies(df, columns=["gender"], prefix="gender", drop_first=True)
        .pipe(pd.get_dummies, columns=["birth_country"], prefix="birth_country", drop_first=True)
    )

def drop_nulls(df: pd.DataFrame) -> pd.DataFrame:
    return df.dropna(subset=["age", "topicality"])


def fix_nom_counts(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(n_nom_all=df["n_nom"] + df["n_win"]).drop(columns=["n_nom"])


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
            tot_cohort_nom=("n_nom_all", "sum"),
            tot_cohort_awards=("n_win", "sum"),
        )
        .reset_index(),
        on="year",
        how="outer",
    )


def parse_dates(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        month=pd.to_datetime(df["pubDate"], errors="coerce").dt.month
    ).drop(columns=["pubDate"])


def linprepare_df(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.pipe(one_hot_encoding)
        .pipe(drop_nulls)
        .pipe(fix_nom_counts)
        .pipe(prepare_bestseller_stats)
        .pipe(parse_dates)
        .pipe(compute_cohort_stats)
    )
