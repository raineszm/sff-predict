"""
Data preparation utilities for science fiction and fantasy award prediction.

This module provides functions for cleaning and transforming raw award data
into a format suitable for machine learning models.
"""

import pandas as pd

AUTHOR_COLUMNS = ["age", "gender", "awards_as_of_year", "birth_country"]
"""Columns related to author information that need to be aggregated."""


def reduce_authors(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate author information for works with multiple authors.

    For works with multiple authors, this function aggregates author-related
    statistics by taking means, lists, or sums as appropriate. This prevents
    data duplication while preserving author information.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing author data with columns from AUTHOR_COLUMNS
        and a 'work_qid' column to identify unique works.

    Returns
    -------
    pd.DataFrame
        DataFrame with aggregated author statistics per work:
        - mean_age: Average age of all authors
        - gender: List of all author genders
        - birth_country: List of all author birth countries
        - awards_as_of_year: Sum of all author awards

    """
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


def fix_nom_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Combine nomination and win counts into total nominations.

    Creates a new column 'n_nom_all' that represents the total number of
    nominations (including wins) and removes the original 'n_nom' column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'n_nom' and 'n_win' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with 'n_nom_all' column (n_nom + n_win) and 'n_nom' removed.

    Example
    -------
    >>> df = pd.DataFrame({'n_nom': [3, 1], 'n_win': [1, 0]})
    >>> result = fix_nom_counts(df)
    >>> print(result['n_nom_all'].tolist())
    [4, 1]
    """
    return df.assign(n_nom_all=df["n_nom"] + df["n_win"]).drop(columns=["n_nom"])


def prepare_bestseller_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and rename bestseller statistics columns.

    Fills missing values with 0 and renames columns to be more descriptive.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing bestseller columns: 'max', 'count', 'median'.

    Returns
    -------
    pd.DataFrame
        DataFrame with cleaned bestseller columns:
        - max_bestseller_rank: Highest bestseller rank achieved
        - months_on_bestseller: Number of months on bestseller list
        - median_bestseller_rank: Median bestseller rank

    """
    return df.fillna({"max": 0, "count": 0, "median": 0}).rename(
        columns={
            "max": "max_bestseller_rank",
            "count": "months_on_bestseller",
            "median": "median_bestseller_rank",
        }
    )


def parse_dates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract month from publication dates.

    Converts the 'pubDate' column to datetime and extracts the month (1-12),
    then removes the original 'pubDate' column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a 'pubDate' column with date strings.

    Returns
    -------
    pd.DataFrame
        DataFrame with 'month' column (1-12) and 'pubDate' removed.

    """
    return df.assign(month=pd.to_datetime(df.pubDate).dt.month).drop(
        columns=["pubDate"]
    )


def prepare_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply complete data preparation pipeline.

    Parameters
    ----------
    df : pd.DataFrame
        Raw input DataFrame with award and author data.

    Returns
    -------
    pd.DataFrame
        Fully prepared DataFrame ready for machine learning modeling.

    """
    return (
        df.pipe(reduce_authors)
        .pipe(fix_nom_counts)
        .pipe(prepare_bestseller_stats)
        .pipe(parse_dates)
    )
