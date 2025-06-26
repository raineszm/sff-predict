"""
Baseline model for award winner prediction.

This module provides a simple baseline model that predicts award winners
based on nomination counts using multinomial sampling.
"""

import pandas as pd
import pandas.api.typing as pdt
import numpy as np


def naive_cohort_win_counts(df_year: pdt.DataFrameGroupBy) -> np.ndarray:
    """
    Generate random winner predictions based on nomination counts for a year cohort.

    For a given year, this function predicts winners by sampling from a
    multinomial distribution where the probability of winning is proportional
    to the number of nominations each work has received.

    Parameters
    ----------
    df_year : pdt.DataFrameGroupBy
        DataFrame grouped by year containing 'n_nom_all' and 'n_win' columns.

    Returns
    -------
    np.ndarray
        Array of predicted winner counts (0 or 1) for each nominee in the cohort.


    """
    tot_cohort_nom = df_year["n_nom_all"].sum()
    p = df_year["n_nom_all"] / tot_cohort_nom
    k = df_year["n_win"].sum()
    return np.random.multinomial(k, p)


def naive_win_counts(X: pd.DataFrame) -> pd.Series:
    """
    Predict winners using multinomial sampling proportional to nominations.

    Parameters
    ----------
    X : pd.DataFrame
        DataFrame containing 'year', 'n_nom_all', and 'n_win' columns.
        The 'n_win' column should contain the actual number of winners
        for each year (used to determine how many winners to predict).

    Returns
    -------
    pd.Series
        Series of predicted winner counts (0 or 1) for each row, indexed
        to match the input DataFrame.

    """
    return (
        X.groupby("year").apply(naive_cohort_win_counts, include_groups=False).explode()
    )
