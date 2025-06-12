import pandas as pd
import pandas.api.typing as pdt
import numpy as np


def naive_cohort_win_counts(df_year: pdt.DataFrameGroupBy) -> np.ndarray:
    p = df_year["n_nom_all"] / df_year["tot_cohort_nom"]
    k = df_year["tot_cohort_awards"].iat[0]
    return np.random.multinomial(k, p)


def naive_win_counts(X: pd.DataFrame) -> pd.Series:
    """
    Model the winners by doing a multinomial draw from the nominees
    with a weight proportional to the number of nominations
    """
    return (
        X.groupby("year").apply(naive_cohort_win_counts, include_groups=False).explode()
    )
