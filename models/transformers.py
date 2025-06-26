"""
Custom scikit-learn transformers for feature engineering.
"""

import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OneHotEncoder


class RowCountEncoder(TransformerMixin, BaseEstimator):
    """
    Encode list columns by exploding them, one hot encoding, and aggregating counts.

    Parameters
    ----------
    columns : list
        List of column names to encode. These columns should contain lists.
    key_column : str
        Column name to use for grouping after exploding the list columns.

    Attributes
    ----------
    encoder : OneHotEncoder
        Fitted one-hot encoder for the exploded values.
    other_columns : list
        Columns that are not being encoded.
    n_features_in_ : int
        Number of features seen during fit.

    """

    def __init__(self, columns, key_column):
        self.key_column = key_column
        self.columns = columns
        self.encoder = OneHotEncoder(
            sparse_output=False,
            handle_unknown="infrequent_if_exist",
            min_frequency=0.01,
        )
        self.encoder.set_output(transform="pandas")

    def _explode(self, X: pd.DataFrame) -> pd.DataFrame:
        return X[[self.key_column] + self.columns].explode(self.columns)

    def fit(self, X, y=None):
        self.n_features_in_ = len(X.columns)
        self.other_columns = [c for c in X.columns if c not in self.columns]
        self.encoder.fit(self._explode(X).drop(columns=[self.key_column]))
        return self

    def get_feature_names_out(self, input_features=None):
        return self.other_columns + self.encoder.get_feature_names_out(input_features)

    def transform(self, X: pd.DataFrame, y=None):
        exploded = self._explode(X)

        return pd.merge(
            X.drop(columns=self.columns),
            pd.concat(
                [
                    exploded[self.key_column],
                    self.encoder.transform(exploded.drop(columns=[self.key_column])),
                ]
            )
            .groupby(self.key_column)
            .sum(),
            on=self.key_column,
        )


def impute_topicality(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fill missing topicality scores with year-specific means.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'topicality' and 'year' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with imputed topicality scores.

    """
    return df.assign(
        topicality=df.groupby("year")["topicality"].transform(
            lambda x: x.fillna(x.mean())
        )
    )


def cohort_zscore(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """
    Compute z-scores within year cohorts.

    For each specified column, computes z-scores relative to the mean
    and standard deviation within the same year.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with 'year' column and columns to z-score.
    columns : list[str]
        List of column names to compute z-scores for.

    Returns
    -------
    pd.DataFrame
        DataFrame with z-scored columns (original columns replaced).
    """
    return df.assign(
        **{
            column: (df[column] - df.groupby("year")[column].transform("mean"))
            / df.groupby("year")[column].transform("std")
            for column in columns
        }
    )


def cohort_proportion(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """
    Compute proportions within year cohorts.

    For each specified column, computes the proportion of the total
    within each year.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with 'year' column and columns to compute proportions for.
    columns : list[str]
        List of column names to compute proportions for.

    Returns
    -------
    pd.DataFrame
        DataFrame with proportion columns (original columns replaced).
    """
    return df.assign(
        **{
            column: df[column] / df.groupby("year")[column].transform("sum")
            for column in columns
        }
    )
