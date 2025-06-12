import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OneHotEncoder


class RowCountEncoder(TransformerMixin, BaseEstimator):
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
