from sentence_transformers import SentenceTransformer
from sklearn.svm import LinearSVC
import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from typing import Optional


def learn_nullspace_normal_vector(
    train_desc: pd.DataFrame, model_name: str = "all-MiniLM-L6-v2"
):
    model = SentenceTransformer(model_name)
    X_train = model.encode(train_desc["description"].tolist())
    svc = LinearSVC(class_weight="balanced")
    svc.fit(X_train, train_desc["says_award"])
    w = svc.coef_ / np.linalg.norm(svc.coef_)
    return w


class NullspaceProjector(BaseEstimator, TransformerMixin):
    def __init__(self, normal_vector: np.ndarray):
        self.w = normal_vector
        self.nullspace_projector = np.eye(self.w.shape[0]) - np.outer(self.w, self.w)

    def fit(self, X: np.ndarray, y: Optional[np.ndarray] = None):
        if self.w.shape[1] != X.shape[1]:
            raise ValueError(
                "Feature dimension mismatch: normal vector has {:d} features, but data has {:d} features".format(
                    self.w.shape[0], X.shape[1]
                )
            )
        return self

    def transform(self, X: np.ndarray):
        return X @ self.nullspace_projector
