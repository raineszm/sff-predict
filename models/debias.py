"""
Bias mitigation for text embeddings.
"""

from sentence_transformers import SentenceTransformer
from sklearn.svm import LinearSVC
import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from typing import Optional


def learn_nullspace_normal_vector(
    train_desc: pd.DataFrame, model_name: str = "all-MiniLM-L6-v2"
):
    """
    Learn a normal vector to the bias subspace using SVM.

    Parameters
    ----------
    train_desc : pd.DataFrame
        Training data with columns:
        - 'description': Text descriptions to embed
        - 'says_award': Binary labels indicating bias presence
    model_name : str, default="all-MiniLM-L6-v2"
        Name of the SentenceTransformer model to use for embeddings.

    Returns
    -------
    np.ndarray
        Normalized normal vector to the bias subspace.
    """
    model = SentenceTransformer(model_name)
    X_train = model.encode(train_desc["description"].tolist())
    svc = LinearSVC(class_weight="balanced")
    svc.fit(X_train, train_desc["says_award"])
    w = svc.coef_ / np.linalg.norm(svc.coef_)
    return w


class NullspaceProjector(BaseEstimator, TransformerMixin):
    """
    Project embeddings to remove bias using nullspace projection.

    This transformer removes bias from embeddings by projecting them
    onto the nullspace of the bias direction.

    Parameters
    ----------
    normal_vector : np.ndarray
        Normal vector to the bias subspace (learned from training data).

    Attributes
    ----------
    w : np.ndarray
        Normal vector to bias subspace.
    nullspace_projector : np.ndarray
        Projection matrix for removing bias.

    """

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
