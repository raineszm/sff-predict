from sentence_transformers import SentenceTransformer
from sklearn.svm import LinearSVC
import pandas as pd
import numpy as np


class DebiasingModel:
    def __init__(self, model_name: str = "all-MiniLM-L6-v2"):
        self.model = SentenceTransformer(model_name)
        self.svc = LinearSVC(class_weight="balanced")
        self.nullspace_projector = None

    def fit(self, train_desc: pd.DataFrame):
        X_train = self.model.encode(train_desc["description"].tolist())
        self.svc.fit(X_train, train_desc["says_award"])
        w = self.svc.coef_ / np.linalg.norm(self.svc.coef_)
        self.nullspace_projector = np.eye(X_train.shape[1]) - np.outer(w, w)

    def transform(self, X: np.ndarray):
        return X @ self.nullspace_projector
