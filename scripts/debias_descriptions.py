from sklearn.svm import LinearSVC
import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
from snakemake.script import snakemake


def train_nullspace_projector(train_desc: pd.DataFrame):
    """
    Trains a projector for LEACE debiasing.

    Since some descriptions include mention of the awards that the novel has won, we want to remove this information from the descriptions.

    We do this by training a linear SVM to classify whether a description mentions an award, and then projecting the descriptions onto the nullspace of the SVM's weight vector.
    """
    model = SentenceTransformer("all-MiniLM-L6-v2")
    svc = LinearSVC(class_weight="balanced")
    X_train = model.encode(train_desc["description"].tolist())
    svc.fit(X_train, train_desc["says_award"])
    w = svc.coef_ / np.linalg.norm(svc.coef_)
    return np.eye(X_train.shape[1]) - np.outer(w, w)


print("Training nullspace projector...")
train_desc = pd.read_csv(snakemake.input.train_desc, index_col="work_qid")
proj = train_nullspace_projector(train_desc)

print("Debiasing descriptions")
description_embeddings = pd.read_parquet(snakemake.input.description_embeddings)
X_data = np.vstack(description_embeddings["embedding"].values)
X_debiased = np.dot(X_data, proj)

print("Saving debiased descriptions")
df_debiased = pd.DataFrame(
    {"embedding": list(X_debiased)}, index=description_embeddings.index
)
df_debiased.to_parquet(
    snakemake.output.descriptions_debiased, compression="zstd", compression_level=15
)
