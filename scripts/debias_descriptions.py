import pandas as pd
import numpy as np
from snakemake.script import snakemake
from models.debias import NullspaceProjector, learn_nullspace_normal_vector

print("Training Debiasing Model...")
train_desc = pd.read_csv(snakemake.input.train_desc, index_col="work_qid")
w = learn_nullspace_normal_vector(train_desc, model_name=snakemake.params.model_name)
debiasing_model = NullspaceProjector(w)

print("Debiasing descriptions")
description_embeddings = pd.read_parquet(snakemake.input.description_embeddings)
X_data = np.vstack(description_embeddings["embedding"].tolist())
X_debiased = debiasing_model.fit_transform(X_data)

print("Saving debiased descriptions")
df_debiased = pd.DataFrame(
    {"embedding": list(X_debiased)}, index=description_embeddings.index
)
df_debiased.to_parquet(
    snakemake.output.descriptions_debiased, compression="zstd", compression_level=15
)
