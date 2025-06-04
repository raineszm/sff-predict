import pandas as pd
import numpy as np
from snakemake.script import snakemake
from models.debias import DebiasingModel

print("Training Debiasing Model...")
debiasing_model = DebiasingModel(model_name="all-MiniLM-L6-v2")
train_desc = pd.read_csv(snakemake.input.train_desc, index_col="work_qid")
debiasing_model.fit(train_desc)

print("Debiasing descriptions")
description_embeddings = pd.read_parquet(snakemake.input.description_embeddings)
X_data = np.vstack(description_embeddings["embedding"].values)
X_debiased = debiasing_model.transform(X_data)

print("Saving debiased descriptions")
df_debiased = pd.DataFrame(
    {"embedding": list(X_debiased)}, index=description_embeddings.index
)
df_debiased.to_parquet(
    snakemake.output.descriptions_debiased, compression="zstd", compression_level=15
)
