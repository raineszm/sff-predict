from snakemake.script import snakemake
import pandas as pd
from sentence_transformers import SentenceTransformer

# Initialize the model
model = SentenceTransformer(snakemake.params.embedding_model)

print("Embedding descriptions")

# Read the descriptions
df = pd.read_csv(snakemake.input["descriptions"], index_col="work_qid")

# Get embeddings
embeddings = model.encode(
    df["description"].tolist(), convert_to_numpy=True, show_progress_bar=True
).astype("float32")

# Create output dataframe
result_df = pd.DataFrame({"embedding": list(embeddings)}, index=df.index)

# Save to parquet
result_df.to_parquet(
    snakemake.output["description_embeddings"], compression="zstd", compression_level=12
)
