from snakemake.script import snakemake
import duckdb
from utils.embedding import register_embedding_function
from tqdm.auto import tqdm

print("Embedding headlines (This may take a while)")

n_headlines = duckdb.sql(
    f"SELECT COUNT(*) FROM '{snakemake.input['headlines']}'"
).fetchone()[0]

with tqdm(total=n_headlines, desc="Embedding headlines") as total_pbar:
    register_embedding_function(total_pbar=total_pbar)

    duckdb.sql(f"""
    COPY (
    SELECT id, embed_text(headline || ' ' || abstract) AS embedding
    FROM '{snakemake.input["headlines"]}'
    )
    TO '{snakemake.output["headline_embeddings"]}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
