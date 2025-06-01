from snakemake.script import snakemake
import duckdb
from utils.embedding import register_embedding_function

register_embedding_function()

print("Embedding headlines (This may take a while)")

duckdb.sql(f"""
COPY (
SELECT id, embed_text(headline || ' ' || abstract) AS embedding
FROM '{snakemake.input["headlines"]}'
)
TO '{snakemake.output["headline_embeddings"]}' (FORMAT PARQUET, COMPRESSION ZSTD)
""")
