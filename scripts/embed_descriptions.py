from snakemake.script import snakemake
import duckdb
from utils.embedding import register_embedding_function

register_embedding_function()

print("Embedding descriptions")

duckdb.sql(f"""
COPY (
SELECT work_qid, embed_text(description) AS embedding
FROM '{snakemake.input["descriptions"]}'
)
TO '{snakemake.output["description_embeddings"]}' (FORMAT PARQUET, COMPRESSION ZSTD)
""")
