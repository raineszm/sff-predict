import duckdb
import sys

sys.path.append(".")

from utils.embedding import register_embedding_function
from utils.sentiment import register_sentiment_function
from tqdm.auto import tqdm


def embed_headlines(input_headlines, output_headline_embeddings):
    print("Embedding headlines (This may take a while)")

    n_headlines = duckdb.sql(f"SELECT COUNT(*) FROM '{input_headlines}'").fetchone()[0]

    with tqdm(total=n_headlines, desc="Embedding headlines") as total_pbar:
        register_sentiment_function()
        register_embedding_function(total_pbar=total_pbar)

        # Strip the nyt://article/ prefix from the id
        # and combine the headline and abstract for embedding
        duckdb.sql(f"""
        CREATE VIEW concat_headline_abstract AS (
            SELECT id.substr(strlen('nyt://article/') + 1) as id, headline || ' ' || abstract AS headline_abstract
            FROM '{input_headlines}'
        );
        """)

        # Do vectorized semanticembedding and sentiment analysis
        # of the headlines
        duckdb.sql(f"""
        COPY (
        SELECT id, embed_text(headline_abstract) AS embedding, sentiment(headline_abstract) AS sentiment
        FROM concat_headline_abstract
        )
        TO '{output_headline_embeddings}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """)


if __name__ == "__main__":
    from snakemake.script import snakemake

    embed_headlines(
        snakemake.input["headlines"], snakemake.output["headline_embeddings"]
    )
