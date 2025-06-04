import pandas as pd
from sentence_transformers import SentenceTransformer
from snakemake.script import snakemake

# Initialize the models
embedding_model = SentenceTransformer(snakemake.params.embedding_model)


def embed_headlines(input_headlines, output_headline_embeddings):
    print("Embedding headlines (This may take a while)")

    # Read the input headlines
    df = pd.read_parquet(input_headlines)

    # Combine headline and abstract
    df["headline_abstract"] = (
        df["headline"].str.removeprefix("LEAD:").str.lower() + " " + df["abstract"]
    )

    # Get embeddings
    embeddings = embedding_model.encode(
        df["headline_abstract"].tolist(), convert_to_numpy=True, show_progress_bar=True
    ).astype("float32")

    # Create output dataframe
    result_df = pd.DataFrame(
        {
            "id": df["id"],
            "embedding": list(embeddings),
        }
    )

    # Save to parquet
    result_df.to_parquet(
        output_headline_embeddings, compression="zstd", compression_level=12
    )


if __name__ == "__main__":
    embed_headlines(
        snakemake.input["headlines"], snakemake.output["headline_embeddings"]
    )
