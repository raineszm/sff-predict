import duckdb
import duckdb.typing
import pyarrow as pa
import numpy as np
from sentence_transformers import SentenceTransformer

# A good all-rounder model
model = SentenceTransformer("all-MiniLM-L6-v2")


def embed_text(input_text: pa.Array) -> pa.ListArray:
    texts = input_text.to_pylist()
    embeddings = model.encode(
        texts, convert_to_numpy=True, show_progress_bar=True
    ).astype(np.float32)
    D = embeddings.shape[1]
    N = embeddings.shape[0]
    offsets = pa.array([i * D for i in range(N + 1)], type=pa.int32())
    values = pa.array(embeddings.reshape(-1), type=pa.float32())
    return pa.ListArray.from_arrays(offsets, values)


def register_embedding_function():
    # Register embedding function with duckdb
    duckdb.create_function(
        "embed_text",
        embed_text,
        parameters=[duckdb.typing.VARCHAR],
        return_type=list[duckdb.typing.FLOAT],
        type="arrow",
    )
