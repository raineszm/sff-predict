import duckdb
import pyarrow as pa
from transformers import pipeline
from typing import Optional
from tqdm.auto import tqdm


def assign_sentiment(
    text: pa.Array, classifier: pipeline, total_pbar: Optional[tqdm] = None
) -> pa.Array:
    sentiments = [t["label"] for t in classifier(text.to_pylist())]
    if total_pbar is not None:
        total_pbar.update(len(sentiments))
    return pa.array(sentiments, type=pa.string())


def register_sentiment_function(total_pbar: Optional[tqdm] = None):
    # We use distilbert with transformers to get a sentiment score
    classifier = pipeline(
        "sentiment-analysis",
        model="distilbert-base-uncased-finetuned-sst-2-english",
        revision="714eb0f",
    )

    duckdb.create_function(
        "sentiment",
        lambda x: assign_sentiment(x, classifier, total_pbar=total_pbar),
        parameters=[duckdb.typing.VARCHAR],
        return_type=duckdb.typing.VARCHAR,
        type="arrow",
    )
