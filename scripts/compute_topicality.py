from snakemake.script import snakemake
import pandas as pd
import numpy as np
from models.topics import (
    CompoundFilter,
    HighKurtosisFilter,
    NonStationaryFilter,
    TopicExtractionModel,
)
from sentence_transformers import SentenceTransformer
from tqdm.auto import tqdm

# Fetch debiased embeddings
desc_embeddings = pd.read_parquet(snakemake.input.debiased_embeddings)
# Fetch novel data
novels = pd.read_csv(snakemake.input.novels)

# Combine so we can group by year
desc_merged = desc_embeddings.merge(
    novels[["work_qid", "year"]].groupby("work_qid").agg({"year": "first"}),
    on="work_qid",
    how="left",
)

desc_merged = desc_merged[desc_merged.year >= 1959]


# Load headline embeddings
headline_embeddings = pd.read_parquet(snakemake.input.headline_embeddings)
# Load headline text
headlines = pd.read_json(snakemake.input.headlines, lines=True)
headlines_merged = headlines.merge(headline_embeddings, on="id", how="left")


def compute_topicality(year: int, df_year: pd.DataFrame) -> np.array:
    # load and train the topic model
    yearly_headlines = headlines_merged[headlines_merged.year == year]

    tm = TopicExtractionModel(
        embedding_model=SentenceTransformer("all-MiniLM-L6-v2"),
        topic_filter=CompoundFilter(
            [
                HighKurtosisFilter(kurtosis_threshold=2),
                NonStationaryFilter(pvalue_threshold=0.1, n_periods=3),
            ]
        ),
    )

    tm.fit(
        yearly_headlines.headline.tolist(),
        np.stack(yearly_headlines.embedding.values),
        yearly_headlines.month.tolist(),
    )

    return pd.DataFrame({"work_qid": df_year.index, "topicality": df_year.topicality})


for year, df_year in tqdm(desc_merged.groupby("year"), desc="Computing topicality"):
    topicality_scores = compute_topicality(year, df_year).reset_index(drop=True)
    topicality_scores.to_csv(snakemake.output.topicality_scores, index=False, mode="a")
