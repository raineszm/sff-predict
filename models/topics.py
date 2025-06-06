import numpy as np
from sentence_transformers import SentenceTransformer
from bertopic import BERTopic
from bertopic.representation import KeyBERTInspired
from sklearn.metrics.pairwise import pairwise_kernels
from scipy.stats import kurtosis
from statsmodels.tsa.stattools import adfuller
from abc import ABC, abstractmethod
import pandas as pd
from typing import Optional


def fill_missing_months(ts: pd.Series, months: pd.Series) -> pd.Series:
    ts.index = months
    return ts.reindex(range(1, 13), fill_value=0).sort_index()


MIN_MONTHLY_FREQUENCY = 10


class TopicFilter(ABC):
    @abstractmethod
    def is_topic_selected(self, topic_data: pd.DataFrame) -> bool:
        raise NotImplementedError

    def filter_topics(self, topics_over_time: pd.DataFrame) -> list:
        self.topics = []
        for topic, group in topics_over_time.groupby("Topic"):
            if all(group["Frequency"] < MIN_MONTHLY_FREQUENCY):
                continue
            if self.is_topic_selected(
                fill_missing_months(group["Frequency"], group["Timestamp"])
            ):
                self.topics.append(topic)
        return self.topics

    def fit(self, topics_over_time: pd.DataFrame) -> "TopicFilter":
        self.topics = self.filter_topics(topics_over_time).tolist()
        return self


class HighKurtosisFilter(TopicFilter):
    def __init__(self, kurtosis_threshold: float = 2, **kwargs):
        super().__init__(**kwargs)
        self.kurtosis_threshold = kurtosis_threshold

    def is_topic_selected(self, ts: pd.Series) -> bool:
        return kurtosis(ts, bias=False, fisher=True) > self.kurtosis_threshold


class NonStationaryFilter(TopicFilter):
    def __init__(self, pvalue_threshold: float = 0.05, n_periods: int = 3, **kwargs):
        super().__init__(**kwargs)
        self.pvalue_threshold = pvalue_threshold
        self.n_periods = n_periods

    def is_topic_selected(self, ts: pd.Series) -> bool:
        periodic = np.tile(ts.sort_index(), self.n_periods * len(ts))
        return adfuller(periodic, regression="ctt")[1] < self.pvalue_threshold


class CompoundFilter(TopicFilter):
    def __init__(self, filters: list[TopicFilter]):
        self.filters = filters

    def is_topic_selected(self, ts: pd.Series) -> bool:
        return all(filter.is_topic_selected(ts) for filter in self.filters)


class TopicExtractionModel:
    def __init__(
        self,
        embedding_model: SentenceTransformer,
        topic_filter: Optional[TopicFilter] = None,
    ):
        self.topic_filter = topic_filter or CompoundFilter(
            [
                HighKurtosisFilter(kurtosis_threshold=2),
                NonStationaryFilter(pvalue_threshold=0.1, n_periods=3),
            ]
        )
        self.model = BERTopic(
            language="english",
            representation_model=KeyBERTInspired(),
            embedding_model=embedding_model,
            calculate_probabilities=False,
        )

    def _recompute_centroids(self):
        hot_topics = self.topic_filter.filter_topics(self.topics_over_time)

        self.centroids = self.model.topic_embeddings_[hot_topics]
        topic_sizes = np.array([self.model.topic_sizes_[i] for i in hot_topics])
        self.weights = topic_sizes / topic_sizes.sum()

    def fit(self, texts: list[str], embeddings: np.ndarray, months: list[int]):
        self.model.fit(texts, embeddings)

        self.topics_over_time = self.model.topics_over_time(
            texts,
            months,
            global_tuning=False,
            evolution_tuning=False,
        )
        self._recompute_centroids()
        return self

    def compute_alignment(
        self,
        embeddings: np.ndarray,
        weighted: bool = False,
        metric: str = "cosine",
        p: int = 1,
        mod_sim: bool = False,
        **metric_kwargs,
    ) -> np.array:
        similarities = pairwise_kernels(
            embeddings, self.centroids, metric=metric, **metric_kwargs
        )
        if mod_sim:
            similarities = np.abs(similarities)

        if p == np.inf:
            return np.max(similarities, axis=1)

        similarities = similarities**p
        if weighted:
            return (similarities.dot(self.weights)) ** (1 / p)
        else:
            return (similarities.sum(axis=1)) ** (1 / p)
