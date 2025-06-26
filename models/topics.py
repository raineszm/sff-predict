"""
Topic modeling and topicality scoring using BERTopic.
"""

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
    """
    Fill months with missing counts with 0s.


    Parameters
    ----------
    ts : pd.Series
        Time series data with frequency values.
    months : pd.Series
        Month indices corresponding to the time series data.

    Returns
    -------
    pd.Series
        Complete time series with all months (1-12) and missing values filled with 0.

    Example
    -------
    >>> ts = pd.Series([10, 15, 8], index=[1, 3, 5])
    >>> months = pd.Series([1, 3, 5])
    >>> result = fill_missing_months(ts, months)
    >>> print(result[2])  # Missing month filled with 0
    0
    """
    ts.index = months
    return ts.reindex(range(1, 13), fill_value=0).sort_index()


MIN_MONTHLY_FREQUENCY = 10
"""Minimum monthly frequency threshold for topic filtering."""


class TopicFilter(ABC):
    """
    Abstract base class for topic filtering strategies.
    """

    @abstractmethod
    def is_topic_selected(self, topic_data: pd.DataFrame) -> bool:
        """
        Determine if a topic should be selected based on filtering criteria.

        Parameters
        ----------
        topic_data : pd.DataFrame
            Topic frequency data over time.

        Returns
        -------
        bool
            True if the topic should be selected, False otherwise.
        """
        raise NotImplementedError

    def filter_topics(self, topics_over_time: pd.DataFrame) -> list:
        """
        Filter topics based on the implemented selection criteria.

        Parameters
        ----------
        topics_over_time : pd.DataFrame
            DataFrame with columns 'Topic', 'Frequency', and 'Timestamp'.

        Returns
        -------
        list
            List of topic IDs that pass the filtering criteria.
        """
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
        """
        Fit the filter to the topics over time data.

        Parameters
        ----------
        topics_over_time : pd.DataFrame
            Topics over time data to fit the filter to.

        Returns
        -------
        TopicFilter
            Self for method chaining.
        """
        self.topics = self.filter_topics(topics_over_time).tolist()
        return self


class HighKurtosisFilter(TopicFilter):
    """
    Filter topics based on kurtosis of frequency distribution.

    This filter selects topics that have high kurtosis, indicating they
    have peaked frequency distributions rather than uniform ones.

    Parameters
    ----------
    kurtosis_threshold : float, default=2
        Minimum kurtosis value required for topic selection.
    """

    def __init__(self, kurtosis_threshold: float = 2, **kwargs):
        super().__init__(**kwargs)
        self.kurtosis_threshold = kurtosis_threshold

    def is_topic_selected(self, ts: pd.Series) -> bool:
        return kurtosis(ts, bias=False, fisher=True) > self.kurtosis_threshold


class NonStationaryFilter(TopicFilter):
    """
    Filter topics based on stationarity tests.

    This filter selects topics that show non-stationary behavior,
    indicating they have temporal trends rather than random fluctuations.

    Parameters
    ----------
    pvalue_threshold : float, default=0.05
        Maximum p-value for stationarity test.
    n_periods : int, default=3
        Number of periods to repeat for testing.
    """

    def __init__(self, pvalue_threshold: float = 0.05, n_periods: int = 3, **kwargs):
        super().__init__(**kwargs)
        self.pvalue_threshold = pvalue_threshold
        self.n_periods = n_periods

    def is_topic_selected(self, ts: pd.Series) -> bool:
        periodic = np.tile(ts.sort_index(), self.n_periods * len(ts))
        return adfuller(periodic, regression="ctt")[1] < self.pvalue_threshold


class CompoundFilter(TopicFilter):
    """
    Combine multiple topic filters with AND logic.

    A topic is selected only if it passes all the individual filters.

    Parameters
    ----------
    filters : list[TopicFilter]
        List of topic filters to combine.
    """

    def __init__(self, filters: list[TopicFilter]):
        self.filters = filters

    def is_topic_selected(self, ts: pd.Series) -> bool:
        return all(filter.is_topic_selected(ts) for filter in self.filters)


class TopicExtractionModel:
    """
    Modeling class using BERTopic for topicality analysis.

    Parameters
    ----------
    embedding_model : SentenceTransformer
        Pre-trained sentence transformer model for text embeddings.
    topic_filter : TopicFilter, optional
        Filter for selecting "hot" topics. Defaults to a compound filter
        combining high kurtosis and non-stationary criteria.
    """

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
        """
        Recompute topic centroids and weights based on filtered topics.

        Updates the centroids and weights attributes based on the
        currently selected "hot" topics.
        """
        hot_topics = self.topic_filter.filter_topics(self.topics_over_time)

        self.centroids = self.model.topic_embeddings_[hot_topics]
        topic_sizes = np.array([self.model.topic_sizes_[i] for i in hot_topics])
        self.weights = topic_sizes / topic_sizes.sum()

    def fit(self, texts: list[str], embeddings: np.ndarray, months: list[int]):
        """
        Fit the topic model to training data.

        Parameters
        ----------
        texts : list[str]
            List of text documents (e.g., news headlines).
        embeddings : np.ndarray
            Pre-computed embeddings for the texts.
        months : list[int]
            Month indices (1-12) corresponding to each text.

        Returns
        -------
        TopicExtractionModel
            Self for method chaining.
        """
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
        """
        Compute topicality alignment scores for new embeddings.

        Calculates how well each embedding aligns with the "hot" topics
        identified during training.

        Parameters
        ----------
        embeddings : np.ndarray
            Embeddings to compute alignment for (e.g., book descriptions).
        weighted : bool, default=False
            Whether to weight topic similarities by topic size.
        metric : str, default="cosine"
            Similarity metric to use ('cosine', 'euclidean', etc.).
        p : int, default=1
            Power for L-p norm aggregation (use np.inf for max).
        mod_sim : bool, default=False
            Whether to use absolute similarity values.
        **metric_kwargs
            Additional arguments for the similarity metric.

        Returns
        -------
        np.array
            Array of alignment scores for each embedding.

        Example
        -------
        >>> model = TopicExtractionModel(embedding_model)
        >>> model.fit(headlines, headline_embeddings, headline_months)
        >>> book_alignments = model.compute_alignment(book_embeddings)
        >>> print(book_alignments.mean())  # Average topicality score
        0.45
        """
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
