import numpy as np
from sentence_transformers import SentenceTransformer
from bertopic import BERTopic
from bertopic.representation import KeyBERTInspired
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import kurtosis


class TopicExtractionModel:
    def __init__(
        self,
        embedding_model: SentenceTransformer,
        monthly_threshold_frequency: int = 10,
        kurtosis_threshold: float = 2,
    ):
        self.monthly_threshold_frequency = monthly_threshold_frequency
        self.kurtosis_threshold = kurtosis_threshold
        self.model = BERTopic(
            language="english",
            representation_model=KeyBERTInspired(),
            embedding_model=embedding_model,
            calculate_probabilities=False,
        )

    def fit(self, texts: list[str], embeddings: np.ndarray, months: list[int]):
        self.model.fit(texts, embeddings)

        self.topics_over_time = self.model.topics_over_time(
            texts,
            months,
            global_tuning=False,
            evolution_tuning=False,
        )

        kt = (
            self.topics_over_time[
                self.topics_over_time["Frequency"] >= self.monthly_threshold_frequency
            ]
            .sort_values("Timestamp")
            .groupby("Topic")
            .Frequency.agg(lambda x: kurtosis(x, bias=False, fisher=True))
        )
        hot_topics = kt[kt > self.kurtosis_threshold].index

        self.centroids = self.model.topic_embeddings_[hot_topics]
        topic_sizes = np.array([self.model.topic_sizes_[i] for i in hot_topics])
        self.weights = topic_sizes / topic_sizes.sum()
        return self

    def compute_alignment(self, embeddings: np.ndarray) -> np.array:
        similarities = cosine_similarity(embeddings, self.centroids)
        return (1 + similarities).sum(axis=1)
