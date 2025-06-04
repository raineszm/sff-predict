import numpy as np
from sentence_transformers import SentenceTransformer
from bertopic import BERTopic
from bertopic.representation import KeyBERTInspired
from sklearn.metrics.pairwise import cosine_similarity


class TopicExtractionModel:
    def __init__(self, embedding_model: SentenceTransformer, top_n_topics: int = 20):
        self.model = BERTopic(
            language="english",
            representation_model=KeyBERTInspired(),
            embedding_model=embedding_model,
            calculate_probabilities=False,
        )
        self.top_n_topics = top_n_topics
        self.top_n_topic_embeddings = None

    def fit(self, texts: list[str], embeddings: np.ndarray):
        self.model.fit(texts, embeddings)
        self.centroids = self.model.topic_embeddings_[1 : self.top_n_topics + 1, :]
        topic_sizes = np.array(
            [self.model.topic_sizes_[i] for i in range(self.top_n_topics)]
        )
        self.weights = topic_sizes / topic_sizes.sum()
        return self

    def compute_alignment(self, embeddings: np.ndarray) -> np.array:
        return cosine_similarity(embeddings, self.centroids).dot(self.weights)
