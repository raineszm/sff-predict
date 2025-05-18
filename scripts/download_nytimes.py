import httpx
import os
import json
from sentence_transformers import SentenceTransformer

API_BASE = "https://api.nytimes.com/svc/archive/v1/"


class NYTimesArchiveAPI:
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.client = httpx.Client(base_url=API_BASE)

    def get_articles(self, year: int, month: int):
        url = f"{year}/{month}.json"
        params = {"api-key": self.api_key}
        response = self.client.get(url, params=params)
        response.raise_for_status()
        return response.json()

    @classmethod
    def from_envfile(cls):
        api_key = os.getenv("NYTIMES_API_KEY")
        if not api_key:
            raise ValueError("NYTIMES_API_KEY not found in environment variables.")
        return cls(api_key)


api = NYTimesArchiveAPI.from_envfile()

with open("data/raw/nytimes_articles.json", "w") as f:
    # model = SentenceTransformer("all-MiniLM-L6-v2")
    # headlines = [
    #     article["headline"]["main"]
    #     for article in api.get_articles(1959, 1)["response"]["docs"]
    # ]
    # embeddings = model.encode(headlines)
    # json.dump(embeddings.mean(axis=0).tolist(), f)
    json.dump(api.get_articles(1959, 1)["response"]["docs"], f)
