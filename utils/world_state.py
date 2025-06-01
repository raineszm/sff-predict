from .api import CachedApi, path_to_cache_key
import os
from typing import Iterator, TypedDict
from httpx_ratelimiter import LimiterTransport


class NYTArticle(TypedDict):
    id: str
    headline: str
    abstract: str
    year: int
    month: int


class NYTimesArchiveAPI:
    API_BASE = "https://api.nytimes.com/svc/archive/v1/"

    def __init__(self, api_key: str):
        self.api_key = api_key
        self._api = CachedApi(
            "nytimes",
            self.API_BASE,
            key_fn=path_to_cache_key,
            base_transport=LimiterTransport(
                per_minute=5,
                per_day=500,
                max_delay=60_000,
            ),
        )

    def get_articles(self, year: int, month: int) -> dict:
        return self._api.call_endpoint(
            f"{year}/{month}.json", {"api-key": self.api_key}
        )

    def get_headlines(self, year: int, month: int) -> Iterator[NYTArticle]:
        articles = self.get_articles(year, month)
        for d in articles["response"]["docs"]:
            if d.get("print_page") == "1":
                yield {
                    "id": d["_id"],
                    "headline": d["headline"]["print_headline"],
                    "abstract": d["abstract"],
                    "year": year,
                    "month": month,
                }

    @classmethod
    def from_envfile(cls):
        api_key = os.getenv("NYTIMES_API_KEY")
        if not api_key:
            raise ValueError("NYTIMES_API_KEY not found in environment variables.")
        return cls(api_key)
