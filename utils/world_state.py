from .api import CachedApi, path_to_cache_key
import os
from typing import Iterator, TypedDict
from httpx import HTTPTransport, Response, Request
from time import sleep
from loguru import logger
from httpx import HTTPStatusError


class NYTArticle(TypedDict):
    id: str
    headline: str
    abstract: str
    year: int
    month: int


class NYTLimiterTransport(HTTPTransport):
    """
    A transport that limits the number of requests to the NYT API.

    The NYT API has a rate limit of 5 requests per second, and they're very strict about it.
    This transport will sleep for 12 seconds between requests in order to avoid being rate limited.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def handle_request(self, request: Request, **kwargs) -> Response:
        logger.debug("Okay, calling the API.")
        response = super().handle_request(request, **kwargs)
        logger.debug("Sleeping for 12 seconds to avoid being rate limited")
        sleep(12)  # Sleep for 12 seconds to avoid being rate limited
        logger.debug("Done sleeping.")
        return response


class NYTimesArchiveAPI:
    API_BASE = "https://api.nytimes.com/svc/archive/v1/"

    def __init__(self, api_key: str):
        self.api_key = api_key
        self._api = CachedApi(
            "nytimes",
            self.API_BASE,
            key_fn=path_to_cache_key,
            transport=NYTLimiterTransport(http2=True),
            diskcache_settings={"eviction_policy": "none"},
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
