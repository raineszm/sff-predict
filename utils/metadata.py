from typing import Optional

import wikipediaapi
import httpx
from httpx_retries import RetryTransport, Retry
from httpx_ratelimiter import LimiterTransport
from hishel import CacheTransport, Controller, FileStorage

USER_AGENT = "scifi-fantasy/0.1 (dev@zmraines.com)"


def create_caching_client(api_root: str) -> httpx.Client:
    """
    Create an HTTP client with caching and retry capabilities.
    """
    # Enforce google api rate limits
    limiter = LimiterTransport(
        per_minute=100,
        per_day=1000,
    )

    retry = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429],  # Retry on rate limit errors
    )

    # Add retry capability to transport layer
    retry_transport = RetryTransport(
        transport=limiter,
        retry=retry,
    )

    # Add file-based caching on top of retry transport
    cache_transport = CacheTransport(
        transport=retry_transport,
        storage=FileStorage(),
        controller=Controller(
            always_revalidate=False,
            force_cache=True,  # Prefer cache to reduce API load
        ),
    )

    return httpx.Client(
        http2=True,  # Enable HTTP/2 for better performance
        headers={"User-Agent": USER_AGENT},
        transport=cache_transport,
        base_url=api_root,
    )


class WikipediaDescriptionProvider:
    """
    Attempt to extract a short summary of the plot of a work from Wikipedia.
    """

    def __init__(self):
        self.client = wikipediaapi.Wikipedia(user_agent=USER_AGENT)

    def _article_name(self, wikipedia_url: str) -> str:
        """
        Extract the name of the article from the Wikipedia URL.
        """
        return wikipedia_url.removeprefix("https://en.wikipedia.org/wiki/")

    def get_description(self, wikipedia_url: str) -> Optional[str]:
        page = self.client.page(self._article_name(wikipedia_url))
        if not page.exists():
            return None

        return next(
            map(
                lambda s: s.text,
                filter(lambda s: s.title.startswith("Plot"), page.sections),
            ),
            None,
        )


class GoogleBooksDescriptionProvider:
    BASE_URL = "https://www.googleapis.com/books/v1"

    def __init__(self, api_key: str):
        self.api_key = api_key
        self.client = create_caching_client(self.BASE_URL)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.client.__exit__(exc_type, exc_value, traceback)

    def _query(self, query: str) -> dict:
        resp = self.client.get(
            "/volumes",
            params={
                "q": query,
                "key": self.api_key,
                "projection": "LITE",
                "printType": "BOOKS",
                "fields": "totalItems,items(volumeInfo/description)",
            },
        )
        resp.raise_for_status()
        return resp.json()

    def get_description(self, isbn: str) -> Optional[str]:
        data = self._query(f"isbn:{isbn}")
        if data.get("totalItems", 0) == 0:
            return None
        return data.get("items", [{}])[0].get("volumeInfo", {}).get("description")
