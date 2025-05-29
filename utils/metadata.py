from typing import Optional

import wikipediaapi
import httpx
from httpx_retries import RetryTransport, Retry
from hishel import CacheTransport, Controller, FileStorage

USER_AGENT = "scifi-fantasy/0.1 (dev@zmraines.com)"


def create_caching_client(api_root: str) -> httpx.Client:
    """
    Create an HTTP client with caching and retry capabilities.
    """
    retry = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429],  # Retry on rate limit errors
    )

    # Add retry capability to transport layer
    retry_transport = RetryTransport(
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
