from typing import Optional, Iterable, Tuple
from contextlib import AbstractContextManager
from bs4 import BeautifulSoup
import httpx
from httpx_retries import RetryTransport, Retry
from httpx_ratelimiter import LimiterTransport
from hishel import CacheTransport, Controller, FileStorage

USER_AGENT = "scifi-fantasy/0.1 (dev@zmraines.com)"


def create_caching_client(
    api_root: str, base_transport: httpx.BaseTransport = None
) -> httpx.Client:
    """
    Create an HTTP client with caching and retry capabilities.
    """
    # Add retry capability to transport layer
    retry_transport = RetryTransport(
        transport=base_transport or httpx.HTTPTransport(http2=True),
        retry=Retry(
            total=3,
            backoff_factor=0.5,
            status_forcelist=[429],  # Retry on rate limit errors
        ),
    )

    # Add file-based caching on top of retry transport
    cache_transport = CacheTransport(
        transport=retry_transport,
        storage=FileStorage(ttl=1_000_000),
        controller=Controller(
            allow_stale=True,
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


class WikipediaDescriptionProvider(AbstractContextManager):
    """
    Attempt to extract a short summary of the plot of a work from Wikipedia
    by calling the MediaWiki API and stripping HTML tags.
    """

    WIKI_API_URL = "https://en.wikipedia.org/w/"

    def __init__(self):
        # HTTP client preconfigured with caching and retries
        self.client = create_caching_client(self.WIKI_API_URL)

    def __exit__(self, *args):
        self.client.__exit__(*args)

    def _article_name(self, wikipedia_url: str) -> str:
        """
        Extract the article title from the Wikipedia URL.
        """
        return wikipedia_url.removeprefix("https://en.wikipedia.org/wiki/")

    def get_sections(self, title: str) -> Iterable[Tuple[int, str]]:
        """
        Fetch the sections of the article.
        """
        resp = self.client.get(
            "api.php",
            params={
                "action": "parse",
                "page": title,
                "prop": "sections",
                "format": "json",
            },
        )
        resp.raise_for_status()
        data = resp.json().get("parse", {}).get("sections", [])
        for s in data:
            yield s["index"], s["line"]

    def get_section_html(self, article_title: str, section_idx: int) -> Optional[str]:
        """
        Fetch the HTML for a section of the article.
        """
        resp = self.client.get(
            "api.php",
            params={
                "action": "parse",
                "page": article_title,
                "prop": "text",
                "section": section_idx,
                "format": "json",
            },
        )
        resp.raise_for_status()
        return resp.json().get("parse", {}).get("text", {}).get("*", None)

    def get_description(self, wikipedia_url: str) -> Optional[str]:
        """
        Fetch the first 'Plot' section from the page and return its plain-text.
        """
        article_title = self._article_name(wikipedia_url)

        sections = self.get_sections(article_title)

        if not (
            plot_idx := next((s[0] for s in sections if s[1].startswith("Plot")), None)
        ):
            return None

        if not (html := self.get_section_html(article_title, plot_idx)):
            return None
        return BeautifulSoup(html, "html.parser").get_text(separator="\n").strip()


class GoogleBooksDescriptionProvider(AbstractContextManager):
    BASE_URL = "https://www.googleapis.com/books/v1"

    def __init__(self, api_key: str):
        self.api_key = api_key
        self.client = create_caching_client(
            self.BASE_URL,
            base_transport=LimiterTransport(
                per_minute=100,
                per_day=1000,
                max_delay=6000,
            ),
        )

    def __exit__(self, *args):
        self.client.__exit__(*args)

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
