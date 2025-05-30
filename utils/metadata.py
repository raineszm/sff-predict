from typing import Optional, Iterable, Tuple
from contextlib import AbstractContextManager
from bs4 import BeautifulSoup
import httpx
from httpx_retries import RetryTransport, Retry
from httpx_ratelimiter import LimiterTransport
import os
from diskcache import Cache
import hashlib

USER_AGENT = "scifi-fantasy/0.1 (dev@zmraines.com)"


class CachedApi:
    def __init__(
        self, cache_name: str, api_root: str, base_transport: httpx.BaseTransport = None
    ):
        self.api_root = api_root
        self.cache = Cache(directory=os.path.join(".cache", cache_name))
        self.client = self.create_caching_client(base_transport)

    def __exit__(self, *args):
        self.client.__exit__(*args)

    @staticmethod
    def create_caching_client(
        base_transport: httpx.BaseTransport = None,
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

        return httpx.Client(
            http2=True,  # Enable HTTP/2 for better performance
            headers={"User-Agent": USER_AGENT},
            transport=retry_transport,
        )

    def _cache_key(self, request: httpx.Request) -> str:
        return hashlib.sha256(request.url.query).hexdigest()

    def _check_cache(self, request: httpx.Request) -> Optional[dict]:
        return self.cache.get(self._cache_key(request), default=None)

    def call_endpoint(self, params: dict) -> dict:
        request = httpx.Request(
            "GET",
            self.api_root,
            params=params,
        )
        if cached_response := self._check_cache(request):
            return cached_response
        resp = self.client.send(request)
        resp.raise_for_status()
        data = resp.json()
        self.cache.set(self._cache_key(request), data, expire=None)
        return data


class WikipediaDescriptionProvider(AbstractContextManager):
    """
    Attempt to extract a short summary of the plot of a work from Wikipedia
    by calling the MediaWiki API and stripping HTML tags.
    """

    WIKI_API_URL = "https://en.wikipedia.org/w/api.php"

    def __init__(self):
        # HTTP client preconfigured with caching and retries
        self.client = CachedApi("wikipedia", self.WIKI_API_URL)

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
        resp = self.client.call_endpoint(
            params={
                "action": "parse",
                "page": title,
                "prop": "sections",
                "format": "json",
            },
        )
        data = resp.get("parse", {}).get("sections", [])
        for s in data:
            yield s["index"], s["line"]

    def get_section_html(self, article_title: str, section_idx: int) -> Optional[str]:
        """
        Fetch the HTML for a section of the article.
        """
        resp = self.client.call_endpoint(
            params={
                "action": "parse",
                "page": article_title,
                "prop": "text",
                "section": section_idx,
                "format": "json",
            },
        )
        return resp.get("parse", {}).get("text", {}).get("*", None)

    @staticmethod
    def is_synopsis(section: str) -> bool:
        """
        Heuristic to determine if a section is a story synopsis.
        """
        return any(
            map(
                section.casefold().startswith,
                [
                    "plot",
                    "synopsis",
                    "premise",
                    "summary",
                    "overview",
                    "setting",
                    "back story",
                    "backstory",
                ],
            )
        )

    def get_summary(self, article_title: str) -> Optional[str]:
        """
        Get summary of the article.
        """
        resp = self.client.call_endpoint(
            params={
                "action": "query",
                "titles": article_title,
                "prop": "extracts",
                "exsectionformat": "plain",
                "explaintext": True,
                "exintro": True,
                "format": "json",
            },
        )
        return list(resp.get("query", {}).get("pages", {}).values())[0].get("extract")

    def get_description(self, wikipedia_url: str) -> Optional[str]:
        """
        Fetch the first 'Plot' section from the page and return its plain-text.
        """
        article_title = self._article_name(wikipedia_url)

        sections = self.get_sections(article_title)

        if not (
            plot_idx := next((s[0] for s in sections if self.is_synopsis(s[1])), None)
        ):
            if summary := self.get_summary(article_title):
                return summary
            return None

        if not (html := self.get_section_html(article_title, plot_idx)):
            return None
        return BeautifulSoup(html, "html.parser").get_text(separator="\n").strip()


class GoogleBooksDescriptionProvider(AbstractContextManager):
    BASE_URL = "https://www.googleapis.com/books/v1/volumes"

    def __init__(self, api_key: str):
        self.api_key = api_key
        self.client = CachedApi(
            "googlebooks",
            self.BASE_URL,
            base_transport=LimiterTransport(
                per_minute=100,
                per_day=1000,
                max_delay=60_000,
            ),
        )

    def __exit__(self, *args):
        self.client.__exit__(*args)

    def _query(self, query: str) -> dict:
        return self.client.call_endpoint(
            params={
                "q": query,
                "key": self.api_key,
                "projection": "LITE",
                "printType": "BOOKS",
                "fields": "totalItems,items(volumeInfo/description)",
            },
        )

    def get_description(self, isbn: str) -> Optional[str]:
        data = self._query(f"isbn:{isbn}")
        if data.get("totalItems", 0) == 0:
            return None
        return data.get("items", [{}])[0].get("volumeInfo", {}).get("description")
