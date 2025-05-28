"""
Metadata providers for looking up work identifiers from various sources.

This module provides classes for resolving book/work identifiers across different
databases and services, including local lookups and external API calls to
OpenLibrary.
"""

from typing import Optional
import httpx
from httpx_retries import RetryTransport, Retry
from hishel import CacheTransport, Controller, FileStorage


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
        headers={
            "User-Agent": "scifi-fantasy/0.1 (dev@zmraines.com)",
        },
        transport=cache_transport,
        base_url=api_root,
    )


class OpenLibraryClient:
    """
    Provides lookup functionality using the OpenLibrary API.
    """

    def __init__(self):
        self.client = create_caching_client("https://openlibrary.org/")

    def __search(self, query: str) -> Optional[dict]:
        """Search OpenLibrary for a single result matching the query."""
        response = self.client.get(
            "search.json",
            params={
                "q": query,
                "fields": "key,description",  # only get the fields we need
                "limit": 1,
            },
        )
        response.raise_for_status()
        json = response.json()
        if json["numFound"] != 0:
            return json["docs"][0]

    def lookup_by_id(self, id: str) -> Optional[dict]:
        """Look up a work by its OpenLibrary identifier."""
        return self.__search(f"key:/works/{id}")

    def lookup_by_isbn(self, isbn: str) -> Optional[dict]:
        """Look up a work by its ISBN."""
        return self.__search(f"isbn:{isbn}")

    def lookup_by_title_authors(self, title: str, authors: list[str]) -> Optional[dict]:
        """Look up a work by its title and authors."""
        return self.__search(f"intitle:{title} inauthors:{','.join(authors)}")
