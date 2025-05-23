"""
Metadata providers for looking up work identifiers from various sources.

This module provides classes for resolving book/work identifiers across different
databases and services, including local lookups and external API calls to
OpenLibrary.
"""

from dataclasses import dataclass, field
from typing import Optional
import pandas as pd
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


@dataclass
class AwardIdentifiers:
    """
    Container for various identifiers that can be used to look up a work.

    Attributes:
        work_id: Goodreads work identifier
        isdfb: Internet Speculative Fiction Database identifier
        isbn: List of ISBN numbers (both ISBN-10 and ISBN-13)
        openlibrary: List of OpenLibrary identifiers
    """

    work_id: Optional[int] = None
    isdfb: Optional[str] = None
    isbn: list[str] = field(default_factory=list)
    openlibrary: list[str] = field(default_factory=list)

    def __post_init__(self):
        if self.work_id is not None:
            self.work_id = int(self.work_id)

        self.isbn = [isbn.replace("-", "").strip() for isbn in self.isbn]

    @classmethod
    def from_award(cls, award: pd.Series) -> "AwardIdentifiers":
        # Handle missing or empty work_id values
        work_id = (
            None if (award.work_id is pd.NA or award.work_id == "") else award.work_id
        )
        return cls(
            work_id=work_id,
            isdfb=award.isdfb_ids,
            isbn=award.isbns.split(";"),
            openlibrary=award.ol_ids.split(";"),
        )


class LocalIdentifierMetadataProvider:
    """
    Provides lookup functionality for work identifiers using the identifiers
    data file.
    """

    def __init__(self, identifiers: pd.DataFrame):
        self._identifiers = identifiers
        self._valid_work_ids = set(identifiers.work_id)

    def is_valid(self, work_id: int) -> bool:
        """Check if a work ID exists in the local database."""
        return work_id in self._valid_work_ids

    def lookup_by_kind(self, kind: str, value: str) -> Optional[int]:
        """Look up a work ID by identifier type and value."""
        by_kind = self._identifiers.loc[
            (self._identifiers.kind == kind) & (self._identifiers.value == value)
        ].work_id
        if not by_kind.empty:
            return by_kind.iloc[0]
        return None

    def as_work_id(self, work_or_book_id: int) -> Optional[int]:
        """Convert a work or book ID to a canonical work ID."""
        # Check if it's already a valid work ID
        if work_or_book_id in self._valid_work_ids:
            return work_or_book_id
        # Try to resolve as a book ID to its parent work
        if (work_id := self.lookup_by_kind("book_id", work_or_book_id)) is not None:
            return work_id
        return None

    def lookup(self, award_ids: AwardIdentifiers) -> Optional[int]:
        """Look up a work ID using any of the identifiers in AwardIdentifiers."""
        # Try direct work ID first (most reliable)
        if award_ids.work_id is not None:
            return self.as_work_id(award_ids.work_id)

        # Try ISBN lookup if available
        if award_ids.isbn:
            for isbn in award_ids.isbn:
                # Determine ISBN type based on length
                kind = "isbn13" if len(isbn) == 13 else "isbn"
                if (work_id := self.lookup_by_kind(kind, isbn)) is not None:
                    return work_id
        return None


class OpenLibraryMetadataProvider:
    """
    Provides lookup functionality using the OpenLibrary API.
    """

    def __init__(self, local_provider: LocalIdentifierMetadataProvider):
        self.client = create_caching_client("https://openlibrary.org/")
        self._local_provider = local_provider

    def __search(self, query: str) -> Optional[dict]:
        """Search OpenLibrary for a single result matching the query."""
        response = self.client.get(
            "search.json",
            params={
                "q": query,
                "fields": "id_goodreads,isbn,id_amazon",  # only get the fields we need
                "limit": 1,
            },
        )
        response.raise_for_status()
        json = response.json()
        if json["numFound"] != 0:
            return json["docs"][0]

    def _lookup_by_ol_id(self, ol_id: str) -> Optional[int]:
        """Look up a work ID using an OpenLibrary identifier."""
        doc = self.__search(ol_id)
        if doc is None:
            return None
        # Try each type of ID in sequence
        id_lookups = [
            ("id_goodreads", self._local_provider.as_work_id),
            ("isbn", lambda x: self._local_provider.lookup_by_kind("isbn", x)),
            ("id_amazon", lambda x: self._local_provider.lookup_by_kind("amazon", x)),
        ]

        for doc_key, lookup_fn in id_lookups:
            if work_id := try_ids(doc.get(doc_key, []), lookup_fn):
                return work_id

        return None

    def lookup(self, award_ids: AwardIdentifiers) -> Optional[int]:
        """Look up a work ID using OpenLibrary identifiers."""
        if award_ids.openlibrary:
            for ol_id in award_ids.openlibrary:
                if (work_id := self._lookup_by_ol_id(ol_id)) is not None:
                    return work_id
        return None


def try_ids(ids: list, lookup_method) -> Optional[int]:
    """Helper function to try a lookup method on a list of identifiers."""
    for id_value in ids:
        if work_id := lookup_method(id_value):
            return work_id
    return None
