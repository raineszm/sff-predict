from dataclasses import dataclass, field
from typing import Optional
import pandas as pd
import httpx
from httpx_retries import RetryTransport, Retry
from hishel import CacheTransport, Controller, FileStorage


def create_caching_client(api_root: str) -> httpx.Client:
    # Be a good citizen
    retry = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429],
    )

    retry_transport = RetryTransport(
        retry=retry,
    )

    cache_transport = CacheTransport(
        transport=retry_transport,
        storage=FileStorage(),
        controller=Controller(
            always_revalidate=False,
            force_cache=True,
        ),
    )

    return httpx.Client(
        http2=True,
        headers={
            "User-Agent": "scifi-fantasy/0.1 (dev@zmraines.com)",
        },
        transport=cache_transport,
        base_url=api_root,
    )


@dataclass
class AwardIdentifiers:
    work_id: Optional[int] = None
    isdfb: Optional[str] = None
    isbn: list[str] = field(default_factory=list)
    openlibrary: list[str] = field(default_factory=list)

    def __post_init__(self):
        if self.work_id is not None:
            self.work_id = int(self.work_id)

    @classmethod
    def from_award(cls, award: pd.Series) -> "AwardIdentifiers":
        if award.work_id is pd.NA or award.work_id == "":
            work_id = None
        return cls(
            work_id=work_id,
            isdfb=award.isdfb_ids,
            isbn=award.isbns.split(";"),
            openlibrary=award.ol_ids.split(";"),
        )


class LocalIdentifierMetadataProvider:
    def __init__(self, identifiers: pd.DataFrame):
        self._identifiers = identifiers
        self._valid_work_ids = set(identifiers.work_id)

    def is_valid(self, work_id: int) -> bool:
        return work_id in self._valid_work_ids

    def lookup_by_kind(self, kind: str, value: str) -> Optional[int]:
        by_kind = self._identifiers.loc[
            (self._identifiers.kind == kind) & (self._identifiers.value == value)
        ].work_id
        if not by_kind.empty:
            return by_kind.iloc[0]
        return None

    def as_work_id(self, work_or_book_id: int) -> Optional[int]:
        if work_or_book_id in self._valid_work_ids:
            return work_or_book_id
        if (work_id := self.lookup_by_kind("book_id", work_or_book_id)) is not None:
            return work_id
        return None

    def lookup(self, award_ids: AwardIdentifiers) -> Optional[int]:
        if award_ids.work_id is not None:
            return self.as_work_id(award_ids.work_id)
        if award_ids.isbn is not None:
            for isbn in award_ids.isbn:
                isbn = isbn.replace("-", "")
                kind = "isbn" if len(isbn) == 10 else "isbn13"
                if (work_id := self.lookup_by_kind(kind, isbn)) is not None:
                    return work_id
        return None


class OpenLibraryMetadataProvider:
    def __init__(self, local_provider: LocalIdentifierMetadataProvider):
        self.client = create_caching_client("https://openlibrary.org/")
        self._local_provider = local_provider

    def __search(self, query: str) -> Optional[dict]:
        response = self.client.get(
            "search.json",
            params={
                "q": query,
                "fields": "id_goodreads,isbn,id_amazon",
                "limit": 1,
            },
        )
        response.raise_for_status()
        json = response.json()
        if json["numFound"] != 0:
            return json["docs"][0]

    def _lookup_by_ol_id(self, ol_id: str) -> Optional[int]:
        doc = self.__search(ol_id)
        if doc is None:
            return None

        for id in doc.get("id_goodreads", []):
            work_id = self._local_provider.as_work_id(id)
            if work_id is not None:
                return work_id

        for isbn in doc.get("isbn", []):
            if (
                work_id := self._local_provider.lookup_by_kind("isbn", isbn)
            ) is not None:
                return work_id

        for id in doc.get("id_amazon", []):
            work_id = self._local_provider.lookup_by_kind("amazon", id)
            if work_id is not None:
                return work_id

        return None

    def lookup(self, award_ids: AwardIdentifiers) -> Optional[int]:
        if award_ids.openlibrary is not None:
            for ol_id in award_ids.openlibrary:
                if (work_id := self._lookup_by_ol_id(ol_id)) is not None:
                    return work_id
        return None
