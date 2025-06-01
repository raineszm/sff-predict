import httpx
import urllib.parse
from diskcache import Cache
from httpx_retries import Retry, RetryTransport
from typing import Callable, Optional

import hashlib
import os

USER_AGENT = "scifi-fantasy/0.1 (dev@zmraines.com)"


class CachedApi:
    def __init__(
        self,
        cache_name: str,
        api_root: str,
        base_transport: httpx.BaseTransport = None,
        key_fn: Callable[[httpx.Request], str] = None,
    ):
        self.api_root = api_root
        self.cache = Cache(directory=os.path.join(".cache", cache_name))
        self.client = self.create_caching_client(base_transport)
        self.key_fn = key_fn or query_to_cache_key

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
                total=5,
                backoff_factor=0.5,
                status_forcelist=[429],  # Retry on rate limit errors
            ),
        )

        return httpx.Client(
            http2=True,  # Enable HTTP/2 for better performance
            headers={"User-Agent": USER_AGENT},
            transport=retry_transport,
        )

    def _check_cache(self, request: httpx.Request) -> Optional[dict]:
        return self.cache.get(self.key_fn(request), default=None)

    def _send_request(self, request: httpx.Request) -> dict:
        if cached_response := self._check_cache(request):
            return cached_response
        resp = self.client.send(request)
        resp.raise_for_status()
        data = resp.json()
        self.cache.set(self.key_fn(request), data, expire=None)
        return data

    def call_endpoint(self, path: str = "", params: Optional[dict] = None) -> dict:
        if params is None:
            params = {}
        return self._send_request(
            httpx.Request(
                "GET",
                urllib.parse.urljoin(self.api_root, path),
                params=params,
            )
        )


def query_to_cache_key(request: httpx.Request) -> str:
    return hashlib.sha256(request.url.query).hexdigest()


def path_to_cache_key(request: httpx.Request) -> str:
    return hashlib.sha256(request.url.path.encode()).hexdigest()
