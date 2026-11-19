"""Read-only, rate-limited client for the IGVF portal API.

Only GET is implemented: this module cannot write to the portal.
Rate is capped at HARD_MAX_RPS regardless of what the caller asks for.
"""
from __future__ import annotations

import sys
import threading
import time
from urllib.parse import urlencode

import requests

API = "https://api.data.igvf.org"
HARD_MAX_RPS = 10.0       # buildprompt: at most 10 requests/s, never exceeded
DEFAULT_RPS = 4.0
USER_AGENT = "igvfsim-research/0.1 (read-only; python-requests)"


class Portal:
    def __init__(self, rps: float = DEFAULT_RPS, timeout: float = 300, retries: int = 4):
        self.min_interval = 1.0 / min(rps, HARD_MAX_RPS)
        self.timeout, self.retries = timeout, retries
        self._last = 0.0
        self._lock = threading.Lock()
        self.n_requests = 0
        self.s = requests.Session()
        self.s.headers.update({"Accept": "application/json", "User-Agent": USER_AGENT})

    def _wait(self):
        with self._lock:
            dt = time.monotonic() - self._last
            if dt < self.min_interval:
                time.sleep(self.min_interval - dt)
            self._last = time.monotonic()

    def get(self, path: str, params: list[tuple[str, str]] | None = None) -> dict:
        url = API + path + (("?" + urlencode(params, safe="@.")) if params else "")
        for attempt in range(self.retries + 1):
            self._wait()
            self.n_requests += 1
            try:
                r = self.s.get(url, timeout=self.timeout)
                if r.status_code in (429, 500, 502, 503, 504):
                    raise requests.HTTPError(f"HTTP {r.status_code}", response=r)
                r.raise_for_status()
                return r.json()
            except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as e:
                status = getattr(getattr(e, "response", None), "status_code", None)
                if attempt == self.retries or (status and status < 500 and status != 429):
                    raise
                back = 2 ** attempt * 2
                print(f"  retry {attempt + 1} in {back}s: {e}", file=sys.stderr)
                time.sleep(back)
        raise RuntimeError("unreachable")

    def search(self, item_type: str, fields: list[str], *, frm: int = 0, limit: int | str = 500,
               status: str = "released", sort: str = "accession") -> dict:
        params = [("type", item_type), ("status", status), ("format", "json"),
                  ("limit", str(limit))]
        if limit != "all":
            params.append(("from", str(frm)))
        if sort:
            params.append(("sort", sort))
        params += [("field", f) for f in fields]
        return self.get("/search/", params)
