from __future__ import annotations

import http.client
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

from src.data.provenance import file_sha256


def resumable_download(
    url: str,
    destination: Path,
    size: int,
    sha256: str,
    *,
    retries: int = 6,
    timeout: float = 60,
) -> dict:
    """Resume a pinned file, promoting the partial only after a full SHA-256 check."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.stat().st_size != size or file_sha256(destination) != sha256:
            raise ValueError(f"Verified cache path contains unexpected bytes: {destination}")
        return {"mode": "verified_cache", "sha256": sha256, "size": size}

    partial = destination.with_name(destination.name + ".part")
    checkpoint = partial.with_suffix(".json")
    request_pin = {"url": url, "size": size, "sha256": sha256}
    state = json.loads(checkpoint.read_text(encoding="utf-8")) if checkpoint.exists() else request_pin
    if any(state.get(key) != value for key, value in request_pin.items()):
        raise ValueError(f"Existing partial belongs to a different source or pin: {partial}")

    def offset() -> int:
        return partial.stat().st_size if partial.exists() else 0

    initial_offset = offset()
    if initial_offset > size:
        raise ValueError(f"Partial is larger than the pinned size: {partial}")
    for attempt in range(retries + 1):
        start = offset()
        if start == size:
            break
        headers = {"Accept-Encoding": "identity"}
        if start:
            headers["Range"] = f"bytes={start}-"
            validator = state.get("etag") or state.get("last_modified")
            if validator:
                headers["If-Range"] = validator
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if response.headers.get("Content-Encoding", "identity").lower() != "identity":
                    raise ValueError("Encoded HTTP body cannot be resumed by byte offset")
                if response.status == 206:
                    content_range = response.headers.get("Content-Range", "")
                    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
                    if match is None or (int(match[1]), int(match[3])) != (start, size):
                        raise ValueError(f"Unexpected resume response: {content_range!r}")
                    end = int(match[2]) + 1
                elif start == 0 and response.status == 200:
                    end = size
                else:
                    raise ValueError("The server refused the pinned range. Preserved the partial file.")
                if not start < end <= size:
                    raise ValueError("Response range exceeds the pinned file size")
                etag = response.headers.get("ETag")
                validators = {
                    "etag": etag if etag and not etag.startswith("W/") else None,
                    "last_modified": response.headers.get("Last-Modified"),
                }
                for key, value in validators.items():
                    if state.get(key) and state[key] != value:
                        raise ValueError(f"Remote {key} changed. Preserved the partial file.")
                state = {**request_pin, **validators}
                checkpoint.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
                if offset() != start:
                    raise ValueError("Concurrent write to the partial file")
                with partial.open("ab") as output:
                    written = start
                    while written < end:
                        chunk = response.read(min(1024 * 1024, end - written))
                        if not chunk:
                            raise ConnectionError("HTTP response ended before its promised size")
                        output.write(chunk)
                        written += len(chunk)
                if end < size:
                    raise ConnectionError("The response covered only part of the requested suffix")
        except urllib.error.HTTPError as error:
            if error.code not in {408, 429, 500, 502, 503, 504} or attempt == retries:
                raise
        except (OSError, http.client.HTTPException):
            if attempt == retries:
                raise
        if offset() < size:
            time.sleep(min(2 ** attempt, 30))

    if not partial.exists() or offset() != size:
        raise ValueError(f"Incomplete download for {destination}")
    observed = file_sha256(partial)
    if observed != sha256:
        raise ValueError(f"SHA-256 mismatch ({observed}); preserved the partial file {partial}")
    if destination.exists():
        raise FileExistsError(destination)
    partial.rename(destination)
    checkpoint.unlink(missing_ok=True)
    return {
        "mode": "download" if initial_offset < size else "completed_partial",
        "sha256": observed,
        "size": size,
        "initial_offset": initial_offset,
        "etag": state.get("etag"),
        "last_modified": state.get("last_modified"),
    }
