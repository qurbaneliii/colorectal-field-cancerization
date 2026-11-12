from __future__ import annotations

import gzip
import hashlib
import io
import tarfile
from email.message import Message

import pytest

from src.data.geo_archive import compare_archives, inspect_cel_archive
from src.data.resumable_download import resumable_download


def test_archive_distinguishes_repacking_recompression_and_changed_payload(tmp_path):
    name = "GSM1001_sample.CEL.gz"
    payload = b"HG-U219 header and intensity bytes"
    spec = {
        "accession": "test", "url": "", "bytes": 0, "historical_sha256": "",
        "chip": "HG-U219", "sample_by_file": {name: "GSM1001"}, "historical_order": [name],
    }

    def build(label, cel, tar_mtime, gzip_mtime):
        path = tmp_path / f"{label}.tar"
        compressed = gzip.compress(cel, mtime=gzip_mtime)
        with tarfile.open(path, "w") as out:
            member = tarfile.TarInfo(name)
            member.size = len(compressed)
            member.mtime = tar_mtime
            out.addfile(member, io.BytesIO(compressed))
        return inspect_cel_archive(path, spec)[1]

    original = build("original", payload, 100, 100)
    repack = compare_archives(build("repack", payload, 120, 100), original)
    assert repack["tar_header_differences"] == [name]
    assert repack["gzip_differences"] == repack["cel_differences"] == []
    recompression = compare_archives(build("recompress", payload, 120, 120), original)
    assert recompression["gzip_differences"] == [name]
    assert recompression["cel_differences"] == []
    changed = compare_archives(build("changed", payload + b"!", 120, 120), original)
    assert changed["cel_differences"] == [name]


def test_resumed_download_validates_range_then_reuses_verified_cache(tmp_path, monkeypatch):
    data = b"known input bytes"
    target = tmp_path / "download.tar"
    partial = tmp_path / "download.tar.part"
    partial.write_bytes(data[:5])
    digest = hashlib.sha256(data).hexdigest()
    calls = []

    class RangeResponse(io.BytesIO):
        status = 206
        headers = Message()
        headers["Content-Range"] = f"bytes 5-{len(data) - 1}/{len(data)}"

    def open_request(request, **_):
        calls.append(request.get_header("Range"))
        return RangeResponse(data[5:])

    monkeypatch.setattr("urllib.request.urlopen", open_request)
    result = resumable_download("https://example.invalid/data", target, len(data), digest)
    assert result["initial_offset"] == 5 and target.read_bytes() == data
    assert not partial.exists()
    assert resumable_download("https://example.invalid/data", target, len(data), digest)["mode"] == (
        "verified_cache"
    )
    assert calls == ["bytes=5-"]


def test_download_preserves_full_partial_on_checksum_mismatch(tmp_path):
    target = tmp_path / "download.tar"
    partial = tmp_path / "download.tar.part"
    partial.write_bytes(b"wrong")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        resumable_download("https://example.invalid/data", target, 5, "0" * 64)
    assert partial.read_bytes() == b"wrong"
    assert not target.exists()


def test_download_refuses_a_resume_at_the_wrong_offset(tmp_path, monkeypatch):
    target = tmp_path / "download.tar"
    partial = tmp_path / "download.tar.part"
    partial.write_bytes(b"abc")

    class WrongRange(io.BytesIO):
        status = 206
        headers = Message()
        headers["Content-Range"] = "bytes 2-5/6"

    monkeypatch.setattr("urllib.request.urlopen", lambda *_, **__: WrongRange(b"cdef"))
    with pytest.raises(ValueError, match="Unexpected resume response"):
        resumable_download("https://example.invalid/data", target, 6, "0" * 64)
    assert partial.read_bytes() == b"abc"
    assert not target.exists()
