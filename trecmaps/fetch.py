"""Downloads: the only code in trecmaps that touches the network (used by scripts/0*_get_*.py)."""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

USER_AGENT = "trec-maps/0.1 (https://github.com/savvas-paragkamian; TREC transect maps)"


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            size = resp.headers.get("Content-Length")
            return int(size) if size else None
    except Exception:
        return None


def download(url: str, dest: Path, *, dry_run: bool = False, log=print) -> bool:
    """Stream `url` to `dest` via `dest.part`; skip when `dest` exists. Returns True when the file is present."""
    dest = Path(dest)
    if dest.exists():
        log(f"have {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB), skipped")
        return True
    if dry_run:
        size = head_size(url)
        log(f"would download {url}" + (f" ({size / 2**20:.1f} MiB)" if size else ""))
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=300) as resp, open(part, "wb") as f:
            shutil.copyfileobj(resp, f, length=1 << 20)
    except urllib.error.HTTPError as exc:
        part.unlink(missing_ok=True)
        if exc.code == 404:
            log(f"404 for {url} (an all-ocean DEM tile is not published), skipped")
            return False
        raise
    part.rename(dest)
    log(f"saved {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB)")
    return True
