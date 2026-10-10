#!/usr/bin/env python3
"""Download Open Targets Platform datasets (Parquet) with resume.

  uv run python scripts/11_fetch_opentargets.py --list      # sizes only, downloads nothing
  uv run python scripts/11_fetch_opentargets.py             # download the default datasets
  uv run python scripts/11_fetch_opentargets.py --extras    # also the optional extras
  uv run python scripts/11_fetch_opentargets.py --datasets disease target

Ctrl-C is safe: re-run the same command and it continues (finished files are skipped,
partial *.part files resume from where they stopped). Files go to data/raw/opentargets/<release>/<dataset>/.
"""
import argparse
import json
import os
import re
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx

BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
DEFAULT = ["disease", "target", "drug_molecule", "drug_mechanism_of_action", "clinical_indication",
           "clinical_report", "clinical_target", "association_overall_indirect"]
EXTRAS = ["target_tractability", "chemical_probe", "drug_warning", "association_overall_direct"]
UA = "cell-interaction-map/0.1 (personal research project; github.com/alirezahekmati/cell-interaction-map)"
lock = threading.Lock()


def gb(n):
    return f"{n / 1e9:.2f} GB"


def get_text(client, url, tries=6):
    for i in range(tries):
        try:
            r = client.get(url)
            if r.status_code == 404:
                return ""
            r.raise_for_status()
            return r.text
        except httpx.HTTPError:
            time.sleep(min(60, 3 * 2 ** i))
    raise RuntimeError(f"cannot list {url}")


def list_files(client, url, depth=0):
    """Relative paths of the .parquet files under a directory listing."""
    out = []
    for h in re.findall(r'href="([^"?#]+)"', get_text(client, url)):
        if h.startswith(("/", "http", ".", "_")):
            continue
        if h.endswith("/"):
            if depth < 3:
                out += [h + f for f in list_files(client, url + h, depth + 1)]
        elif h.endswith(".parquet"):
            out.append(h)
    return out


def remote_size(client, url, cache):
    if url in cache:
        return cache[url]
    for i in range(6):
        try:
            r = client.head(url)
            n = r.headers.get("content-length")
            size = int(n) if (r.status_code == 200 and n) else None
            if size is None:
                r = client.get(url, headers={"Range": "bytes=0-0"})
                m = re.search(r"/(\d+)$", r.headers.get("content-range", ""))
                size = int(m.group(1)) if m else None
            if size is not None:
                cache[url] = size
                return size
        except httpx.HTTPError:
            pass
        time.sleep(min(60, 3 * 2 ** i))
    raise RuntimeError(f"cannot get size of {url}")


def complete(dest, size):
    return dest.exists() and dest.stat().st_size == size


def fetch(client, url, dest, size, pause):
    time.sleep(pause)
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    t0, start_have = time.time(), (part.stat().st_size if part.exists() else 0)
    for attempt in range(1, 9):
        try:
            have = part.stat().st_size if part.exists() else 0
            if have > size:
                part.unlink()
                have = 0
            if have < size:
                headers = {"Range": f"bytes={have}-"} if have else {}
                with client.stream("GET", url, headers=headers) as r:
                    if r.status_code not in (200, 206):
                        raise httpx.HTTPError(f"HTTP {r.status_code}")
                    mode = "ab" if (have and r.status_code == 206) else "wb"
                    with open(part, mode) as f:
                        for chunk in r.iter_bytes(1 << 20):
                            f.write(chunk)
            if part.stat().st_size == size:
                part.rename(dest)
                mb = (size - start_have) / 1e6
                with lock:
                    print(f"    ok  {dest.name[:48]:48s} {size / 1e6:8.1f} MB  {mb / max(time.time() - t0, 1e-6):6.2f} MB/s", flush=True)
                return True
        except (httpx.HTTPError, OSError) as e:
            wait = min(300, 5 * 2 ** (attempt - 1))
            with lock:
                print(f"    retry {attempt}/8 {dest.name[:40]}: {e!r}; waiting {wait}s", flush=True)
            time.sleep(wait)
    with lock:
        print(f"    FAILED {dest.name}", flush=True)
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default="26.09")
    ap.add_argument("--datasets", nargs="+")
    ap.add_argument("--extras", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--pause", type=float, default=0.3)
    a = ap.parse_args()

    root = Path(os.environ.get("BIO_ROOT", ".")) / "data/raw/opentargets" / a.release
    root.mkdir(parents=True, exist_ok=True)
    cache_f, log_f = root / "sizes.json", root / "done.jsonl"
    cache = json.loads(cache_f.read_text()) if cache_f.exists() else {}
    names = a.datasets or DEFAULT + (EXTRAS if a.extras else [])
    client = httpx.Client(follow_redirects=True, timeout=httpx.Timeout(30.0, read=180.0), headers={"User-Agent": UA})

    print(f"Open Targets release {a.release}", flush=True)
    plan = {}
    for ds in names:
        base = f"{BASE}/{a.release}/output/{ds}/"
        files = list_files(client, base)
        if not files:
            print(f"  {ds:34s} NOT FOUND or empty, skipped", flush=True)
            continue
        with ThreadPoolExecutor(a.workers) as ex:
            sizes = list(ex.map(lambda f: remote_size(client, base + f, cache), files))
        cache_f.write_text(json.dumps(cache))
        plan[ds] = [(base + f, root / ds / f, s) for f, s in zip(files, sizes)]
        have = sum(s for _u, d, s in plan[ds] if complete(d, s))
        print(f"  {ds:34s} {len(files):4d} files {gb(sum(sizes)):>10s}   (already have {gb(have)})", flush=True)

    need = sum(s for items in plan.values() for _u, d, s in items if not complete(d, s))
    free = shutil.disk_usage(root).free
    print(f"\nstill to download: {gb(need)} | free disk: {gb(free)}", flush=True)
    if a.list:
        return 0
    if need > free * 0.9:
        print("Not enough free disk space; free some space or use --datasets to choose fewer.")
        return 2

    failed, t0 = [], time.time()
    for ds, items in plan.items():
        todo = [(u, d, s) for u, d, s in items if not complete(d, s)]
        print(f"\n== {ds}: {len(todo)} of {len(items)} files to download ==", flush=True)
        with ThreadPoolExecutor(a.workers) as ex:
            futs = {ex.submit(fetch, client, u, d, s, a.pause): (d, s) for u, d, s in todo}
            for fu in as_completed(futs):
                d, s = futs[fu]
                if fu.result():
                    with lock, log_f.open("a") as f:
                        f.write(json.dumps({"path": str(d.relative_to(root)), "size": s, "t": int(time.time())}) + "\n")
                else:
                    failed.append(str(d))
        print(f"== {ds} finished | elapsed {(time.time() - t0) / 60:.1f} min | failures so far: {len(failed)} ==", flush=True)

    print("\nDONE." if not failed else f"\nDONE with {len(failed)} failed files; re-run the same command to retry them:")
    for f in failed:
        print("  ", f)
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\ninterrupted; re-run the same command to resume", flush=True)
        os._exit(130)
