"""Download raw sources into data/raw/ and record a manifest (url, date, size, sha256)."""
import hashlib, json, zipfile
from datetime import datetime, timezone
from pathlib import Path
import httpx

RAW = Path("data/raw")
MANIFEST = RAW / "MANIFEST.json"

GEM = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model"
SOURCES = {
    "hpa/subcellular_location.tsv.zip":
        "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip",
    "uniprot/HUMAN_9606_idmapping.dat.gz":
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
    "humangem/metabolites.tsv": f"{GEM}/metabolites.tsv",
    "humangem/reactions.tsv": f"{GEM}/reactions.tsv",
    "humangem/genes.tsv": f"{GEM}/genes.tsv",
}

manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}

def record(rel, url, path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    manifest[rel] = {"url": url, "bytes": path.stat().st_size, "sha256": h.hexdigest(),
                     "fetched": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    MANIFEST.write_text(json.dumps(manifest, indent=2))

def download(rel, url):
    path = RAW / rel
    if path.exists():
        print(f"[skip] {rel} (already there)")
        return True
    tmp = path.with_suffix(path.suffix + ".part")
    try:
        with httpx.stream("GET", url, follow_redirects=True, timeout=120) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for chunk in r.iter_bytes(1 << 20):
                    f.write(chunk)
        tmp.rename(path)
        record(rel, url, path)
        print(f"[ok]   {rel}  ({path.stat().st_size/1e6:.1f} MB)")
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as z:
                z.extractall(path.parent)
        return True
    except Exception as e:
        tmp.unlink(missing_ok=True)
        print(f"[FAIL] {rel}: {e}\n       {url}")
        return False

def fetch_omnipath():
    out = RAW / "omnipath" / "interactions.parquet"
    if out.exists():
        print("[skip] omnipath/interactions.parquet (already there)")
        return True
    try:
        import omnipath as op
        df = op.interactions.OmniPath.get(genesymbols=True)
        df.to_parquet(out)
        record("omnipath/interactions.parquet", "omnipath python client", out)
        print(f"[ok]   omnipath/interactions.parquet  ({len(df):,} rows)")
        return True
    except Exception as e:
        print(f"[FAIL] omnipath: {e}")
        return False

if __name__ == "__main__":
    results = [download(rel, url) for rel, url in SOURCES.items()]
    results.append(fetch_omnipath())
    if not (RAW / "signor").glob("*.tsv") or not list((RAW / "signor").glob("*.tsv")):
        print("[TODO] SIGNOR: download manually, see instructions")
    print(f"\n{sum(results)}/{len(results)} automatic downloads succeeded")
