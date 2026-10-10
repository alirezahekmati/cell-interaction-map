import os
from pathlib import Path
import polars as pl

root = Path(os.environ.get("BIO_ROOT", ".")) / "data/raw/opentargets"
rel = sorted(p.name for p in root.iterdir() if p.is_dir())[-1]
print("release dir:", rel)
for d in sorted(p for p in (root / rel).iterdir() if p.is_dir()):
    files = list(d.rglob("*.parquet"))
    print(f"\n===== {d.name} ({len(files)} files) =====")
    if not files:
        print("no parquet files yet")
        continue
    if list(d.rglob("*.part")):
        print("still downloading, skipped")
        continue
    lf = pl.scan_parquet(str(d / "**" / "*.parquet"))
    print("rows:", lf.select(pl.len()).collect().item())
    for c, t in lf.collect_schema().items():
        print(f"  {c}: {str(t)[:150]}")
    for row in lf.head(1).collect().to_dicts():
        print("sample:", str(row)[:900])
