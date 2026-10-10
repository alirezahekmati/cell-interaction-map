"""Step 6e: disease -> target scores for the viewer (Open Targets, CC0).
Writes data/export/disease_index.json  [[id, name, n_targets, synonyms], ...] (largest first)
       data/export/disease_targets_b64.json  base64 of 3-byte records: node index (u16 LE), score*255 (u8),
       per disease in index order, highest score first.
  uv run python scripts/14_build_ot_diseases.py [--min-score 0.1] [--cap 2000]
"""
import argparse
import base64
import gzip
import json

import numpy as np
import polars as pl
from otlib import EXP, PROC, ens_maps, sc

ap = argparse.ArgumentParser()
ap.add_argument("--min-score", type=float, default=0.1)
ap.add_argument("--cap", type=int, default=2000)
a = ap.parse_args()
pl.Config.set_tbl_rows(30)
pl.Config.set_tbl_width_chars(200)

ens2acc, acc2idx, ens2sym = ens_maps()
pairs = [(e, acc2idx[x]) for e, accs in ens2acc.items() for x in accs if x in acc2idx]
mp = pl.DataFrame({"ensg": [p[0] for p in pairs], "idx": [p[1] for p in pairs]},
                  schema={"ensg": pl.String, "idx": pl.Int64})
print("Ensembl genes mapped to graph nodes:", mp["ensg"].n_unique(), "| node indexes:", mp["idx"].n_unique())
assert mp["idx"].max() < 65536

LOW = min(0.05, a.min_score)
df = (sc("association_overall_indirect").select("diseaseId", "targetId", "associationScore")
      .filter((pl.col("associationScore") >= LOW) & pl.col("targetId").is_in(mp["ensg"].to_list()))
      .join(mp.lazy(), left_on="targetId", right_on="ensg")
      .group_by("diseaseId", "idx").agg(pl.col("associationScore").max())
      .collect())
print("pairs at the lowest cutoff considered:", f"{df.height:,}")
print("size of alternatives (3 bytes per pair, before base64):")
for m, cap in [(0.05, 3000), (0.1, 2000), (0.1, 3000), (0.2, 2000)]:
    n = (df.filter(pl.col("associationScore") >= m).group_by("diseaseId").len()
         .select(pl.col("len").clip(upper_bound=cap).sum()).item())
    print(f"  min {m}, cap {cap}: {n:,} pairs, about {n * 3 / 1e6:.1f} MB")

df = (df.filter(pl.col("associationScore") >= a.min_score)
      .sort(["diseaseId", "associationScore", "idx"], descending=[False, True, False]))
df = df.group_by("diseaseId", maintain_order=True).head(a.cap)
cnt = (df.group_by("diseaseId", maintain_order=True).len()
       .sort(["len", "diseaseId"], descending=[True, False]).with_row_index("ord"))
df = (df.join(cnt.select("diseaseId", "ord"), on="diseaseId")
      .sort(["ord", "associationScore", "idx"], descending=[False, True, False]))

dis = sc("disease").select("id", "name", "exactSynonyms").collect()
name, syn = dict(zip(dis["id"], dis["name"])), dict(zip(dis["id"], dis["exactSynonyms"].to_list()))
index = []
for d, n in zip(cnt["diseaseId"], cnt["len"]):
    s = [x for x in (syn.get(d) or []) if x][:6]
    index.append([d, name.get(d) or d, int(n), "|".join(s)[:300]])

N = df.height
idx16 = df["idx"].to_numpy().astype("<u2")
sc8 = np.clip(np.rint(df["associationScore"].to_numpy() * 255), 1, 255).astype("u1")
buf = np.empty((N, 3), dtype=np.uint8)
buf[:, :2] = idx16.view(np.uint8).reshape(-1, 2)
buf[:, 2] = sc8
blob = buf.tobytes()

(EXP / "disease_index.json").write_text(json.dumps(index, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
b64 = json.dumps(base64.b64encode(blob).decode("ascii"))
(EXP / "disease_targets_b64.json").write_text(b64, encoding="utf-8")
ij = json.dumps(index, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
print(f"\nmin score {a.min_score}, cap {a.cap}: {N:,} pairs for {len(index):,} diseases")
print(f"disease_index.json: {len(ij):,} bytes raw / {len(gzip.compress(ij)):,} gzip")
print(f"disease_targets_b64.json: {len(b64):,} bytes raw / {len(gzip.compress(b64.encode())):,} gzip")
print("diseases that hit the cap:", sum(1 for r in index if r[2] >= a.cap),
      "| targets per disease: median", int(np.median([r[2] for r in index])), "max", max(r[2] for r in index))

raw = np.frombuffer(base64.b64decode(json.loads(b64)), dtype=np.uint8).reshape(-1, 3)
assert len(raw) == sum(r[2] for r in index) == N
label = pl.read_parquet(PROC / "layout.parquet")["label"].to_list()
off = 0
for r in index:
    if r[0] == "MONDO_0005148":
        top = [(label[int(raw[off + k, 0]) | (int(raw[off + k, 1]) << 8)], round(int(raw[off + k, 2]) / 255, 2)) for k in range(10)]
        print(f"\nround trip OK. {r[1]} ({r[0]}): {r[2]} targets; top 10 decoded from the blob:\n  {top}")
        break
    off += r[2]
else:
    print("round trip OK (type 2 diabetes mellitus not in the index)")
