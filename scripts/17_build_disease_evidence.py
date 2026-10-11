"""Step 6h: evidence breakdown per (disease, target) pair (Open Targets, CC0).
Reads data/export/disease_index.json and disease_targets_b64.json (script 14).
Writes data/export/disease_evidence_b64.json: base64 of 4 bytes per pair, in the SAME order as disease_targets4_b64.json:
  8 datatypes x 4 bits (score*15), datatype k at byte k//2, low nibble when k is even.
  order: genetic_association, genetic_literature, somatic_mutation, known_drug, affected_pathway, rna_expression, literature, animal_model
"""
import base64
import gzip
import json
import warnings

import numpy as np
import polars as pl
from otlib import EXP, PROC, ens_maps, sc

warnings.filterwarnings("ignore")
DT = ["genetic_association", "genetic_literature", "somatic_mutation", "known_drug",
      "affected_pathway", "rna_expression", "literature", "animal_model"]

index = json.loads((EXP / "disease_index.json").read_text(encoding="utf-8"))
raw = np.frombuffer(base64.b64decode(json.loads((EXP / "disease_targets_b64.json").read_text())), dtype=np.uint8).reshape(-1, 3)
N = len(raw)
assert N == sum(r[2] for r in index), "index and blob disagree"

ens2acc, acc2idx, _ = ens_maps()
pairs = [(e, acc2idx[x]) for e, accs in ens2acc.items() for x in accs if x in acc2idx]
mp = pl.DataFrame({"ensg": [p[0] for p in pairs], "idx": [p[1] for p in pairs]}, schema={"ensg": pl.String, "idx": pl.Int32})
lut = pl.DataFrame({"diseaseId": [r[0] for r in index], "code": list(range(len(index)))}, schema={"diseaseId": pl.String, "code": pl.Int32})
code = np.repeat(np.arange(len(index), dtype=np.int32), [r[2] for r in index])
idx = (raw[:, 0].astype(np.int32) | (raw[:, 1].astype(np.int32) << 8))
pe = (pl.DataFrame({"code": code, "idx": idx, "row": np.arange(N, dtype=np.int64)},
                   schema={"code": pl.Int32, "idx": pl.Int32, "row": pl.Int64}).join(mp, on="idx"))
print("pairs:", f"{N:,}", "| pair x gene rows:", f"{pe.height:,}")

ev = (sc("association_by_datatype_indirect")
      .select("diseaseId", "targetId", pl.col("aggregationValue").alias("dt"), "associationScore")
      .filter(pl.col("dt").is_in(DT) & pl.col("targetId").is_in(mp["ensg"].to_list()))
      .join(lut.lazy(), on="diseaseId")
      .join(pe.lazy(), left_on=["code", "targetId"], right_on=["code", "ensg"])
      .select("row", "dt", "associationScore").collect())
print("evidence rows matched:", f"{ev.height:,}")
rows = ev["row"].to_numpy()
ks = ev["dt"].replace_strict({n: i for i, n in enumerate(DT)}, return_dtype=pl.Int32).to_numpy()
A = np.zeros((N, 8), dtype=np.float32)
np.maximum.at(A, (rows, ks), ev["associationScore"].to_numpy().astype(np.float32))
Q = np.clip(np.rint(A * 15), 0, 15).astype(np.uint8)
buf = np.empty((N, 4), dtype=np.uint8)
for j in range(4):
    buf[:, j] = Q[:, 2 * j] | (Q[:, 2 * j + 1] << 4)
b64 = json.dumps(base64.b64encode(buf.tobytes()).decode("ascii"))
(EXP / "disease_evidence_b64.json").write_text(b64, encoding="utf-8")
print(f"disease_evidence_b64.json: {len(b64):,} bytes raw / {len(gzip.compress(b64.encode())):,} gzip")
print("pairs with no evidence rows at all:", int((Q.sum(axis=1) == 0).sum()))
print("share of pairs with level >= 3 (score >= 0.17):")
for k, n in enumerate(DT):
    print(f"  {n:22s} {100 * float((Q[:, k] >= 3).mean()):5.1f}%")

lab = pl.read_parquet(PROC / "layout.parquet")["label"].to_list()
off = 0
for r in index:
    if r[0] == "MONDO_0005148":
        print(f"\n{r[1]}: top 10 targets, evidence levels 0-15")
        for j in range(10):
            p = off + j
            print(" ", lab[int(idx[p])], {DT[k][:6]: int(Q[p, k]) for k in range(8) if Q[p, k]})
        break
    off += r[2]
