"""Step 6i: drug_names.json {chembl id: 'trade names and synonyms, | separated'} for the drugs shown in the viewer
(Open Targets, CC0). Used only for searching."""
import gzip
import json
import warnings

import polars as pl
from otlib import EXP, sc

warnings.filterwarnings("ignore")
nd = json.loads((EXP / "node_drugs_ot.json").read_text(encoding="utf-8"))
want, primary = set(), {}
for r in nd.values():
    for x in r.get("d", []):
        if len(x) > 8:
            want.add(x[8])
            primary[x[8]] = x[0].lower()
dm = sc("drug_molecule").select("id", "name", "parentId", "tradeNames", "synonyms").collect()
names = {}
for r in dm.iter_rows(named=True):
    c = r["parentId"] or r["id"]
    if c not in want:
        continue
    s = names.setdefault(c, [])
    cand = [r["name"]]
    for key in ("tradeNames", "synonyms"):
        for it in (r[key] or []):
            if isinstance(it, dict):
                cand.append(it.get("label"))
    for x in cand:
        x = (x or "").strip().lower()
        if x and len(x) <= 40 and x != primary[c] and x not in s:
            s.append(x)
out = {c: "|".join(v[:14]) for c, v in names.items() if v}
blob = json.dumps(out, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "drug_names.json").write_bytes(blob)
print(f"drug_names.json: {len(out)} of {len(want)} drugs have extra names, {len(blob):,} bytes raw / {len(gzip.compress(blob)):,} gzip")
look = {v: k for k, v in primary.items()}
for g in ["sirolimus", "semaglutide", "metformin", "osimertinib", "capivasertib", "tocilizumab"]:
    print(" ", g, "->", (out.get(look.get(g)) or "(no extra names)")[:160])
