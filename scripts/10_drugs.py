"""Step 6c: node_drugs.json, keyed by UniProt accession (proteins in the graph only).
t = Pharos target development level, f = IDG family,
d = [[drug, action, agencies]] mechanism-of-action drugs from DrugCentral (CC BY-SA 4.0);
agencies: F=FDA E=EMA P=PMDA ('-' = in DrugCentral but on none of those three lists); n = total drugs.
"""
import gzip
import io
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
RAW, PROC, EXP = ROOT / "data/raw", ROOT / "data/processed", ROOT / "data/export"
EXP.mkdir(parents=True, exist_ok=True)
CAP = 25


def unq(df):
    df = df.rename({c: c.strip('"') for c in df.columns})
    return df.with_columns(pl.all().cast(pl.String).str.strip_chars('"'))


nodes = pl.read_parquet(PROC / "nodes.parquet").filter(pl.col("kind") == "protein").select("key", "label")
acc2lab = {k[2:]: l for k, l in nodes.iter_rows()}

ph = pl.read_csv(RAW / "pharos/pharos400_tdls_canonicals.csv", infer_schema_length=0,
                 columns=["uniprot_id", "tdl", "idg_family"])
tdl = {a: (t, f) for a, t, f in ph.iter_rows() if a in acc2lab and t}

agency = defaultdict(set)
for tag, f in [("F", "FDA"), ("E", "EMA"), ("P", "PMDA")]:
    a = pl.read_csv(RAW / f"drugs/{f}_Approved.csv", has_header=False, new_columns=["sid", "name"],
                    infer_schema_length=0, encoding="utf8-lossy")
    for sid in a["sid"].to_list():
        agency[str(sid).strip()].add(tag)
print("approved-list sizes:", {t: sum(t in v for v in agency.values()) for t in "FEP"})

raw = gzip.open(RAW / "drugs/drug.target.interaction.tsv.gz", "rb").read()
dt = unq(pl.read_csv(io.BytesIO(raw), separator="\t", infer_schema_length=0, quote_char=None,
                     truncate_ragged_lines=True, encoding="utf8-lossy"))
print("rows:", dt.height, "| rows with MOA_SOURCE == '1' (possible column shift):",
      dt.filter(pl.col("MOA_SOURCE") == "1").height)
print(dt.filter(pl.col("MOA_SOURCE") == "1").select("DRUG_NAME", "TARGET_NAME", "MOA", "MOA_SOURCE", "ACTION_TYPE").head(3))
moa = dt.filter((pl.col("MOA") == "1") & (pl.col("ORGANISM") == "Homo sapiens") & pl.col("ACCESSION").is_not_null())
print("mechanism-of-action rows (human):", moa.height, "| distinct drugs:", moa["STRUCT_ID"].n_unique())
print("  distinct drugs on none of FDA/EMA/PMDA:",
      sum(1 for s in moa["STRUCT_ID"].unique().to_list() if not agency.get(str(s).strip())))

ACC = re.compile(r"^[A-Z][0-9A-Z]{5,9}$")
per = defaultdict(dict)
for drug, sid, accs, act in moa.select("DRUG_NAME", "STRUCT_ID", "ACCESSION", "ACTION_TYPE").iter_rows():
    ag = set(agency.get(str(sid).strip(), ()))
    for acc in re.split(r"[|;,\s]+", accs or ""):
        if ACC.match(acc) and acc in acc2lab:
            key = (drug, (act or "").lower())
            old = set(per[acc].get(key, ""))
            per[acc][key] = "".join(sorted(old | ag, key="FEP".index))

out = {}
for acc in acc2lab:
    rec = {}
    if acc in tdl:
        rec["t"] = tdl[acc][0]
        if tdl[acc][1]:
            rec["f"] = tdl[acc][1]
    if acc in per:
        items = sorted(per[acc].items(), key=lambda kv: (-len(kv[1]), kv[0][0].lower()))
        rec["n"] = len(items)
        rec["d"] = [[d, a, ag or "-"] for (d, a), ag in items[:CAP]]
    if rec:
        out[acc] = rec

blob = json.dumps(out, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "node_drugs.json").write_bytes(blob)
import gzip as _gz
print(f"\nnode_drugs.json: {len(out)} proteins, {len(blob):,} bytes raw / {len(_gz.compress(blob)):,} gzip")
print("level counts:", dict(Counter(r.get("t") for r in out.values())))
print("proteins with a drug list:", sum("d" in r for r in out.values()))

tclin = {a for a, r in out.items() if r.get("t") == "Tclin"}
withd = {a for a, r in out.items() if "d" in r}
print(f"\nTclin: {len(tclin)} | with drug list: {len(withd)} | both: {len(tclin & withd)}")
print("Tclin but NO drug list (sample):", [acc2lab[a] for a in sorted(tclin - withd)][:15])
print("drug list but NOT Tclin (sample):", [(acc2lab[a], out[a].get("t")) for a in sorted(withd - tclin)][:15])
big = sorted(withd, key=lambda a: -out[a]["n"])[:6]
print("largest drug lists:", [(acc2lab[a], out[a]["n"]) for a in big])

lab2acc = {l: a for a, l in acc2lab.items()}
print("\n== samples ==")
for g in ["MTOR", "EGFR", "AKT1", "PRKAA1", "TSC2", "TP53", "STK11"]:
    a = lab2acc.get(g)
    print(f"\n--- {g} ({a})")
    print(json.dumps(out.get(a), ensure_ascii=False)[:700] if a else "not in graph")
