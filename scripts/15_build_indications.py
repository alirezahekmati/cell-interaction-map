"""Step 6f: indication flags per (disease, target) pair, and indications per drug (Open Targets, CC0).
Reads  data/export/disease_index.json, disease_targets_b64.json (script 14), node_drugs_ot.json (script 13)
Writes data/export/disease_targets4_b64.json : base64 of 4-byte records (node index u16 LE, score*255 u8, flag u8)
         flag 1 = a drug acting on the target is APPROVED for this disease (or a subtype / broader non-area parent)
         flag 2 = no approved drug for it, but one is in trials for it
       data/export/drug_indications.json : {chembl: {"a": [approved diseases], "na": n, "t": [[disease, stage]], "nt": n}}
"""
import base64
import gzip
import json
from collections import defaultdict

import numpy as np
import polars as pl
from otlib import EXP, PROC, ens_maps, sc

pl.Config.set_tbl_rows(30)
pl.Config.set_fmt_str_lengths(60)
pl.Config.set_tbl_width_chars(200)

RANK = {"APPROVAL": 10, "PREAPPROVAL": 9, "PHASE_4": 8, "PHASE_3": 7, "PHASE_2_3": 6, "PHASE_2": 5,
        "PHASE_1_2": 4, "PHASE_1": 3, "EARLY_PHASE_1": 2}
INV = {v: k for k, v in RANK.items()}
TRIAL = [k for k in RANK if k != "APPROVAL"]

index = json.loads((EXP / "disease_index.json").read_text(encoding="utf-8"))
raw = np.frombuffer(base64.b64decode(json.loads((EXP / "disease_targets_b64.json").read_text())), dtype=np.uint8).reshape(-1, 3)
N = len(raw)
assert N == sum(r[2] for r in index), "index and blob disagree"
nd = json.loads((EXP / "node_drugs_ot.json").read_text(encoding="utf-8"))
want = {x[8] for r in nd.values() for x in r.get("d", []) if len(x) > 8}
if not want:
    raise SystemExit("node_drugs_ot.json rows have no drug id: run Block 1 first")

dis = (sc("disease").select("id", "name", "ancestors", "descendants",
                            pl.col("ontology").struct.field("isTherapeuticArea").alias("ta")).collect())
name = dict(zip(dis["id"], dis["name"]))
did = {d: i for i, d in enumerate(dis["id"].to_list())}
ta = dis.filter(pl.col("ta").fill_null(False))["id"].to_list()
print(f"disease terms: {dis.height} | therapeutic-area terms excluded as 'broader parents': {len(ta)}")


def explode_rel(col):
    return dis.select(pl.col("id").alias("D"), pl.col(col).alias("I")).explode(col if False else "I").drop_nulls("I")


anc = explode_rel("ancestors").filter(~pl.col("I").is_in(ta))
desc = explode_rel("descendants")
selfp = dis.select(pl.col("id").alias("D"), pl.col("id").alias("I"))
rel = pl.concat([selfp, desc, anc]).unique().filter(pl.col("D").is_in([r[0] for r in index]))
rel = rel.filter(pl.col("I").is_in(list(did)))
relA = pl.DataFrame({"d": [did[x] for x in rel["D"]], "i": [did[x] for x in rel["I"]]},
                    schema={"d": pl.Int32, "i": pl.Int32})
print("related-disease pairs:", f"{relA.height:,}")

dm = sc("drug_molecule").select("id", "name", "parentId").collect()
canon_of = {i: (p or i) for i, p in zip(dm["id"], dm["parentId"])}
ci = sc("clinical_indication").select("drugId", "diseaseId", "maxClinicalStage").collect()
cd = [canon_of.get(x, x) for x in ci["drugId"]]
bit = np.array([2 if s == "APPROVAL" else (1 if s in TRIAL else 0) for s in ci["maxClinicalStage"]], dtype=np.int8)
dcode = {d: i for i, d in enumerate(sorted(set(cd)))}
cin = (pl.DataFrame({"i": [did.get(x, -1) for x in ci["diseaseId"]], "drug": [dcode[x] for x in cd], "bit": bit},
                    schema={"i": pl.Int32, "drug": pl.Int32, "bit": pl.Int8})
       .filter((pl.col("i") >= 0) & (pl.col("bit") > 0)).group_by("i", "drug").agg(pl.col("bit").max()))
print("drug-disease indication pairs (approved or trial):", f"{cin.height:,}")

ens2acc, acc2idx, _ = ens_maps()
keys = pl.read_parquet(PROC / "layout.parquet")["key"].to_list()
labels = pl.read_parquet(PROC / "layout.parquet")["label"].to_list()
ct = sc("clinical_target").select("drugId", "targetId", "maxClinicalStage").collect().filter(pl.col("maxClinicalStage") != "WITHDRAWAL")
td = set()
for drug, ens in zip(ct["drugId"], ct["targetId"]):
    c = canon_of.get(drug, drug)
    if c in dcode:
        for acc in ens2acc.get(ens, ()):
            if acc in acc2idx:
                td.add((dcode[c], acc2idx[acc]))
mtor = next((i for i, (k, l) in enumerate(zip(keys, labels)) if l == "MTOR" and k.startswith("P:")), None)
rapa = dm.filter(pl.col("name").str.to_lowercase().is_in(["sirolimus", "everolimus", "temsirolimus"]))["id"].to_list()
for x in rapa:
    if mtor is not None and canon_of.get(x, x) in dcode:
        td.add((dcode[canon_of.get(x, x)], mtor))
tdf = pl.DataFrame({"drug": [a for a, _ in td], "idx": [b for _, b in td]}, schema={"drug": pl.Int32, "idx": pl.Int32})
print("drug-target links used:", f"{tdf.height:,}")

x = (relA.join(cin, on="i").group_by("d", "drug").agg(pl.col("bit").max())
     .join(tdf, on="drug").group_by("d", "idx").agg(pl.col("bit").max()))
print("(disease, target) pairs with an approved or trial drug for that disease:", f"{x.height:,}")

dcol = np.repeat(np.array([did.get(r[0], -1) for r in index], dtype=np.int32), [r[2] for r in index])
idx = raw[:, 0].astype(np.int32) | (raw[:, 1].astype(np.int32) << 8)
pairs = pl.DataFrame({"d": dcol, "idx": idx, "row": np.arange(N, dtype=np.int64)},
                     schema={"d": pl.Int32, "idx": pl.Int32, "row": pl.Int64})
res = pairs.join(x, on=["d", "idx"], how="inner")
flag = np.zeros(N, dtype=np.uint8)
flag[res["row"].to_numpy()] = np.where(res["bit"].to_numpy() == 2, 1, 2)
buf = np.empty((N, 4), dtype=np.uint8)
buf[:, :3] = raw
buf[:, 3] = flag
b64 = json.dumps(base64.b64encode(buf.tobytes()).decode("ascii"))
(EXP / "disease_targets4_b64.json").write_text(b64, encoding="utf-8")
print(f"\ndisease_targets4_b64.json: {len(b64):,} bytes raw / {len(gzip.compress(b64.encode())):,} gzip")
print("flags over all pairs: approved-for-disease:", int((flag == 1).sum()), "| trial-for-disease:", int((flag == 2).sum()))

# ---- per-drug indications
ind = defaultdict(dict)
for drug, dz, st in ci.iter_rows():
    c = canon_of.get(drug, drug)
    r = RANK.get(st)
    if c in want and r is not None and r > dict.get(ind[c], dz, -1):
        ind[c][dz] = r
dout = {}
for c, dd in ind.items():
    app = sorted({name.get(z, z) for z, r in dd.items() if r == 10})
    trl = sorted({(name.get(z, z), INV[r]) for z, r in dd.items() if r < 10 and name.get(z, z) not in app},
                 key=lambda t: (-RANK[t[1]], t[0]))
    dout[c] = {"a": app[:8], "na": len(app), "t": [list(t) for t in trl[:6]], "nt": len(trl)}
blob = json.dumps(dout, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "drug_indications.json").write_bytes(blob)
print(f"drug_indications.json: {len(dout)} drugs, {len(blob):,} bytes raw / {len(gzip.compress(blob)):,} gzip")

# ---- sanity: what the disease summary card would say
off, lab = 0, labels
offs = {}
for r in index:
    offs[r[1].lower()] = (r[0], off, r[2])
    off += r[2]
acc_of = {i: k[2:] for i, k in enumerate(keys) if k.startswith("P:")}


def summarize(nm, minscore=0.3):
    if nm not in offs:
        print(f"\n{nm}: not in the index"); return
    d, o, n = offs[nm]
    cats = defaultdict(list)
    for j in range(o, o + n):
        sc_, f, i = buf[j, 2] / 255, buf[j, 3], int(idx[j])
        if sc_ < minscore:
            break
        r = nd.get(acc_of.get(i, ""), {})
        if f == 1:
            c = "approved drug for this disease"
        elif f == 2:
            c = "drug in trials for this disease"
        elif r.get("s") == "APPROVAL":
            c = "approved drug, but for other diseases (repurposing)"
        elif r.get("d"):
            c = "drug in trials, only for other diseases"
        else:
            c = "no clinical-stage drug"
        cats[c].append((lab[i], round(float(sc_), 2)))
    tot = sum(len(v) for v in cats.values())
    print(f"\n--- {nm} ({d}): {n} stored targets, {tot} with score >= {minscore}")
    for c, v in sorted(cats.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(v):4d}  {c}   e.g. {v[:6]}")


summarize("type 2 diabetes mellitus")
summarize("asthma")
