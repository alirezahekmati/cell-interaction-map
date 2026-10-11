"""Step 6d: node_drugs_ot.json from Open Targets (CC0) plus the Pharos development level.
Keyed by UniProt accession (proteins in the graph only).
  t, f : Pharos level and family
  d    : drugs acting on the protein, one row per parent drug:
         [name, drugType, action, stage, n_ongoing_trials, ongoing_stage, approval_sources]
  s    : best stage reached by any drug (withdrawn and unknown ignored); o : best stage among ongoing trials
  n    : number of drugs
"""
import gzip
import json
import re
from collections import Counter, defaultdict

import polars as pl
from otlib import EXP, PROC, RAW, ens_maps, sc

pl.Config.set_tbl_rows(30)
pl.Config.set_fmt_str_lengths(60)
pl.Config.set_tbl_width_chars(200)

RANK = {"APPROVAL": 10, "PREAPPROVAL": 9, "PHASE_4": 8, "PHASE_3": 7, "PHASE_2_3": 6, "PHASE_2": 5,
        "PHASE_1_2": 4, "PHASE_1": 3, "EARLY_PHASE_1": 2, "IND": 1, "PRECLINICAL": 0, "UNKNOWN": -1, "WITHDRAWAL": -2}
INV = {v: k for k, v in RANK.items()}
ONGOING = ["RECRUITING", "ACTIVE_NOT_RECRUITING", "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION"]

ens2acc, acc2idx, ens2sym = ens_maps()
nodes = pl.read_parquet(PROC / "nodes.parquet").filter(pl.col("kind") == "protein").select("key", "label")
acc2lab = {k[2:]: l for k, l in nodes.iter_rows()}
lab2acc = {l: a for a, l in acc2lab.items()}

ph = pl.read_csv(RAW / "pharos/pharos400_tdls_canonicals.csv", infer_schema_length=0,
                 columns=["uniprot_id", "tdl", "idg_family"])
tdl = {a: (t, f) for a, t, f in ph.iter_rows() if a in acc2lab and t}

dm = sc("drug_molecule").select("id", "name", "drugType", "parentId", "maximumClinicalStage").collect()
canon_of = {i: (p or i) for i, p in zip(dm["id"], dm["parentId"])}
name_of = dict(zip(dm["id"], dm["name"]))
type_of = dict(zip(dm["id"], dm["drugType"]))
stage_of = dict(zip(dm["id"], dm["maximumClinicalStage"]))


def canon(i):
    return canon_of.get(i, i)


def disp(n):
    n = (n or "").strip()
    return n if re.search(r"\d", n) else n.lower()


crx = (sc("clinical_report").select("clinicalStage", "origin", "source", "trialOverallStatus", "drugs")
       .explode("drugs")
       .with_columns(drug=pl.col("drugs").struct.field("drugId")).drop("drugs")
       .filter(pl.col("drug").is_not_null())
       .with_columns(drug=pl.col("drug").replace(canon_of)))
print("== approval reports by origin and source ==")
print(crx.filter(pl.col("clinicalStage") == "APPROVAL").group_by("origin", "source").len()
      .sort("len", descending=True).head(15).collect())

og = (crx.filter((pl.col("origin") == "CLINICAL_TRIAL") & pl.col("trialOverallStatus").is_in(ONGOING))
      .with_columns(rank=pl.col("clinicalStage").replace_strict(RANK, default=-1, return_dtype=pl.Int32))
      .group_by("drug").agg(pl.len().alias("n"), pl.col("rank").max().alias("rmax")).collect())
ongoing = {d: (int(n), INV.get(int(r), "UNKNOWN")) for d, n, r in og.iter_rows()}
ap = (crx.filter(pl.col("clinicalStage") == "APPROVAL").group_by("drug")
      .agg(pl.col("source").drop_nulls().unique().sort().alias("srcs")).collect())
approvals = dict(ap.iter_rows())

ct = sc("clinical_target").select("drugId", "targetId", "maxClinicalStage").collect()
moa = (sc("drug_mechanism_of_action").select("chemblIds", "targets", "actionType").collect()
       .explode("chemblIds").explode("targets").drop_nulls(["chemblIds", "targets"]))
acts, moa_t = defaultdict(set), defaultdict(set)
for c, t, a in moa.iter_rows():
    moa_t[canon(c)].add(t)
    if a:
        acts[(canon(c), t)].add(a.lower())

per = defaultdict(dict)
for drug, ens, stage in ct.iter_rows():
    cd = canon(drug)
    for acc in ens2acc.get(ens, ()):
        if acc not in acc2lab:
            continue
        rec = per[acc].setdefault(cd, {"rank": -99, "stage": "UNKNOWN", "ens": ens, "w": 0})
        if stage == "WITHDRAWAL":
            rec["w"] = 1
            continue
        r = RANK.get(stage, -1)
        if r > rec["rank"]:
            rec["rank"], rec["stage"] = r, stage

# ChEMBL links the rapalogs to FKBP1A only; they inhibit mTORC1 through the FKBP12 complex (curated by hand)
MANUAL = {"MTOR": ["sirolimus", "everolimus", "temsirolimus"]}
dm_l = dm.with_columns(nl=pl.col("name").str.to_lowercase())
for g, names in MANUAL.items():
    acc = lab2acc.get(g)
    for cd in {canon(i) for i in dm_l.filter(pl.col("nl").is_in(names))["id"].to_list()}:
        if acc and cd not in per[acc]:
            st = stage_of.get(cd) or "UNKNOWN"
            per[acc][cd] = {"rank": RANK.get(st, -1), "stage": st, "ens": None, "w": 0, "manual": True}
out, n_rows = {}, 0
for acc in acc2lab:
    rec = {}
    if acc in tdl:
        rec["t"] = tdl[acc][0]
        if tdl[acc][1]:
            rec["f"] = tdl[acc][1]
    if acc in per:
        rows = []
        for cd, r in per[acc].items():
            on_n, on_stage = ongoing.get(cd, (0, ""))
            act = " / ".join(sorted(acts.get((cd, r["ens"]), ())))
            srcs = ",".join((approvals.get(cd) or [])[:5])
            st, rk = r["stage"], r["rank"]
            if rk == -99:
                st = stage_of.get(cd) or "UNKNOWN"
                rk = RANK.get(st, -1)
            if r.get("manual"):
                act, srcs = "allosteric inhibitor (acts via FKBP12)", (srcs + ",manual").strip(",")
            rows.append([disp(name_of.get(cd, cd)), type_of.get(cd) or "", act, st, on_n, on_stage, srcs, (2 if r["rank"] == -99 else r["w"]), rk, cd])
        rows.sort(key=lambda x: (-x[8], -RANK.get(x[5], -9), x[0]))
        rec["n"] = len(rows)
        n_rows += len(rows)
        rec["d"] = [x[:8] + [x[9]] for x in rows]
        known = [x[8] for x in rows if x[8] >= 0 and x[7] != 2]
        if known:
            rec["s"] = INV[max(known)]
        og_r = [RANK[x[5]] for x in rows if x[4] and x[5] in RANK]
        if og_r:
            rec["o"] = INV[max(og_r)]
    if rec:
        out[acc] = rec

blob = json.dumps(out, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "node_drugs_ot.json").write_bytes(blob)
print(f"\nnode_drugs_ot.json: {len(out)} proteins, {n_rows} drug rows, {len(blob):,} bytes raw / {len(gzip.compress(blob)):,} gzip")
has = [r for r in out.values() if "d" in r]
appr = {a for a, r in out.items() if r.get("s") == "APPROVAL"}
act_ = [a for a, r in out.items() if "o" in r]
only = [a for a in act_ if out[a].get("s") != "APPROVAL"]
print("proteins with any drug:", len(has), "| best stage APPROVAL:", len(appr),
      "| with an ongoing trial:", len(act_), "| ongoing trial but not approved:", len(only))
print("best stage counts:", dict(Counter(r.get("s") for r in has).most_common()))
print("drug types in rows:", dict(Counter(x[1] for r in has for x in r["d"]).most_common()))

oldf = EXP / "node_drugs.json"
if oldf.exists():
    old = json.loads(oldf.read_text(encoding="utf-8"))
    dc = {a for a, r in old.items() if any(x[2] != "-" for x in r.get("d", []))}
    print(f"\nDrugCentral approved-drug targets: {len(dc)} | Open Targets APPROVAL: {len(appr)} | both: {len(dc & appr)}")
    print("DrugCentral-approved but not OT-approved (sample):", [acc2lab[a] for a in sorted(dc - appr)][:25])
    print("OT-approved but not DrugCentral (sample):", [acc2lab[a] for a in sorted(appr - dc)][:25])

print("\n== samples ==")
for g in ["MTOR", "AKT1", "EGFR", "PRKAA1", "TSC2", "TP53", "KCNJ11", "GLP1R", "PPARG", "DPP4", "GPD2", "NDUFS2"]:
    a = lab2acc.get(g)
    r = out.get(a, {})
    print(f"\n--- {g} ({a}) level={r.get('t')} best={r.get('s')} ongoing={r.get('o')} drugs={r.get('n', 0)}")
    for x in r.get("d", [])[:7]:
        print("   ", x)

print("\n== where do these drugs link? ==")
for nm in ["sirolimus", "everolimus", "temsirolimus", "ridaforolimus", "perhexiline", "metformin"]:
    for r in dm.filter(pl.col("name").str.to_lowercase().str.starts_with(nm)).head(3).iter_rows(named=True):
        cd = canon(r["id"])
        tg_ct = sorted({ens2sym.get(t, t) + ":" + s for d, t, s in ct.iter_rows() if canon(d) == cd})
        tg_moa = sorted({ens2sym.get(t, t) for t in moa_t.get(cd, ())})
        print(f"  {r['id']} parent={r['parentId']} {r['name']} [{r['drugType']}]")
        print("     clinical_target ->", tg_ct[:8], "| mechanism targets ->", tg_moa[:8])
