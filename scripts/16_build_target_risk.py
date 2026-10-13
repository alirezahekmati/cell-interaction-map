"""Step 6g: target tractability, safety, tool compounds and essentiality (Open Targets, CC0).
Writes data/export/node_tract.json : {acc: [sm, ab, pr, oc, ess, n_safety, n_probes, n_hq_probes, pct_dep, n_lines]}  (startup)
       data/export/node_risk.json  : {acc: {"s": [[event, sources]], "ns": n, "p": [[name, hq, url]], "t": [[tissue, pct, n]]}}  (lazy)
levels: 0 none, 1 predicted (lower confidence), 2 predicted (high confidence), 3 phase 1, 4 advanced clinical, 5 approved drug
ess: 0 no data, 1 not common essential, 2 common essential. pct_dep = % of screened cancer cell lines with gene effect <= -0.5
"""
import gzip
import json
import warnings
from collections import Counter, defaultdict

import numpy as np
import polars as pl
from otlib import EXP, OT, PROC, ens_maps, sc

warnings.filterwarnings("ignore")
pl.Config.set_tbl_rows(80)
pl.Config.set_fmt_str_lengths(60)
pl.Config.set_tbl_width_chars(200)

MODS = ["SM", "AB", "PR", "OC"]
CLIN = {"Approved Drug": 5, "Advanced Clinical": 4, "Phase 1 Clinical": 3}
DEP = -0.5


def level(cat):
    if cat in CLIN:
        return CLIN[cat]
    c = cat.lower()
    return 2 if ("high-quality" in c or "high conf" in c) else 1


ens2acc, acc2idx, ens2sym = ens_maps()
ens_list = list(ens2acc)
nodes = pl.read_parquet(PROC / "nodes.parquet").filter(pl.col("kind") == "protein").select("key", "label")
acc2lab = {k[2:]: l for k, l in nodes.iter_rows()}
lab2acc = {l: a for a, l in acc2lab.items()}

# ---- tractability
tr = (sc("target_tractability").filter(pl.col("value")).filter(pl.col("targetId").is_in(ens_list))
      .select("targetId", "modality", "category").collect())
print("== tractability categories that are TRUE, for targets in your map ==")
print(tr.group_by("modality", "category").len().sort(["modality", "len"], descending=[False, True]))
lv = defaultdict(lambda: [0, 0, 0, 0])
for e, m, c in tr.iter_rows():
    if m in MODS:
        k = MODS.index(m)
        lv[e][k] = max(lv[e][k], level(c))

# ---- safety events
sf = (sc("target_safety_event").filter(pl.col("targetId").is_in(ens_list))
      .select("targetId", "event", "datasource").collect())
ev_by = defaultdict(lambda: defaultdict(set))
for e, ev, ds in sf.iter_rows():
    if ev:
        ev_by[e][ev].add(ds or "?")

# ---- tool compounds
pb = (sc("chemical_probe").filter(pl.col("targetId").is_in(ens_list))
      .select("targetId", "drugFromSource", "isHighQuality", "probesDrugsScore", "urls").collect())
pbs = defaultdict(dict)
for e, nm, hq, score, urls in pb.iter_rows():
    if not nm:
        continue
    url = ""
    for u in (urls or []):
        x = u.get("url") if isinstance(u, dict) else None
        if x and x.startswith(("http://", "https://")):
            url = x
            break
    s = score if (score is not None and score == score) else 0.0
    old = pbs[e].get(nm)
    if old is None or (bool(hq), s) > (old[0], old[1]):
        pbs[e][nm] = (bool(hq), s, url)

# ---- essentiality (one file at a time: the nested screens are large)
parts_b, parts_q = [], []
for f in sorted((OT / "target_essentiality").rglob("*.parquet")):
    lf = pl.scan_parquet(str(f)).filter(pl.col("targetId").is_in(ens_list))
    parts_b.append(lf.select("targetId", "isEssential").collect())
    parts_q.append(
        lf.select("targetId", pl.col("depMapEssentiality").alias("t")).explode("t")
        .select("targetId", pl.col("t").struct.field("tissueName").alias("tissue"), pl.col("t").struct.field("screens").alias("s"))
        .explode("s")
        .select("targetId", "tissue", pl.col("s").struct.field("geneEffect").alias("ge"))
        .drop_nulls("ge")
        .group_by("targetId", "tissue").agg((pl.col("ge") <= DEP).sum().alias("dep"), pl.len().alias("n"))
        .collect())
    print("essentiality file done:", f.name, flush=True)
base, q = pl.concat(parts_b), pl.concat(parts_q)
ess = {e: bool(v) for e, v in base.iter_rows() if v is not None}
tot = q.group_by("targetId").agg(pl.col("dep").sum().alias("dep"), pl.col("n").sum().alias("n"))
pct = {e: (int(round(100 * d / n)), int(n)) for e, d, n in tot.iter_rows() if n}
tis = defaultdict(list)
q2 = (q.filter((pl.col("n") >= 5) & pl.col("tissue").is_not_null())
      .with_columns(p=(pl.col("dep") / pl.col("n") * 100).round())
      .sort(["targetId", "p", "n"], descending=[False, True, True]))
for e, t_, p_, n_ in q2.select("targetId", "tissue", "p", "n").iter_rows():
    if len(tis[e]) < 3:
        tis[e].append([t_, int(p_), int(n_)])

# ---- one record per accession
T = {}
for e, accs in ens2acc.items():
    for acc in accs:
        if acc not in acc2lab:
            continue
        rec = T.setdefault(acc, {"lv": [0, 0, 0, 0], "ess": 0, "ev": {}, "pb": {}, "pct": 0, "n": 0, "tis": []})
        rec["lv"] = [max(a, b) for a, b in zip(rec["lv"], lv.get(e, [0, 0, 0, 0]))]
        if e in ess:
            rec["ess"] = max(rec["ess"], 2 if ess[e] else 1)
        for k, srcs in ev_by.get(e, {}).items():
            rec["ev"].setdefault(k, set()).update(srcs)
        for k, v in pbs.get(e, {}).items():
            if k not in rec["pb"] or (v[0], v[1]) > (rec["pb"][k][0], rec["pb"][k][1]):
                rec["pb"][k] = v
        if e in pct and pct[e][1] > rec["n"]:
            rec["pct"], rec["n"] = pct[e]
            rec["tis"] = tis.get(e, [])

tract_out, risk_out = {}, {}
for acc, rec in T.items():
    nsf, npb = sum(1 for v in rec["ev"].values() if v - {"ToxCast", "AOP-Wiki"}), len(rec["pb"])
    nhq = sum(1 for v in rec["pb"].values() if v[0])
    if not (any(rec["lv"]) or rec["ess"] or nsf or npb):
        continue
    tract_out[acc] = rec["lv"] + [rec["ess"], nsf, npb, nhq, rec["pct"], rec["n"]]
    r = {}
    if nsf:
        evs = sorted(((k, v) for k, v in rec["ev"].items() if v - {"ToxCast", "AOP-Wiki"}), key=lambda kv: (-len([s for s in kv[1] if s != "ToxCast"]), -len(kv[1]), kv[0]))
        r["s"] = [[k, ", ".join(sorted(v))] for k, v in evs[:8]]
        r["ns"] = nsf
    if npb:
        ps = sorted(rec["pb"].items(), key=lambda kv: (-int(kv[1][0]), -kv[1][1], kv[0]))
        r["p"] = [[k, int(v[0]), v[2]] for k, v in ps[:6]]
    if rec["tis"] and rec["pct"] >= 10:
        r["t"] = rec["tis"]
    if r:
        risk_out[acc] = r

b1 = json.dumps(tract_out, separators=(",", ":")).encode()
b2 = json.dumps(risk_out, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "node_tract.json").write_bytes(b1)
(EXP / "node_risk.json").write_bytes(b2)
print(f"\nnode_tract.json: {len(tract_out)} proteins, {len(b1):,} bytes raw / {len(gzip.compress(b1)):,} gzip")
print(f"node_risk.json: {len(risk_out)} proteins, {len(b2):,} bytes raw / {len(gzip.compress(b2)):,} gzip")

print("\n== tractability levels over your proteins (0..5) ==")
for k, m in enumerate(MODS):
    c = Counter(r[k] for r in tract_out.values())
    print(f"  {m}: " + ", ".join(f"{lvl}:{c.get(lvl, 0)}" for lvl in range(6)))
print("\n== essentiality ==")
print("targets with essentiality data:", sum(1 for r in tract_out.values() if r[9]),
      "| common essential:", sum(1 for r in tract_out.values() if r[4] == 2))
pv = [r[8] for r in tract_out.values() if r[9]]
print("percent of cell lines depending (p10/p50/p90/p99):", np.percentile(pv, [10, 50, 90, 99]).round(1).tolist())
for flag, nm in [(2, "common essential"), (1, "not common essential")]:
    xs = [r[8] for r in tract_out.values() if r[4] == flag and r[9]]
    print(f"  mean % of lines depending, {nm}: {np.mean(xs):.1f}  (n={len(xs)})" if xs else f"  {nm}: none")
print("\n== safety ==")
print("proteins with safety events:", sum(1 for r in tract_out.values() if r[5]))
print("top events:", Counter(k for rec in T.values() for k in rec["ev"]).most_common(8))
print("by source:", sf.group_by("datasource").len().sort("len", descending=True).rows())
print("\n== tool compounds ==")
print("proteins with a probe:", sum(1 for r in tract_out.values() if r[6]), "| with a high-quality probe:", sum(1 for r in tract_out.values() if r[7]))

nd = json.loads((EXP / "node_drugs_ot.json").read_text(encoding="utf-8"))
undrugged = {a for a in acc2lab if not nd.get(a, {}).get("d")}
tractable = {a for a, r in tract_out.items() if max(r[0], r[1]) >= 2}
print("\n== triage numbers ==")
print("no clinical-stage drug:", len(undrugged), "| ... and tractable (SM or AB predicted high or better):", len(undrugged & tractable))
print("... and with a high-quality tool compound:", len([a for a in undrugged if tract_out.get(a, [0] * 8)[7]]))
print("... and common essential (likely toxic to hit):", len([a for a in undrugged & tractable if tract_out[a][4] == 2]))

print("\n== samples ==")
for g in ["EGFR", "MTOR", "TSC2", "TP53", "STK11", "GCK", "KRAS", "PLK1"]:
    a = lab2acc.get(g)
    print(f"\n--- {g} ({a})\n  tract: {tract_out.get(a)}\n  risk: {json.dumps(risk_out.get(a), ensure_ascii=False)[:600]}")
