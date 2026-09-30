"""Step 6: export for the viewer -> data/export/
nodes.json, edges.bin (20-byte LE records), edge_details.json (columnar), regions.json, meta.json
"""
import gzip
import json
import os
import re
import shutil
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
PROC, EXP = ROOT / "data/processed", ROOT / "data/export"
EXP.mkdir(parents=True, exist_ok=True)
pl.Config.set_tbl_rows(40)
pl.Config.set_tbl_width_chars(200)

CLS = ["regulation", "binding", "membership", "substrate_product", "transport"]
NET = ["unsigned", "stim", "inhib", "mixed"]
DT = np.dtype([("src", "<u4"), ("tgt", "<u4"), ("cls", "u1"), ("net", "u1"), ("directed", "u1"),
               ("pad", "u1"), ("n_ev", "<u2"), ("n_res", "<u2"), ("score", "<f4")])
assert DT.itemsize == 20
METALS = ("iron", "copper", "zinc", "calcium", "magnesium", "manganese", "cobalt", "nickel",
          "sodium", "potassium", "chromium", "mercury", "lithium")
LABEL_FIX = {"C:l-glutamate": "L-glutamate", "C:l-aspartate": "L-aspartate",
             "C:l-glutamine": "L-glutamine", "C:udp-alpha-d-galactose": "UDP-D-galactose"}
NOISE = re.compile(r"^(calcium|sodium|potassium|chloride|hydron|messenger rna|iron-sulfur)", re.I)
NET_COLORS = {"stim": "#2ecc71", "inhib": "#e74c3c", "mixed": "#f5a623", "unsigned": "#8a8f98"}
REGION_COLORS = {"extracellular": "#9aa5b1", "plasma_membrane": "#e6b800", "cytosol": "#4dabf7",
                 "nucleus": "#9775fa", "mitochondrion": "#ff922b", "er": "#20c997",
                 "golgi": "#f06595", "lysosome": "#fa5252", "vesicle": "#94d82d",
                 "peroxisome": "#15aabf"}


def hdr(t):
    print(f"\n===== {t} =====")


def clean(key, label, kind):
    if key in LABEL_FIX:
        return LABEL_FIX[key]
    if kind not in ("metabolite", "chemical") or label.lower().startswith(METALS):
        return label
    s = re.sub(r"\s*zwitterion\s*$", "", label)
    s = re.sub(r"\(\d+[+-]\)\s*$", "", s).strip()
    return s or label


def items(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [str(x) for x in v if x not in (None, "")]
    return [x.strip() for x in re.split(r"[;,|]", str(v)) if x.strip()]


lay = pl.read_parquet(PROC / "layout.parquet")
nd = pl.read_parquet(PROC / "nodes.parquet").select(
    "key", "uniprot", "chebi", "deg_all", "deg_out_reg", "deg_in_reg", "n_members")
ed = pl.read_parquet(PROC / "edges.parquet")
n = lay.join(nd, on="key", how="left")          # node index = row position of this frame
assert n.height == lay.height

# ---------------------------------------------------------------- nodes
rows, hid, changed = [], [], []
for r in n.iter_rows(named=True):
    lab = clean(r["key"], r["label"], r["kind"])
    if lab != r["label"]:
        changed.append((r["label"], lab))
    h = 0
    if r["kind"] == "chemical":
        h = 1
    elif r["kind"] not in ("protein", "complex", "family") and NOISE.match(lab):
        h = 2
    hid.append(h)
    rows.append({"k": r["key"], "l": lab, "t": r["kind"], "r": r["region"],
                 "x": round(float(r["x"]), 1), "y": round(float(r["y"]), 1), "z": round(float(r["z"]), 1),
                 "d": int(r["deg_all"] or 0), "dr": int((r["deg_out_reg"] or 0) + (r["deg_in_reg"] or 0)),
                 "h": h, "u": r["uniprot"] or "", "c": r["chebi"] or "", "m": int(r["n_members"] or 0)})
hid = np.array(hid, dtype=np.uint8)
(EXP / "nodes.json").write_text(json.dumps(rows, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")

# ---------------------------------------------------------------- edges
idx = {k: i for i, k in enumerate(n["key"])}
known = pl.Series(list(idx))
assert ed.filter(~pl.col("src").is_in(known) | ~pl.col("tgt").is_in(known)).height == 0, "edge endpoint missing"
rec = np.zeros(ed.height, dtype=DT)
rec["src"] = [idx[k] for k in ed["src"]]
rec["tgt"] = [idx[k] for k in ed["tgt"]]
rec["cls"] = [CLS.index(c) for c in ed["cls"]]
rec["net"] = [NET.index(c) for c in ed["net"]]
rec["directed"] = ed["directed"].cast(pl.Int8).fill_null(0).to_numpy()
rec["n_ev"] = np.minimum(ed["n_evidence"].fill_null(0).to_numpy(), 65535).astype("uint16")
rec["n_res"] = np.minimum(ed["n_resources"].fill_null(0).to_numpy(), 65535).astype("uint16")
rec["score"] = ed["max_score"].cast(pl.Float32).fill_null(0).to_numpy()
rec.tofile(EXP / "edges.bin")

det = {"res": [], "mech": [], "pmids": [], "npm": [], "ns": []}
for rs, mech, pm, npm, a, b, c in ed.select("resources", "mechanisms", "pmids", "n_pmids",
                                            "n_stim", "n_inhib", "n_unsigned").iter_rows():
    det["res"].append(";".join(items(rs)))
    det["mech"].append(";".join(items(mech)))
    det["pmids"].append(";".join(items(pm)[:10]))
    det["npm"].append(int(npm or 0))
    det["ns"].append([int(a or 0), int(b or 0), int(c or 0)])
(EXP / "edge_details.json").write_text(json.dumps(det, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")

shutil.copy(PROC / "regions.json", EXP / "regions.json")
meta = {"n_nodes": len(rows), "n_edges": int(ed.height), "cls": CLS, "net": NET,
        "edge_record": {"bytes": 20, "little_endian": True,
                        "fields": [["src", "u32"], ["tgt", "u32"], ["cls", "u8"], ["net", "u8"],
                                   ["directed", "u8"], ["pad", "u8"], ["n_ev", "u16"],
                                   ["n_res", "u16"], ["score", "f32"]]},
        "net_colors": NET_COLORS, "region_colors": REGION_COLORS,
        "hidden_codes": {"1": "chemical", "2": "noise"}}
(EXP / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")

# ---------------------------------------------------------------- report
hdr("Files (bytes raw / gzip)")
for f in sorted(EXP.iterdir()):
    b = f.read_bytes()
    print(f"  {f.name:20s} {len(b):10,d} / {len(gzip.compress(b)):10,d}")

hdr("Round trip: edges.bin vs parquet")
back = np.fromfile(EXP / "edges.bin", dtype=DT)
assert back.shape[0] == ed.height
for i, c in enumerate(CLS):
    print(f"  {c:18s} file {int((back['cls'] == i).sum()):6d}  parquet {int((ed['cls'] == c).sum()):6d}")
vis = int(((hid[back["src"]] == 0) & (hid[back["tgt"]] == 0)).sum())
print(f"  edges with both endpoints visible by default: {vis} of {len(back)}")

hdr("Hidden by default")
print(f"  chemical (h=1): {int((hid == 1).sum())}   noise (h=2): {int((hid == 2).sum())}")
print(pl.DataFrame(rows).filter(pl.col("h") == 2).select("l", "t", "d").sort("d", descending=True))

hdr("Label cleaning")
print(f"  {len(changed)} labels changed; sample:")
for o, c in changed[:12]:
    print(f"    {o!r} -> {c!r}")
m = pl.DataFrame(rows).filter(pl.col("t") == "metabolite")
dup = m.group_by("l").len().filter(pl.col("len") > 1)
print("  duplicate cleaned metabolite labels:", dup.rows() or "none")

hdr("Top visible metabolites by degree")
print(m.filter(pl.col("h") == 0).select("l", "r", "d").sort("d", descending=True).head(12))
print("\nWrote data/export/")
