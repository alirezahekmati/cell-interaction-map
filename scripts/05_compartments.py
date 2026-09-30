"""Step 4: assign every node to one coarse cellular region.

Reads data/processed/{nodes,edges,uniprot_map}.parquet, HPA locations, Human-GEM genes.
Writes data/processed/compartments.parquet:
  key, label, kind, region, region_source, regions_all, deg_all
"""
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
RAW, OUT = ROOT / "data/raw", ROOT / "data/processed"

# tie-break order: specific organelles before the big generic ones
PRIORITY = ["lysosome", "mitochondrion", "peroxisome", "golgi", "er", "vesicle",
            "plasma_membrane", "cytosol", "nucleus", "extracellular"]
# HPA term -> region: first substring match wins, anything else becomes cytosol
RULES = [
    ("extracellular", "extracellular"), ("secreted", "extracellular"),
    ("plasma membrane", "plasma_membrane"), ("cell junction", "plasma_membrane"),
    ("focal adhesion", "plasma_membrane"),
    ("nucle", "nucleus"), ("chromosom", "nucleus"),
    ("mitochond", "mitochondrion"),
    ("endoplasmic", "er"),
    ("golgi", "golgi"),
    ("lysosom", "lysosome"),
    ("peroxisom", "peroxisome"),
    ("endosom", "vesicle"), ("vesicle", "vesicle"),
]
GEM_LETTER = {"c": "cytosol", "e": "extracellular", "r": "er", "m": "mitochondrion",
              "i": "mitochondrion", "x": "peroxisome", "l": "lysosome", "g": "golgi",
              "n": "nucleus"}

pl.Config.set_tbl_rows(70)
pl.Config.set_tbl_cols(14)
pl.Config.set_fmt_str_lengths(40)
pl.Config.set_tbl_width_chars(200)


def hdr(t):
    print(f"\n===== {t} =====")


def read_tsv(p):
    df = pl.read_csv(p, separator="\t", infer_schema_length=0, quote_char=None)
    return df.rename({c: c.strip('"') for c in df.columns})


def to_region(term):
    t = term.lower()
    for k, r in RULES:
        if k in t:
            return r
    return "cytosol"


def terms(s):
    return [x.strip().strip('"') for x in re.split(r"[;|]", s or "") if x.strip()]


def gem_regions(s):
    """Human-GEM compartment field: letters ('ceglmnrx'), or names, separated or not."""
    out = set()
    for tok in re.split(r"[;,|]", s or ""):
        tok = tok.strip().strip('"').lower()
        if not tok:
            continue
        if all(ch in GEM_LETTER for ch in tok):
            out |= {GEM_LETTER[ch] for ch in tok}
        else:
            out.add(to_region(tok))
    return out


# ---------------------------------------------------------------- HPA
hdr("Human Protein Atlas")
hpa = read_tsv(RAW / "hpa/subcellular_location.tsv")
hpa_ens, hpa_name, term_n = {}, {}, Counter()
for r in hpa.iter_rows(named=True):
    main, add = terms(r.get("Main location")), terms(r.get("Additional location"))
    term_n.update(main)
    m = {to_region(t) for t in main} or {to_region(t) for t in add}
    if not m:
        continue
    a = m | {to_region(t) for t in add}
    hpa_ens[r["Gene"]] = (m, a)
    if r.get("Gene name"):
        hpa_name[r["Gene name"].upper()] = (m, a)
print(f"{len(hpa_ens)} genes with a location")
print("Main-location terms -> region   (* = collapsed into cytosol by default)")
for t, n in term_n.most_common():
    r = to_region(t)
    star = "*" if r == "cytosol" and "cytosol" not in t.lower() else " "
    print(f"  {n:5d} {star} {t:38s} {r}")

# ---------------------------------------------------------------- UniProt map
hdr("UniProt map")
um = pl.read_parquet(OUT / "uniprot_map.parquet")
print("columns:", um.columns)
um = um.rename(dict(zip(um.columns[:3], ["acc", "type", "value"])))
acc_ens, acc_name = defaultdict(list), {}
for acc, typ, val in um.select("acc", "type", "value").iter_rows():
    if typ == "Ensembl":
        acc_ens[acc].append(val)
    elif typ == "Gene_Name":
        acc_name.setdefault(acc, val)
print(f"{len(acc_ens)} accessions with Ensembl ids, {len(acc_name)} with gene names")

# ---------------------------------------------------------------- Human-GEM genes
hdr("Human-GEM genes")
gem = read_tsv(RAW / "humangem/genes.tsv")
gem = gem.with_columns(pl.all().str.strip_chars('"'))
print("columns:", gem.columns)
ucol = next((c for c in gem.columns if "uniprot" in c.lower()), None)
ccol = next((c for c in gem.columns if "compartment" in c.lower()), None)
gem_acc = {}
if ucol and ccol:
    for u, c in gem.select(ucol, ccol).iter_rows():
        rs = gem_regions(c)
        if rs:
            for acc in re.split(r"[;,|\s]+", u or ""):
                if acc:
                    gem_acc[acc] = rs
    print(f"{len(gem_acc)} accessions with a compartment; first row:", gem.select(ucol, ccol).row(0))
else:
    print("could not find uniprot/compartment columns; Human-GEM fallback disabled")

# ---------------------------------------------------------------- graph
hdr("Graph")
nodes = pl.read_parquet(OUT / "nodes.parquet").select("key", "label", "kind", "gem_compartments", "deg_all")
edges = pl.read_parquet(OUT / "edges.parquet").select("src", "tgt", "cls")
kind = dict(zip(nodes["key"], nodes["kind"]))
nbr, members = defaultdict(set), defaultdict(set)
for s, t, c in edges.iter_rows():
    if c == "membership":
        members[s].add(t)          # complex/family -> member
    else:
        nbr[s].add(t)
        nbr[t].add(s)
print(f"{nodes.height} nodes, {edges.height} edges")

region, source, cands = {}, {}, {}


def pick(votes, allowed, order):
    allowed = set(allowed) if allowed else set(order)
    v = {r: n for r, n in votes.items() if r in allowed}
    if v:
        top = max(v.values())
        best = {r for r, n in v.items() if n == top}
        if len(best) == 1:
            return next(iter(best)), "vote"
        allowed = best
    return next(r for r in order if r in allowed), "priority"


def votes_from(neigh, snap):
    return Counter(snap[n] for n in neigh if n in snap)


def settle(keys, tag, order=PRIORITY, allowed=None):
    """Neighbour vote in rounds (each round sees the previous one); leftovers -> cytosol."""
    todo = list(keys)
    for _ in range(4):
        snap, new = dict(region), {}
        for k in todo:
            v = votes_from(nbr[k], snap)
            if v:
                new[k] = pick(v, allowed(k) if allowed else None, order)[0]
        for k, r in new.items():
            region[k], source[k] = r, tag
        todo = [k for k in todo if k not in region]
        if not todo:
            break
    for k in todo:
        region[k], source[k] = "cytosol", "default"


# ---- 0b. UniProt subcellular-location text (fallback after HPA and Human-GEM)
hdr("UniProt subcellular location")
UP_RULES = [("cell membrane", "plasma_membrane"), ("cell surface", "plasma_membrane"),
            ("perinuclear", "cytosol")]


def up_region(phrase):
    t = phrase.lower()
    for k, r in UP_RULES + RULES:
        if k in t:
            return r
    if re.search(r"cytoplasm|cytosol|cytoskeleton|cell projection", t):
        return "cytosol"
    return None          # 'Membrane', 'Endomembrane system', ...: no usable signal


def up_regions(text):
    """Ordered, de-duplicated regions from UniProt 'Subcellular location [CC]' text."""
    s = re.sub(r"\{[^}]*\}", "", text or "")
    s = re.sub(r"^\s*SUBCELLULAR LOCATION:\s*", "", s)
    chunks = [c for c in re.split(r"\[Isoform[^\]]*\]:", s) if c.strip()]
    if not chunks:
        return []
    s = chunks[0].split("Note=")[0]
    out = []
    for entry in re.split(r"\.(?:\s+|$)", s):
        r = up_region(entry.split(";")[0])
        if r and r not in out:
            out.append(r)
    return out


up_loc = {}
_upf = RAW / "uniprot/human_subcellular.tsv"
if _upf.exists():
    _u = pl.read_csv(_upf, separator="\t", infer_schema_length=0, quote_char=None)
    for _acc, _txt in zip(_u[_u.columns[0]], _u[_u.columns[2]]):
        _rs = up_regions(_txt)
        if _rs:
            up_loc[_acc] = _rs
    print(f"{len(up_loc)} reviewed entries with a mappable location")
else:
    print("missing data/raw/uniprot/human_subcellular.tsv; UniProt fallback disabled")

# ---- 1. proteins
hdr("Proteins")
multi, none_, n_prot = [], [], 0
for k in nodes.filter(pl.col("kind") == "protein")["key"]:
    n_prot += 1
    acc = k[2:]
    m = a = None
    for g in acc_ens.get(acc, []):
        if g in hpa_ens:
            gm, ga = hpa_ens[g]
            m, a = (gm, ga) if m is None else (m | gm, a | ga)
    nm = acc_name.get(acc, "").upper()
    if m is None and nm in hpa_name:
        m, a = hpa_name[nm]
    src = "hpa"
    if m is None and acc in gem_acc:
        m = a = gem_acc[acc]
        src = "gem"
    if m is None and acc in up_loc:
        m, a = {up_loc[acc][0]}, set(up_loc[acc])
        src = "uniprot"
    if m is None:
        none_.append(k)
        continue
    cands[k] = a
    if len(m) == 1:
        region[k], source[k] = next(iter(m)), src
    else:
        multi.append((k, m, src))
snap = dict(region)
for k, m, src in multi:
    region[k], how = pick(votes_from(nbr[k], snap), m, PRIORITY)
    source[k] = f"{src}_{how}"
OV_LYSO = {"RRAGA", "RRAGB", "RRAGC", "RRAGD",
           "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
           "SLC38A9", "RHEB", "FLCN", "FNIP1", "FNIP2",
           "DEPDC5", "NPRL2", "NPRL3",                      # GATOR1
           "KPTN", "ITFG2", "KICS2", "C12ORF66", "SZT2"}    # KICSTOR
lab = dict(zip(nodes["key"], nodes["label"]))
hit = set()
for k in nodes.filter(pl.col("kind") == "protein")["key"]:
    L = lab[k].upper()
    if L in OV_LYSO or L.startswith("ATP6V"):               # v-ATPase subunits
        region[k], source[k] = "lysosome", "override"
        hit.add(L)
print(f"lysosome overrides: {len(hit)} placed; names not found: {sorted(OV_LYSO - hit)}")
# HPA lists nucleoplasm as the main location for several cytosolic signalling kinases
OV_CYTO = {"PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "PRKAG2", "PRKAG3",
           "STK11", "STRADA", "STRADB", "CAB39", "AKT1", "AKT2", "AKT3"}
moved = []
for k in nodes.filter(pl.col("kind") == "protein")["key"]:
    if lab[k].upper() in OV_CYTO:
        if region.get(k) != "cytosol":
            moved.append(f"{lab[k]}:{region.get(k, 'unplaced')}")
        region[k], source[k] = "cytosol", "override"
print(f"cytosol overrides: moved {moved}")
# mTORC1 acts at the lysosome surface (RHEB, RAG, RPTOR are already there); regions_all keeps the rest
OV_MTOR = {"MTOR", "RPTOR"}
for k in nodes.filter(pl.col("kind") == "protein")["key"]:
    if lab[k].upper() in OV_MTOR:
        region[k], source[k] = "lysosome", "override"
# proteins with no usable UniProt/HPA text (from memory; edit freely)
OV_FIX = {"JAK1": "cytosol", "HBA1": "cytosol", "HBB": "cytosol",
          "ITGA2": "plasma_membrane", "ACVR1": "plasma_membrane", "EGF": "extracellular"}
for k in nodes.filter(pl.col("kind") == "protein")["key"]:
    r = OV_FIX.get(lab[k].upper())
    if r:
        region[k], source[k] = r, "override"
settle([k for k in none_ if k not in region], "nbr")
print(f"{n_prot} proteins: {n_prot - len(multi) - len(none_)} single-region, "
      f"{len(multi)} multi-region, {len(none_)} with no location data")

# ---- 2. complexes and families
todo = [k for k in nodes["key"] if kind[k] in ("complex", "family")]
n_cx = len(todo)
for _ in range(6):
    snap, new = dict(region), {}
    for k in todo:
        v = votes_from(members.get(k, ()), snap)
        if v:
            new[k] = pick(v, None, PRIORITY)[0]
    for k, r in new.items():
        region[k], source[k] = r, "members"
    todo = [k for k in todo if k not in region]
    if not todo:
        break
print(f"\n{n_cx} complexes/families: {n_cx - len(todo)} placed from members, {len(todo)} left for neighbours")
settle(todo, "nbr")
for k in nodes.filter(pl.col("kind") == "complex")["key"]:
    if lab[k] == "mTORC1":
        region[k], source[k] = "lysosome", "override"
    elif lab[k] == "mTORC2":
        region[k], source[k] = "cytosol", "override"

# ---- 3. metabolites (vote restricted to the compartments Human-GEM knows them in)
allowed = {}
for k, g in nodes.filter(pl.col("kind") == "metabolite").select("key", "gem_compartments").iter_rows():
    allowed[k] = gem_regions(g) or None
    if allowed[k]:
        cands[k] = allowed[k]
snap = dict(region)
for k, s in allowed.items():
    if s == {"extracellular"}:
        region[k], source[k] = "extracellular", "gem"
        continue
    v = votes_from(nbr[k], snap)
    if s:
        v = Counter({r: c for r, c in v.items() if r in s})
    top = v.most_common(1)
    if top and top[0][0] != "cytosol" and top[0][1] >= 3 and top[0][1] >= 0.6 * sum(v.values()):
        region[k], source[k] = top[0][0], "nbr_strong"
    elif s is None or "cytosol" in s:
        region[k], source[k] = "cytosol", "gem_default"
    else:
        region[k], source[k] = pick(v, s, PRIORITY)[0], "gem_only"

# ---- 4. everything else, then rules that must not act as votes
settle([k for k in nodes["key"] if k not in region and kind[k] not in ("stimulus", "mirna")], "nbr")
for k in nodes["key"]:
    if kind[k] == "stimulus":
        region[k], source[k] = "extracellular", "rule"
    elif kind[k] == "mirna":
        region[k], source[k] = "cytosol", "rule"

out = pl.DataFrame({
    "key": nodes["key"], "label": nodes["label"], "kind": nodes["kind"],
    "region": [region[k] for k in nodes["key"]],
    "region_source": [source[k] for k in nodes["key"]],
    "regions_all": [",".join(sorted(cands.get(k, ()))) for k in nodes["key"]],
    "deg_all": nodes["deg_all"],
})
out.write_parquet(OUT / "compartments.parquet")

# ---------------------------------------------------------------- report
hdr("Nodes per kind and region")
print(out.group_by("kind", "region").len().pivot(on="region", index="kind", values="len").fill_null(0))

hdr("How each node got its region")
print(out.group_by("kind", "region_source").len().sort("kind", "len", descending=[False, True]))

hdr("Spot check: mTOR figure players")
want = ["TSC2", "TSC1", "TSC", "STK11", "PRKAA1", "PRKAG1", "AMPK", "MTOR", "mTORC1", "mTORC2",
        "RPTOR", "RHEB", "AKT1", "PDPK1", "RRAGA", "LAMTOR1", "SESN2", "SLC38A9", "FLCN",
        "DEPDC5", "NPRL2", "EGFR", "INSR", "ATP", "ADP", "AMP",
        "TNF", "LRP6", "FN1", "ADRA2A", "PPP2CA", "CXCR4", "MTNR1A"]
print(out.filter(pl.col("label").is_in(want) & pl.col("kind").is_in(["protein", "complex", "family", "metabolite"]))
      .select("label", "kind", "region", "region_source", "regions_all").sort("label"))
print("labels not found:", sorted(set(want) - set(out["label"])))

hdr("Well-connected proteins placed without location data")
print(out.filter((pl.col("kind") == "protein") & pl.col("region_source").is_in(["nbr", "default"]))
      .sort("deg_all", descending=True).head(15).select("label", "deg_all", "region", "region_source"))

print("\nWrote data/processed/compartments.parquet")
