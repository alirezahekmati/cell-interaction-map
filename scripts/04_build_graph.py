"""Step 3: normalize SIGNOR + OmniPath into one node table and one merged edge table.

Writes to data/processed/:
  nodes.parquet        one row per entity (protein, complex, family, metabolite, chemical, ...)
  evidence.parquet     one row per supporting statement (kept so every edge can be audited)
  edges.parquet        one row per (source, target, edge class); stimulation/inhibition merged
  uniprot_map.parquet  small slice of the UniProt ID mapping (Gene_Name, Ensembl) for later steps

Edge classes:
  regulation         signed/unsigned causal edges (default layer)
  substrate_product  enzyme makes/consumes a metabolite ("chemical modification", "precursor of")
  transport          relocalization of a metabolite (channels, transporters)
  binding            "form complex" (undirected)
  membership         complex/family -> member protein
"""
import os
import re
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

warnings.filterwarnings("ignore", category=DeprecationWarning)

ROOT = Path(os.environ.get("BIO_ROOT", "."))
RAW, OUT = ROOT / "data/raw", ROOT / "data/processed"
OUT.mkdir(parents=True, exist_ok=True)

ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})")
CHEM_ALIAS = {
    "α-d-glucose": "glucose", "alpha-d-glucose": "glucose", "β-d-glucose": "glucose",
    "d-glucose": "glucose", "nad+": "nad", "nadp+": "nadp",
    "camp": "3',5'-cyclic amp", "cgmp": "3',5'-cyclic gmp",
}
CHEM_ALIAS.update({
    "adenosine 5'-monophosphate": "amp",   # SIGNOR lists AMP under two entities
    "d-glucopyranose": "glucose",
    "udp-d-glucose": "udp-alpha-d-glucose",
    "glutamic acid": "l-glutamate",
    "aspartic acid": "l-aspartate",
    "glutamine": "l-glutamine",
    "udp-d-galactose": "udp-alpha-d-galactose",
})
# metal ions keep their charge: iron(2+) and iron(3+) are different things
METAL_IONS = {"iron", "copper", "zinc", "calcium", "magnesium", "manganese", "cobalt",
              "nickel", "sodium", "potassium", "chromium", "mercury", "lithium"}
CHEM_KINDS = ["metabolite", "chemical"]
SUBSTRATE_MECH = ["chemical modification", "precursor of"]
KEYS = ["src", "tgt", "cls"]
COLS = KEYS + ["sign", "mechanism", "effect", "resource", "pmids", "score", "ev_id"]

pl.Config.set_tbl_rows(60)
pl.Config.set_tbl_cols(12)
pl.Config.set_fmt_str_lengths(40)


def hdr(t):
    print(f"\n===== {t} =====")


def read_tsv(p):
    return pl.read_csv(p, separator="\t", infer_schema_length=0, quote_char=None,
                       truncate_ragged_lines=True)


def read_side(p):
    df = pl.read_csv(p, separator=";", infer_schema_length=0, encoding="utf8-lossy",
                     truncate_ragged_lines=True)
    return df.select(df.columns[:3])


def norm_chem(name):
    """ATP / ATP(4-) -> 'atp'; 'L-glutamine zwitterion' -> 'l-glutamine'."""
    n = name.strip().lower()
    n = re.sub(r"\s+zwitterion$", "", n)
    m = re.match(r"^(.*?)\(\d*[+-]\)$", n)
    if m and m.group(1).strip() not in METAL_IONS:
        n = m.group(1)
    n = re.sub(r"\s+", " ", n).strip()
    return CHEM_ALIAS.get(n, n)


def entity_key(name, typ, ident):
    """-> (node key, kind, uniprot accession, chebi id)"""
    if typ == "protein":
        m = ACC.match(ident)
        acc = m.group(1) if m else ident
        return f"P:{acc}", "protein", acc, ""
    if typ in ("chemical", "smallmolecule"):
        n = norm_chem(name) if name else ident.lower()
        return f"C:{n}", "chemical", "", (ident if ident.startswith("CHEBI:") else "")
    if typ == "complex":
        return ident, "complex", "", ""
    if typ == "proteinfamily":
        return ident, "family", "", ""
    if typ in ("phenotype", "stimulus"):
        return ident, typ, "", ""
    return f"X:{ident or name}", (typ.replace(" ", "_") or "other"), "", ""


# ---------------------------------------------------------------- Human-GEM (metabolite compartments)
hdr("Human-GEM metabolites")
gem = read_tsv(RAW / "humangem/metabolites.tsv")
gem_comp = defaultdict(set)
for mid, chebi in gem.select("mets", "metChEBIID").iter_rows():
    for c in (chebi or "").split(";"):
        if c.strip():
            gem_comp[c.strip()].add(mid[-1])
print(f"{len(gem_comp)} distinct ChEBI ids in Human-GEM")

# ---------------------------------------------------------------- UniProt slice
hdr("UniProt ID mapping")
idm = pl.read_csv(RAW / "uniprot/HUMAN_9606_idmapping.dat.gz", separator="\t", has_header=False,
                  new_columns=["acc", "type", "value"], quote_char=None, infer_schema_length=0)
idm = idm.filter(pl.col("type").is_in(["Gene_Name", "Ensembl"]))
idm.write_parquet(OUT / "uniprot_map.parquet")
gene_name = dict(idm.filter(pl.col("type") == "Gene_Name").unique(subset="acc", keep="first")
                 .select("acc", "value").iter_rows())
print(f"{len(gene_name)} accessions with a gene symbol")

# ---------------------------------------------------------------- SIGNOR entities
hdr("SIGNOR")
sig_file = sorted((RAW / "signor").glob("all_data*.tsv"))[-1]
sig = read_tsv(sig_file)
n0 = sig.height
sig = sig.filter(pl.col("TAX_ID").is_null() | pl.col("TAX_ID").is_in(["9606", "-1"]))
print(f"{sig_file.name}: {n0} rows -> {sig.height} after keeping human / unspecified species")
text_cols = ["ENTITYA", "TYPEA", "IDA", "ENTITYB", "TYPEB", "IDB", "EFFECT", "MECHANISM"]
sig = sig.with_columns([pl.col(c).fill_null("") for c in text_cols])

ents = pl.concat([
    sig.select(name="ENTITYA", type="TYPEA", id="IDA"),
    sig.select(name="ENTITYB", type="TYPEB", id="IDB"),
]).unique()

rows = []
for name, typ, ident in ents.iter_rows():
    key, kind, uni, chebi = entity_key(name, typ, ident)
    rows.append(dict(name=name, type=typ, id=ident, key=key, kind=kind, uniprot=uni,
                     chebi=chebi, is_sm=(typ == "smallmolecule")))
ent_map = pl.DataFrame(rows, infer_schema_length=None)

nodes = {}
for r in rows:
    n = nodes.setdefault(r["key"], dict(key=r["key"], kind=r["kind"], names=Counter(),
                                        uniprot=r["uniprot"], chebi=set(), is_sm=False, gem=[]))
    n["names"][r["name"]] += 1
    if r["chebi"]:
        n["chebi"].add(r["chebi"])
    n["is_sm"] = n["is_sm"] or r["is_sm"]

for n in nodes.values():
    names = [x for x in n["names"] if x]
    if n["kind"] == "chemical":
        n["label"] = min(names, key=lambda s: (len(s), s)) if names else n["key"]
        in_gem = any(c in gem_comp for c in n["chebi"])
        n["kind"] = "metabolite" if (n["is_sm"] or in_gem) else "chemical"
        n["gem"] = sorted(set().union(*[gem_comp.get(c, set()) for c in n["chebi"]]))
    else:
        n["label"] = max(names, key=lambda s: n["names"][s]) if names else n["key"]

kind_of = {k: n["kind"] for k, n in nodes.items()}
chebi_to_key = {c: k for k, n in nodes.items() for c in n["chebi"]}

# ---------------------------------------------------------------- SIGNOR evidence
mapA = ent_map.select(TYPEA="type", IDA="id", ENTITYA="name", src="key")
mapB = ent_map.select(TYPEB="type", IDB="id", ENTITYB="name", tgt="key")
sig_ev = (sig.join(mapA, on=["TYPEA", "IDA", "ENTITYA"], how="left")
             .join(mapB, on=["TYPEB", "IDB", "ENTITYB"], how="left"))
sig_ev = sig_ev.select(
    "src", "tgt",
    sign=(pl.when(pl.col("EFFECT").str.starts_with("up-regulates")).then(1)
            .when(pl.col("EFFECT").str.starts_with("down-regulates")).then(-1)
            .otherwise(0).cast(pl.Int8)),
    mechanism=pl.col("MECHANISM"),
    effect=pl.col("EFFECT"),
    resource=pl.lit("SIGNOR"),
    pmids=pl.col("PMID"),
    score=pl.col("SCORE").cast(pl.Float64, strict=False),
    ev_id=pl.col("SIGNOR_ID"),
)
kinds = pl.DataFrame({"key": list(kind_of), "kind": list(kind_of.values())})
chem = pl.col("src_kind").is_in(CHEM_KINDS) | pl.col("tgt_kind").is_in(CHEM_KINDS)
sig_ev = (sig_ev.join(kinds.rename({"key": "src", "kind": "src_kind"}), on="src", how="left")
                .join(kinds.rename({"key": "tgt", "kind": "tgt_kind"}), on="tgt", how="left")
                .with_columns(cls=(pl.when(pl.col("effect") == "form complex").then(pl.lit("binding"))
                                     .when(chem & pl.col("mechanism").is_in(SUBSTRATE_MECH))
                                     .then(pl.lit("substrate_product"))
                                     .when(chem & (pl.col("mechanism") == "relocalization"))
                                     .then(pl.lit("transport"))
                                     .otherwise(pl.lit("regulation")))))
swap = (pl.col("cls") == "binding") & (pl.col("src") > pl.col("tgt"))   # binding is undirected
sig_ev = sig_ev.with_columns(src=pl.when(swap).then(pl.col("tgt")).otherwise(pl.col("src")),
                             tgt=pl.when(swap).then(pl.col("src")).otherwise(pl.col("tgt")))
sig_ev = sig_ev.filter(pl.col("src").is_not_null() & pl.col("tgt").is_not_null()).select(COLS)
print("SIGNOR evidence rows by class:")
print(sig_ev["cls"].value_counts(sort=True))

# ---------------------------------------------------------------- OmniPath evidence
hdr("OmniPath")
op = pl.read_parquet(RAW / "omnipath/interactions.parquet")
n_op = op.height
op = op.filter(~pl.col("source").str.starts_with("COMPLEX") & ~pl.col("target").str.starts_with("COMPLEX"))
sym = dict(zip(op["source"], op["source_genesymbol"]))
sym.update(zip(op["target"], op["target_genesymbol"]))
op = op.with_columns(
    resource=(pl.col("sources").str.split(";")
              .list.eval(pl.element().filter(~pl.element().str.starts_with("SIGNOR")))
              .list.join(";").fill_null("")))
n_signor_only = op.filter(pl.col("resource") == "").height
op = op.filter(pl.col("resource") != "")
print(f"{n_op} rows; dropped {n_signor_only} that are SIGNOR-only (already loaded directly); {op.height} kept")

base = op.select(src=pl.lit("P:") + pl.col("source"), tgt=pl.lit("P:") + pl.col("target"),
                 resource=pl.col("resource"), pmids=pl.col("references_stripped"),
                 stim=pl.col("is_stimulation"), inhib=pl.col("is_inhibition"),
                 c_stim=pl.col("consensus_stimulation").fill_null(False),
                 c_inhib=pl.col("consensus_inhibition").fill_null(False))
# A record flagged both stimulation and inhibition: when the resources' consensus points at
# exactly one sign, keep only that sign. Keep both when consensus is absent or truly split.
_both = pl.col("stim") & pl.col("inhib")
KEEP_STIM = pl.col("stim") & ~(_both & pl.col("c_inhib") & ~pl.col("c_stim"))
KEEP_INHIB = pl.col("inhib") & ~(_both & pl.col("c_stim") & ~pl.col("c_inhib"))
n_both = base.filter(_both).height
n_resolved = base.filter(_both & (pl.col("c_stim") != pl.col("c_inhib"))).height


def op_frame(df, sign, effect):
    return df.select(
        "src", "tgt", "resource", "pmids",
        cls=pl.lit("regulation"), sign=pl.lit(sign, dtype=pl.Int8), mechanism=pl.lit(""),
        effect=pl.lit(effect), score=pl.lit(None, dtype=pl.Float64), ev_id=pl.lit("OP"),
    ).select(COLS)


op_ev = pl.concat([
    op_frame(base.filter(KEEP_STIM), 1, "stimulation"),
    op_frame(base.filter(KEEP_INHIB), -1, "inhibition"),
    op_frame(base.filter(~pl.col("stim") & ~pl.col("inhib")), 0, ""),
])
print(f"OmniPath evidence rows: {op_ev.height}; {n_both} records flagged both signs, "
      f"{n_resolved} of them resolved to one sign by consensus")

# ---------------------------------------------------------------- membership (complexes / families)
hdr("Complex and family membership")


def member_key(m):
    if re.match(r"^SIGNOR-(C|PF)\d+$", m):
        return m
    if m.startswith("CHEBI:"):
        return chebi_to_key.get(m)
    mm = ACC.match(m)
    return f"P:{mm.group(1)}" if mm else None


parents, mem, unmatched = {}, set(), Counter()
for fname, kind in [("SIGNOR_complexes.csv", "complex"), ("SIGNOR_PF.csv", "family")]:
    for sid, name, members in read_side(RAW / "signor" / fname).iter_rows():
        parents[sid] = (kind, name or sid)
        for m in (members or "").split(","):
            m = m.strip()
            if not m:
                continue
            k = member_key(m)
            if k is None:
                unmatched[re.sub(r"[0-9]+.*$", "", m) or m[:8]] += 1
            else:
                mem.add((sid, k))
print(f"{len(parents)} complexes/families, {len(mem)} membership links; unmatched member ids by prefix: "
      f"{dict(unmatched.most_common(8))}")
mem_ev = pl.DataFrame({"src": [a for a, _ in mem], "tgt": [b for _, b in mem]}).select(
    "src", "tgt", cls=pl.lit("membership"), sign=pl.lit(0, dtype=pl.Int8),
    mechanism=pl.lit("membership"), effect=pl.lit(""), resource=pl.lit("SIGNOR-complexes"),
    pmids=pl.lit(None, dtype=pl.Utf8), score=pl.lit(None, dtype=pl.Float64),
    ev_id=pl.col("src")).select(COLS)

# ---------------------------------------------------------------- nodes that only appear via OmniPath / membership
ev = pl.concat([sig_ev, op_ev, mem_ev]).filter(pl.col("src") != pl.col("tgt"))
for k in (set(ev["src"]) | set(ev["tgt"])) - nodes.keys():
    if k.startswith("P:"):
        acc = k[2:]
        nodes[k] = dict(key=k, kind="protein", label=sym.get(acc) or gene_name.get(acc) or acc,
                        uniprot=acc, chebi=set(), names=Counter(), gem=[])
    elif k in parents:
        nodes[k] = dict(key=k, kind=parents[k][0], label=parents[k][1], uniprot="",
                        chebi=set(), names=Counter(), gem=[])
    else:
        nodes[k] = dict(key=k, kind="other", label=k, uniprot="", chebi=set(), names=Counter(), gem=[])

# ---------------------------------------------------------------- merge to one line per (src, tgt, class)
hdr("Merging")


def uniq_join(col, out, split=None, head=None, count=None):
    s = pl.col(col).str.split(split) if split else pl.col(col)
    d = ev.select(*KEYS, v=s)
    if split:
        d = d.explode("v")
    d = d.filter(pl.col("v").is_not_null() & (pl.col("v") != "")).unique()
    v = pl.col("v").sort()
    if head:
        v = v.head(head)
    agg = [v.str.join(";").alias(out)]
    if count:
        agg.append(pl.len().alias(count))
    return d.group_by(KEYS).agg(agg)


edges = (ev.group_by(KEYS).agg(
            n_stim=(pl.col("sign") == 1).sum(),
            n_inhib=(pl.col("sign") == -1).sum(),
            n_unsigned=(pl.col("sign") == 0).sum(),
            n_evidence=pl.len(),
            max_score=pl.col("score").max())
         .join(uniq_join("resource", "resources", split=";", count="n_resources"), on=KEYS, how="left")
         .join(uniq_join("pmids", "pmids", split=";", head=40, count="n_pmids"), on=KEYS, how="left")
         .join(uniq_join("mechanism", "mechanisms"), on=KEYS, how="left")
         .with_columns(pl.col("n_pmids").fill_null(0), pl.col("n_resources").fill_null(0))
         .with_columns(
             net=(pl.when((pl.col("n_stim") > 0) & (pl.col("n_inhib") > 0)).then(pl.lit("mixed"))
                    .when(pl.col("n_stim") > 0).then(pl.lit("stim"))
                    .when(pl.col("n_inhib") > 0).then(pl.lit("inhib"))
                    .otherwise(pl.lit("unsigned"))),
             directed=pl.col("cls") != "binding"))
print(f"{ev.height} evidence rows -> {edges.height} merged edges")

# ---------------------------------------------------------------- node table
node_df = pl.DataFrame([dict(key=n["key"], label=n["label"], kind=n["kind"], uniprot=n["uniprot"],
                             chebi=";".join(sorted(n["chebi"])), gem_compartments="".join(n["gem"]))
                        for n in nodes.values()], infer_schema_length=None)
reg = edges.filter(pl.col("cls") == "regulation")
node_df = (node_df
           .join(reg.group_by("src").agg(deg_out_reg=pl.len()).rename({"src": "key"}), on="key", how="left")
           .join(reg.group_by("tgt").agg(deg_in_reg=pl.len()).rename({"tgt": "key"}), on="key", how="left")
           .join(pl.concat([edges.select(key="src"), edges.select(key="tgt")])
                 .group_by("key").agg(deg_all=pl.len()), on="key", how="left")
           .join(edges.filter(pl.col("cls") == "membership").group_by("src")
                 .agg(n_members=pl.len()).rename({"src": "key"}), on="key", how="left")
           .with_columns(pl.col(["deg_out_reg", "deg_in_reg", "deg_all", "n_members"]).fill_null(0)))

node_df.write_parquet(OUT / "nodes.parquet")
ev.write_parquet(OUT / "evidence.parquet")
edges.write_parquet(OUT / "edges.parquet")

# ================================================================= report
hdr("Nodes by kind")
print(node_df["kind"].value_counts(sort=True))
hdr("Edges by class and net sign")
print(edges.group_by("cls", "net").len().sort("cls", "net"))

hdr("Metabolite nodes that merged several SIGNOR entries (check these look right)")
multi = sorted([n for n in nodes.values() if n["kind"] == "metabolite" and len(n["names"]) > 1],
               key=lambda n: -len(n["names"]))[:25]
for n in multi:
    print(f"{n['label']:<28} names={sorted(n['names'])}  chebi={sorted(n['chebi'])}  gem={''.join(n['gem'])}")

lab = node_df.select("key", "label")
edges_l = (edges.join(lab.rename({"key": "src", "label": "src_label"}), on="src")
                .join(lab.rename({"key": "tgt", "label": "tgt_label"}), on="tgt"))


def find_keys(terms):
    t = [x.lower() for x in terms]
    return node_df.filter(pl.col("label").str.to_lowercase().is_in(t)
                          | pl.col("key").is_in([f"C:{x}" for x in t]))["key"].to_list()


hdr("Hubs: regulation edges in and out")
for term in ["ATP", "AMP", "ADP", "AMPK", "TSC2", "TSC"]:
    ks = find_keys([term])
    if not ks:
        print(f"\n[{term}] not found")
        continue
    sub = edges_l.filter(pl.col("src").is_in(ks) | pl.col("tgt").is_in(ks))
    by_cls = dict(sub.group_by("cls").len().iter_rows())
    print(f"\n[{term}] nodes={ks} edges by class={by_cls}")
    print(sub.filter(pl.col("cls") == "regulation").sort("n_evidence", descending=True)
             .select("src_label", "tgt_label", "net", "n_stim", "n_inhib", "n_pmids", "resources").head(25))

hdr("Check against the mTOR figure (expected sign from the figure / textbook biology)")
CHECKS = [
    (["AKT", "AKT1", "AKT2", "AKT3"], ["TSC2", "TSC"], "inhib"),
    (["MAPK1", "MAPK3", "ERK1/2"], ["TSC2", "TSC"], "inhib"),
    (["RPS6KA1", "RPS6KA3", "RSK"], ["TSC2", "TSC"], "inhib"),
    (["IKBKB"], ["TSC1", "TSC2", "TSC"], "inhib"),
    (["AMPK", "PRKAA1", "PRKAA2"], ["TSC2", "TSC"], "stim"),
    (["STK11"], ["AMPK", "PRKAA1", "PRKAA2"], "stim"),
    (["DDIT4"], ["TSC2", "TSC"], "stim"),
    (["PDPK1"], ["AKT", "AKT1", "AKT2", "AKT3"], "stim"),
    (["mTORC2"], ["AKT", "AKT1", "AKT2", "AKT3"], "stim"),
    (["TSC2", "TSC"], ["RHEB"], "inhib"),
    (["RHEB"], ["mTORC1", "MTOR"], "stim"),
    (["AMPK", "PRKAA1", "PRKAA2"], ["RPTOR", "mTORC1", "MTOR"], "inhib"),
    (["AMP"], ["AMPK", "PRKAG1", "PRKAA1"], "stim"),
    (["ADP"], ["AMPK", "PRKAG1", "PRKAA1"], "stim"),
    (["ATP"], ["AMPK", "PRKAG1", "PRKAA1"], "inhib"),
]
for srcs, tgts, exp in CHECKS:
    S, T = find_keys(srcs), find_keys(tgts)
    sub = edges_l.filter(pl.col("src").is_in(S) & pl.col("tgt").is_in(T) & (pl.col("cls") == "regulation"))
    tag = f"{'/'.join(srcs[:2])} -> {'/'.join(tgts[:2])}"
    if sub.height == 0:
        print(f"MISS      {tag:<34} expected {exp}")
        continue
    st, inh = sub["n_stim"].sum(), sub["n_inhib"].sum()
    got = "mixed" if st and inh else "stim" if st else "inhib" if inh else "unsigned"
    status = "OK" if got == exp else "MIXED" if got == "mixed" else "OPPOSITE" if got in ("stim", "inhib") else "UNSIGNED"
    pairs = ", ".join(sorted({f"{a}>{b}" for a, b in zip(sub["src_label"], sub["tgt_label"])})[:4])
    print(f"{status:<9} {tag:<34} expected {exp:<6} found {got:<8} stim={st} inhib={inh}  [{pairs}]")

print(f"\nWrote nodes/edges/evidence parquet to {OUT}")
