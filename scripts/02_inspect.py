from pathlib import Path
import polars as pl

R = Path("data/raw")
pl.Config.set_tbl_rows(40); pl.Config.set_tbl_cols(12); pl.Config.set_fmt_str_lengths(40)

def tsv(p, **kw):
    return pl.read_csv(p, separator="\t", infer_schema_length=0, quote_char=None,
                       truncate_ragged_lines=True, **kw)

def hdr(t): print(f"\n{'='*12} {t} {'='*12}")

# ---- SIGNOR ----
hdr("SIGNOR files")
files = sorted(p for p in (R/"signor").iterdir() if p.is_file())
for p in files: print(p.name, f"{p.stat().st_size/1e6:.1f} MB")

rel = None
for p in files:
    try:
        d = tsv(p)
        if "ENTITYA" in d.columns:
            rel = d; print("relations file:", p.name); break
    except Exception as e:
        print("skip", p.name, e)

if rel is not None:
    print(rel.shape)
    for c in ["EFFECT", "MECHANISM", "TYPEA", "TYPEB", "DATABASEA", "TAX_ID"]:
        hdr(f"SIGNOR {c}")
        print(rel[c].value_counts(sort=True).head(25))
    chem = rel.filter((pl.col("TYPEA") == "chemical") | (pl.col("TYPEB") == "chemical"))
    hdr(f"SIGNOR edges touching chemicals: {chem.height}")
    names = pl.concat([
        chem.filter(pl.col("TYPEA") == "chemical").select(n="ENTITYA", i="IDA"),
        chem.filter(pl.col("TYPEB") == "chemical").select(n="ENTITYB", i="IDB")])
    print(names.group_by(["n", "i"]).len().sort("len", descending=True).head(25))
    hdr("SIGNOR edges with ATP / AMP / ADP (by name)")
    m = rel.filter(pl.col("ENTITYA").str.to_uppercase().is_in(["ATP", "AMP", "ADP"]) |
                   pl.col("ENTITYB").str.to_uppercase().is_in(["ATP", "AMP", "ADP"]))
    print(m.select("ENTITYA", "IDA", "ENTITYB", "IDB", "EFFECT", "MECHANISM"))
    hdr("SIGNOR TSC2 in/out")
    t = rel.filter((pl.col("ENTITYA") == "TSC2") | (pl.col("ENTITYB") == "TSC2"))
    print("edges:", t.height)
    print(t.select("ENTITYA", "ENTITYB", "EFFECT", "MECHANISM").head(15))
    cx = rel.filter(pl.col("IDA").str.starts_with("SIGNOR-") | pl.col("IDB").str.starts_with("SIGNOR-"))
    hdr(f"SIGNOR edges touching SIGNOR-* entities (complexes/families/etc.): {cx.height}")

# ---- OmniPath ----
hdr("OmniPath")
op = pl.read_parquet(R/"omnipath/interactions.parquet")
print(op.shape); print(op.schema)
print(op.head(2))
for c in ["is_stimulation", "is_inhibition", "is_directed"]:
    if c in op.columns: print(c, op[c].value_counts().to_dicts())

# ---- Human-GEM ----
hdr("Human-GEM metabolites")
mets = tsv(R/"humangem/metabolites.tsv"); print(mets.shape); print(mets.columns)
print(mets.head(2))
print("compartment letter counts:", mets["mets"].str.slice(-1).value_counts(sort=True).to_dicts())
atp = mets.filter(pl.any_horizontal([pl.col(c).str.to_uppercase() == "ATP" for c in mets.columns]))
print("ATP rows:"); print(atp)
hdr("Human-GEM reactions"); rx = tsv(R/"humangem/reactions.tsv"); print(rx.shape); print(rx.columns); print(rx.head(2))
hdr("Human-GEM genes"); g = tsv(R/"humangem/genes.tsv"); print(g.shape); print(g.columns); print(g.head(2))

# ---- HPA ----
hdr("HPA subcellular")
hpa = tsv(R/"hpa/subcellular_location.tsv"); print(hpa.shape)
loc = hpa["Main location"].drop_nulls().str.split(";").explode()
print(loc.value_counts(sort=True).head(40))

# ---- UniProt ID mapping ----
hdr("UniProt idmapping")
idm = pl.read_csv(R/"uniprot/HUMAN_9606_idmapping.dat.gz", separator="\t", has_header=False,
                  new_columns=["acc", "type", "value"], quote_char=None, infer_schema_length=0)
print(idm.shape)
print(idm["type"].value_counts(sort=True).head(30))
