from pathlib import Path
import polars as pl

R = Path("data/raw")
pl.Config.set_tbl_rows(70); pl.Config.set_tbl_cols(8); pl.Config.set_fmt_str_lengths(45)
def tsv(p):
    return pl.read_csv(p, separator="\t", infer_schema_length=0, quote_char=None, truncate_ragged_lines=True)
def hdr(t): print(f"\n===== {t} =====")

# SIGNOR side files: raw text, since the separator is unknown
for name in ["SIGNOR_complexes.csv", "SIGNOR_PF.csv", "SIGNOR-PH.csv", "SIGNOR-ST.csv"]:
    hdr(name)
    lines = (R/"signor"/name).read_text(errors="replace").splitlines()
    print(len(lines), "lines")
    for l in lines[:3]: print(l[:300])

# Metabolite-like entities and their edges in SIGNOR
rel = tsv(R/"signor/all_data_30_09_26.tsv")
ents = pl.concat([rel.select(name="ENTITYA", id="IDA", type="TYPEA"),
                  rel.select(name="ENTITYB", id="IDB", type="TYPEB")])
hdr("entity types (unique entities)")
print(ents.unique().group_by("type").len().sort("len", descending=True))
pat = r"(?i)(^atp|^adp|^amp|^gtp|^gdp|^nad|cyclic|acetyl-coa|calcium|hydron|glucose|phosphate|adenosine)"
hits = (ents.filter(pl.col("name").str.contains(pat) & ~pl.col("type").is_in(["protein", "complex", "proteinfamily"]))
            .group_by(["name", "id", "type"]).len().sort("len", descending=True))
hdr("metabolite-like entities (name, id, type, n edges)")
print(hits.head(60))
ids = hits["id"].to_list()
touch = rel.filter(pl.col("IDA").is_in(ids) | pl.col("IDB").is_in(ids))
hdr(f"edges touching those: {touch.height}; by mechanism")
print(touch.group_by("MECHANISM").len().sort("len", descending=True))
hdr("regulation-type edges only (not modification/precursor)")
reg = touch.filter(~pl.col("MECHANISM").is_in(["chemical modification", "precursor of"]))
print(reg.height)
print(reg.select("ENTITYA", "ENTITYB", "EFFECT", "MECHANISM").head(70))

# Human-GEM
for f in ["reactions", "genes"]:
    d = tsv(R/f"humangem/{f}.tsv")
    hdr(f"Human-GEM {f}")
    print(d.shape)
    print({k: (v[:70] if v else v) for k, v in d.row(0, named=True).items()})

# HPA
hpa = tsv(R/"hpa/subcellular_location.tsv")
hdr("HPA")
print(hpa.shape, hpa.columns)
print(hpa["Main location"].drop_nulls().str.split(";").explode().value_counts(sort=True))

# UniProt ID mapping
idm = pl.read_csv(R/"uniprot/HUMAN_9606_idmapping.dat.gz", separator="\t", has_header=False,
                  new_columns=["acc", "type", "value"], quote_char=None, infer_schema_length=0)
hdr("UniProt idmapping")
print(idm.shape)
print(idm["type"].value_counts(sort=True).head(30))
