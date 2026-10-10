"""Shared helpers for the Open Targets build scripts."""
import os
from collections import defaultdict
from pathlib import Path

import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
RAW, PROC, EXP = ROOT / "data/raw", ROOT / "data/processed", ROOT / "data/export"
EXP.mkdir(parents=True, exist_ok=True)
OT = sorted(p for p in (RAW / "opentargets").iterdir() if p.is_dir())[-1]


def sc(name):
    return pl.scan_parquet(str(OT / name / "**" / "*.parquet"))


def ens_maps():
    """-> (ens2acc {ENSG: set of accessions in the graph}, acc2idx {accession: node index}, ens2sym)."""
    keys = pl.read_parquet(PROC / "layout.parquet")["key"].to_list()      # same order as nodes.json
    acc2idx = {k[2:]: i for i, k in enumerate(keys) if k.startswith("P:")}
    ol = list(acc2idx)
    tg = (sc("target").select(pl.col("id").alias("ensg"), "approvedSymbol", "proteinIds")
          .explode("proteinIds")
          .with_columns(src=pl.col("proteinIds").struct.field("source"), acc=pl.col("proteinIds").struct.field("id"))
          .select("ensg", "approvedSymbol", "src", "acc").collect())
    ens2acc = defaultdict(set)
    for e, a in tg.filter((pl.col("src") == "uniprot_swissprot") & pl.col("acc").is_in(ol)).select("ensg", "acc").iter_rows():
        ens2acc[e].add(a)
    have = set(ens2acc)
    for e, a in (tg.filter(pl.col("src").is_in(["uniprot_obsolete", "uniprot_trembl"]) & pl.col("acc").is_in(ol))
                 .select("ensg", "acc").iter_rows()):
        if e not in have:
            ens2acc[e].add(a)
    return ens2acc, acc2idx, dict(zip(tg["ensg"], tg["approvedSymbol"]))
