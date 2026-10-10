"""Step 6b: node_info.json for the viewer, keyed by UniProt accession (proteins in the graph only).
Source: data/raw/uniprot/human_function_disease.tsv (reviewed human entries; UniProt, CC BY 4.0).
"""
import gzip
import json
import os
import re
from pathlib import Path

import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
RAW, PROC, EXP = ROOT / "data/raw", ROOT / "data/processed", ROOT / "data/export"
EXP.mkdir(parents=True, exist_ok=True)


def strip_tags(s):
    return re.sub(r"\{[^}]*\}", "", s or "")


def shorten(s, limit):
    if len(s) <= limit:
        return s
    cut = s[:limit]
    k = cut.rfind(". ")
    if k > 120:
        return cut[:k + 1]
    k = cut.rfind(" ")
    return cut[:k].rstrip(" ,;") + "…"


def clean_func(s):
    s = strip_tags(s)
    s = re.sub(r"^\s*FUNCTION:\s*", "", s)
    s = re.sub(r"\[Isoform[^\]]*\]:.*", "", s, flags=re.S)
    s = s.split("Note=")[0]
    s = re.sub(r"\s*\([^()]*PubMed:[^()]*\)", "", s)
    s = re.sub(r"\s*\(PubMed:[^)]*\)", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\.\s+\.", ".", s)
    return shorten(s, 450)


def split_names(s):
    """'Tuberin (Tuberous sclerosis 2 protein)' -> ('Tuberin', ['Tuberous sclerosis 2 protein'])"""
    main, alts, cur, depth = [], [], [], 0
    for ch in strip_tags(s):
        if ch == "(":
            if depth == 0:
                cur = []
            else:
                cur.append(ch)
            depth += 1
        elif ch == ")" and depth:
            depth -= 1
            if depth == 0:
                alts.append("".join(cur).strip())
            else:
                cur.append(ch)
        elif depth:
            cur.append(ch)
        else:
            main.append(ch)
    name = re.sub(r"\s+", " ", "".join(main)).strip()
    out = []
    for a in alts:
        if a and not a.startswith("EC ") and len(a) <= 60 and a.lower() != name.lower() and a not in out:
            out.append(a)
    return name, out[:3]


def clean_note(s):
    s = re.sub(r"^\s*Note=", "", s)
    return shorten(re.sub(r"\s+", " ", s).strip(" ;"), 220)


def diseases(s):
    named, notes = [], []
    for chunk in re.split(r"(?:^|\s)DISEASE:\s*", strip_tags(s)):
        chunk = chunk.strip()
        if not chunk:
            continue
        if chunk.startswith("Note="):
            notes.append(clean_note(chunk))
            continue
        name, mim = "", ""
        m = re.match(r"(.+?)\s*\[MIM:(\d+)\]", chunk)
        if m and len(m.group(1)) <= 120:
            name, mim = m.group(1), m.group(2)
        else:
            m = re.match(r"(.+?):\s", chunk)
            name = m.group(1) if m and len(m.group(1)) <= 120 else chunk[:80]
        name = re.sub(r"\s+", " ", name).strip()
        if name and all(name != n for n, _ in named):
            named.append((name, mim))
    return named, (notes[0] if notes else "")


nodes = pl.read_parquet(PROC / "nodes.parquet").filter(pl.col("kind") == "protein").select("key", "label")
acc2lab = {k[2:]: l for k, l in nodes.iter_rows()}
u = pl.read_csv(RAW / "uniprot/human_function_disease.tsv", separator="\t",
                infer_schema_length=0, quote_char=None)
u = u.rename(dict(zip(u.columns, ["acc", "gene", "pname", "func", "dis"])))

info = {}
for acc, gene, pname, func, dis in u.iter_rows():
    if acc not in acc2lab:
        continue
    rec = {}
    p, alts = split_names(pname)
    if p:
        rec["p"] = p
    if alts:
        rec["a"] = alts
    f = clean_func(func)
    if f:
        rec["f"] = f
    named, note = diseases(dis)
    if named:
        rec["d"] = [[n, m] for n, m in named[:12]]
        if len(named) > 12:
            rec["dm"] = len(named) - 12
    elif note:
        rec["n"] = note
    if rec:
        info[acc] = rec

blob = json.dumps(info, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
(EXP / "node_info.json").write_bytes(blob)

print(f"proteins in graph: {len(acc2lab)}; with any info: {len(info)}")
print("  with name:", sum("p" in r for r in info.values()), " function:", sum("f" in r for r in info.values()),
      " named disease:", sum("d" in r for r in info.values()), " disease note only:", sum("n" in r for r in info.values()))
print(f"node_info.json: {len(blob):,} bytes raw / {len(gzip.compress(blob)):,} gzip")
print("\n== sanity (all should be 0 except the last) ==")
print("  function text still containing 'PubMed':", sum("PubMed" in r.get("f", "") for r in info.values()))
print("  disease names containing 'Note=' or 'DISEASE':",
      sum(any(("Note=" in n or "DISEASE" in n) for n, _ in r.get("d", [])) for r in info.values()))
print("  names containing '{' :", sum("{" in r.get("p", "") for r in info.values()))
print("  function text cut with an ellipsis (ok if small):", sum(r.get("f", "").endswith("…") for r in info.values()))

lab2acc = {l: a for a, l in acc2lab.items()}
print("\n== samples ==")
for g in ["TSC2", "AKT1", "MTOR", "PRKAA1", "STRADA", "TP53"]:
    a = lab2acc.get(g)
    print(f"\n--- {g} ({a})")
    print(json.dumps(info.get(a), ensure_ascii=False, indent=1)[:1500] if a else "not in graph")
