"""Step 5: precomputed 3D layout (no browser physics).

Reads  data/processed/{compartments,edges}.parquet
Writes data/processed/layout.parquet (key, label, kind, region, x, y, z)
       data/processed/regions.json   (geometry of every region, for the viewer)
"""
import json
import os
import random
import time
from collections import defaultdict
from pathlib import Path

import igraph as ig
import numpy as np
import polars as pl

ROOT = Path(os.environ.get("BIO_ROOT", "."))
OUT = ROOT / "data/processed"
SEED = 7
random.seed(SEED)
ig.set_random_number_generator(random)
rng = np.random.default_rng(SEED)

C0 = 10.0    # length per n^(1/3); bigger = airier regions
ORGANELLES = ["nucleus", "mitochondrion", "er", "golgi", "lysosome", "vesicle", "peroxisome"]
POSTHOC = {"metabolite", "chemical"}


def hdr(t):
    print(f"\n===== {t} =====", flush=True)


comp = pl.read_parquet(OUT / "compartments.parquet").select("key", "label", "kind", "region")
edges = pl.read_parquet(OUT / "edges.parquet").select("src", "tgt", "cls")
keys = comp["key"].to_list()
region = dict(zip(keys, comp["region"]))
kind = dict(zip(keys, comp["kind"]))
label = dict(zip(keys, comp["label"]))
_c = comp.group_by("region").len()
count = dict(zip(_c["region"], _c["len"]))

# ---------------------------------------------------------------- geometry
def pack(radii, R, tries=1500):
    """Centres for balls of given radii inside a ball of radius R; ball 0 pinned at the origin."""
    k, lim = len(radii), 0.88 * R
    if radii[0] > lim:
        return None, False
    c = rng.normal(size=(k, 3))
    c = c / np.linalg.norm(c, axis=1, keepdims=True) * R * 0.5 * rng.random((k, 1))
    c[0] = 0
    for _ in range(tries):
        moved = 0.0
        for i in range(k):
            for j in range(i + 1, k):
                d = c[i] - c[j]
                dist = np.linalg.norm(d) + 1e-9
                need = radii[i] + radii[j] + 0.04 * R
                if dist < need:
                    push = (need - dist) / 2 * d / dist
                    if i == 0:
                        c[j] -= 2 * push
                    else:
                        c[i] += push
                        c[j] -= push
                    moved += need - dist
        for i in range(1, k):
            room = lim - radii[i]
            if room < 0:
                return None, False
            dist = np.linalg.norm(c[i])
            if dist > room:
                c[i] *= room / dist
                moved += dist - room
        if moved < 1e-3:
            return c, True
    return c, False


radii = np.array([C0 * max(count.get(r, 0), 1) ** (1 / 3) for r in ORGANELLES])
n_cell = sum(count.get(r, 0) for r in ORGANELLES + ["cytosol"])
R = 1.3 * C0 * n_cell ** (1 / 3)
for _ in range(60):
    cen, ok = pack(radii, R)
    if ok:
        break
    R *= 1.06
else:
    raise SystemExit("could not pack the organelles")

geo = {}
for r, c, rad in zip(ORGANELLES, cen, radii):
    geo[r] = dict(shape="ball", c=np.asarray(c, dtype=float), r=float(rad))
geo["cytosol"] = dict(shape="cyto", c=np.zeros(3), r=0.93 * R, r_in=1.06 * float(radii[0]))
geo["plasma_membrane"] = dict(shape="sphere", c=np.zeros(3), r=R)
geo["extracellular"] = dict(shape="shell", c=np.zeros(3), r_in=1.25 * R, r=1.75 * R)
ORG_SPHERES = [(geo[r]["c"], geo[r]["r"]) for r in ORGANELLES]


def project(reg, P):
    """Pull points into the legal volume of a region."""
    g = geo[reg]
    v = np.atleast_2d(P) - g["c"]
    r = np.linalg.norm(v, axis=1, keepdims=True) + 1e-9
    if g["shape"] == "sphere":
        return g["c"] + v / r * g["r"]
    lo = g.get("r_in", 0.0)
    hi = g["r"] * (1.0 if g["shape"] == "shell" else 0.95)
    out = g["c"] + v / r * np.clip(r, lo, hi)
    if reg == "cytosol":
        for _ in range(3):
            for oc, orad in ORG_SPHERES:
                w = out - oc
                d = np.linalg.norm(w, axis=1, keepdims=True) + 1e-9
                out = np.where(d < orad * 1.06, oc + w / d * orad * 1.06, out)
            r = np.linalg.norm(out, axis=1, keepdims=True) + 1e-9
            out = out / r * np.clip(r, lo, hi)
    return out


def place(reg, P):
    """Radial equalisation: keep each node's direction, spread radii evenly through the region."""
    n = len(P)
    Q = P - np.median(P, axis=0)
    r = np.linalg.norm(Q, axis=1)
    d = Q / np.maximum(r, 1e-9)[:, None]
    zero = r < 1e-9
    if zero.any():
        z = rng.normal(size=(int(zero.sum()), 3))
        d[zero] = z / np.linalg.norm(z, axis=1, keepdims=True)
    frac = (np.argsort(np.argsort(r)) + 0.5) / n
    g = geo[reg]
    if g["shape"] == "sphere":
        rad = g["r"] * (1 + rng.normal(0, 0.008, n))
    else:
        lo = g.get("r_in", 0.0)
        hi = g["r"] * (1.0 if g["shape"] == "shell" else 0.95)
        rad = (lo ** 3 + frac * (hi ** 3 - lo ** 3)) ** (1 / 3)
    return project(reg, g["c"] + d * rad[:, None])


def force3d(ks, pairs):
    n = len(ks)
    if n < 5:
        return rng.normal(size=(n, 3))
    idx = {k: i for i, k in enumerate(ks)}
    g = ig.Graph(n=n, edges=[(idx[a], idx[b]) for a, b in pairs])
    niter = 500 if n < 2500 else 300
    try:
        lay = g.layout_fruchterman_reingold(dim=3, niter=niter)
    except TypeError:
        lay = g.layout("fr3d", niter=niter)
    return np.array(lay.coords, dtype=float)


hdr("Geometry")
print(f"cell radius {R:.0f}; nucleus r={radii[0]:.0f} at origin")
for r, c, rad in zip(ORGANELLES[1:], cen[1:], radii[1:]):
    print(f"  {r:14s} r={rad:5.0f}  centre=({c[0]:7.0f},{c[1]:7.0f},{c[2]:7.0f})")

# ---------------------------------------------------------------- force layout per region
pairs_by = defaultdict(set)
for s, t in edges.select("src", "tgt").iter_rows():
    if (s != t and s in region and t in region and region[s] == region[t]
            and kind[s] not in POSTHOC and kind[t] not in POSTHOC):
        pairs_by[region[s]].add((min(s, t), max(s, t)))

laid = defaultdict(list)
for k in keys:
    if kind[k] not in POSTHOC:
        laid[region[k]].append(k)

hdr("Force layout per region")
pos = {}
for reg, ks in sorted(laid.items(), key=lambda kv: len(kv[1])):
    t0 = time.time()
    Q = place(reg, force3d(ks, pairs_by[reg]))
    for k, p in zip(ks, Q):
        pos[k] = p
    print(f"{reg:16s} {len(ks):5d} nodes {len(pairs_by[reg]):6d} edges {time.time() - t0:6.1f}s", flush=True)

# ---------------------------------------------------------------- metabolites and chemicals
nb = defaultdict(list)
for s, t in edges.select("src", "tgt").iter_rows():
    nb[s].append(t)
    nb[t].append(s)

hdr("Metabolites and chemicals")
n_same = n_cross = n_none = 0
for k in keys:
    if kind[k] not in POSTHOC:
        continue
    reg, g = region[k], geo[region[k]]
    cand = [m for m in nb[k] if m in pos and kind[m] not in POSTHOC and region[m] == reg]
    if cand:
        n_same += 1
    else:
        cand = [m for m in nb[k] if m in pos and kind[m] not in POSTHOC]
        n_cross += bool(cand)
    if cand:
        pick = rng.choice(len(cand), size=min(5, len(cand)), replace=False)
        p = np.mean([pos[cand[i]] for i in pick], axis=0)
    else:
        n_none += 1
        p = g["c"] + rng.normal(size=3) * g["r"] * 0.4
    p = p + rng.normal(size=3) * g["r"] * 0.02
    pos[k] = project(reg, p)[0]
print(f"placed by same-region partners: {n_same}; by cross-region partners: {n_cross}; no partners: {n_none}")

# ---------------------------------------------------------------- write
X = np.array([pos[k] for k in keys], dtype=float)
assert np.isfinite(X).all()
pl.DataFrame({"key": keys, "label": comp["label"], "kind": comp["kind"], "region": comp["region"],
              "x": X[:, 0].astype("float32"), "y": X[:, 1].astype("float32"),
              "z": X[:, 2].astype("float32")}).write_parquet(OUT / "layout.parquet")
js = {"cell_radius": float(R), "regions": {
    r: {"shape": g["shape"], "center": [float(x) for x in g["c"]], "radius": float(g["r"]),
        "r_in": float(g.get("r_in", 0.0)), "n": int(count.get(r, 0))} for r, g in geo.items()}}
(OUT / "regions.json").write_text(json.dumps(js, indent=1), encoding="utf-8")

# ---------------------------------------------------------------- sanity
hdr("Sanity: edge lengths (units as above)")
same, cross = [], []
for s, t in edges.select("src", "tgt").iter_rows():
    if s in pos and t in pos and s != t and kind[s] not in POSTHOC and kind[t] not in POSTHOC:
        (same if region[s] == region[t] else cross).append(float(np.linalg.norm(pos[s] - pos[t])))
by = defaultdict(list)
for k in keys:
    if kind[k] not in POSTHOC:
        by[region[k]].append(k)
regs = [r for r in by if len(by[r]) > 1]
base = []
for _ in range(50000):
    L = by[regs[rng.integers(len(regs))]]
    i, j = rng.integers(len(L), size=2)
    if i != j:
        base.append(float(np.linalg.norm(pos[L[i]] - pos[L[j]])))
ms, mb = float(np.median(same)), float(np.median(base))
print(f"median length, same-region edges:        {ms:7.1f}  ({len(same)} edges)")
print(f"median length, random same-region pairs: {mb:7.1f}   -> edges are {ms / mb:.2f}x a random pair")
print(f"median length, cross-region edges:       {float(np.median(cross)):7.1f}  ({len(cross)} edges)")

hdr("Sanity: mTOR figure pairs (distance)")
lab2k = {}
for k in keys:
    if kind[k] in ("protein", "complex"):
        lab2k.setdefault(label[k].upper(), k)
for a, b in [("RHEB", "MTORC1"), ("RPTOR", "MTORC1"), ("TSC2", "RHEB"), ("AKT1", "TSC2"),
             ("PRKAA1", "TSC2"), ("EGFR", "AKT1")]:
    if a in lab2k and b in lab2k:
        d = np.linalg.norm(pos[lab2k[a]] - pos[lab2k[b]])
        print(f"  {a:8s} - {b:8s} {d:7.1f}   ({region[lab2k[a]]} / {region[lab2k[b]]})")
print("\nWrote data/processed/layout.parquet and regions.json")
