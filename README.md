# Cell interaction map

An interactive 3D map of signed, directed molecular interactions in a human cell.
Nodes sit in coarse cellular compartments; click a node to draw its edges
(green = stimulating, red = inhibiting ⊣, amber = mixed, gray = unsigned).

**Live viewer:** https://alirezahekmati.github.io/cell-interaction-map/

## Read this first
- Compartments are coarse and heuristic (HPA, UniProt, Human-GEM, plus a few manual overrides
  from the author's knowledge of mTOR signalling). Positions are schematic, not physical.
- A line means "reported in a database", not "experimentally verified".
- Metabolite regulation is sparse in signalling databases (e.g. AMP→AMPK is not curated in SIGNOR).

## Data sources
- SIGNOR (CC BY 4.0): https://signor.uniroma2.it
- OmniPath: https://omnipathdb.org. An aggregator whose resources have different licences.
  The default export contains entries with academic-use restrictions, so treat the merged
  data as academic/non-commercial unless you rebuild with a stricter licence filter.
- Human Protein Atlas, UniProt, Human-GEM: compartments and metabolite classes.
  Check each project's terms before reusing.

## Rebuild
Scripts in `scripts/` (run in order 01 → 07, Python via `uv`) produce `data/export/`.
Copy that into `viewer/public/data/`, then in `viewer/`: `npm install && npm run build`.
The site is built into `docs/` and served by GitHub Pages.
