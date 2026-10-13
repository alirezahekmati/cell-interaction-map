# Cell interaction map

An interactive 3D map of about 9,700 human proteins, plus complexes and metabolites, and their signed, directed interactions. It is overlaid with drug, clinical-trial, tractability, safety and disease-association data. It is built for early drug-discovery triage: which targets matter for a disease, which of them have no drug yet but could be drugged, and what sits around them in the signalling network.

**Live viewer:** https://alirezahekmati.github.io/cell-interaction-map/

> A research and exploration tool, not medical advice. Please read [Limitations](#limitations) before drawing conclusions.

## Contents

- [What the map shows](#what-the-map-shows)
- [Quick start](#quick-start)
- [The screen](#the-screen)
- [Panels](#panels)
- [The details panel](#the-details-panel)
- [Workflows](#workflows)
- [Colours, levels and scores explained](#colours-levels-and-scores-explained)
- [Keyboard and accessibility](#keyboard-and-accessibility)
- [Limitations](#limitations)
- [Data sources and licences](#data-sources-and-licences)
- [Run it locally and rebuild the data](#run-it-locally-and-rebuild-the-data)

## What the map shows

- Each dot is a protein, complex, protein family, metabolite or other entity: about 12,200 nodes and 52,407 merged edges. All evidence for one pair, direction and edge layer is merged into one line.
- Dots sit in coarse cell compartments: extracellular, plasma membrane, cytosol, nucleus, mitochondrion, ER, Golgi, lysosome, vesicle and peroxisome. **Positions are schematic, not physical.**
- Lines are drawn only for what you select, so the map stays readable. Line symbols:
  - green with an arrowhead: stimulating
  - red with a T-bar (⊣): inhibiting
  - amber, arrow plus T-bar: mixed (reports conflict)
  - gray: no direction or sign given
- Edge layers: **regulation** (on by default), binding, membership (complex or family to member), substrate/product and transport.

## Quick start

1. Drag to rotate, right-drag to pan, scroll to zoom. Hover a dot for a quick label, and keep hovering for more.
2. Type a protein (for example `AKT1`) in the top search box and press Enter. Its edges draw, and the details panel opens on the right.
3. Open the **Disease targets** panel, type a disease (for example `asthma`) and pick it. The top targets are selected, bigger dots mean higher scores, and a summary card counts them by drug status.
4. Set **Display → Colour by → drug status** to see at a glance what has a drug.
5. Click **Table** (top right) for the current set as a sortable table, and **Copy link** to share exactly what you are looking at.
6. **Reset** (top right) clears everything. **Back** (top left, or Backspace) returns to the previous view.

## The screen

- **Left sidebar:** search boxes and collapsible panels (click a panel title to open or close it).
- **Centre:** the 3D map. Top left: the Back button, and a chip while a perturbation preview is running. Top right: **Table**, **Copy link** and **Reset**. Bottom left: how many nodes and edges are shown.
- **Right:** the details panel for what you selected.
- **Bottom of the map:** the table, when open.

Clicking a dot selects it. Shift-click adds or removes a dot. Clicking empty space clears the focus.

## Panels

### Search and selection
The top search box finds proteins and metabolites. Choosing a result **adds** it to the selection, so you can search several in a row. Each selected node is a chip in the details panel, and the **×** on a chip removes just that node.

### Drug search
Type a drug by generic name or trade name (for example `sirolimus` or `rapamune`). The card shows its type, best stage, the diseases it is approved for and in trials for, and every target it acts on in this map with the action. If a disease is open, it also says which of the drug's targets are among that disease's targets, with their scores. **Select its targets** adds them to the selection.

### Disease targets
Search a disease (data from Open Targets). The panel selects the highest-scoring targets that are present in this map.
- **Min score** and **Max targets** control the cut-off. Scores below 0.1 are not stored, and the list is capped at 2,000 targets per disease.
- The disease appears as one group chip with a **×** that removes the whole group. Dot size follows the score (or the hit score, see below).
- The **summary card** counts the shown targets by status. Each line is a button that filters the map.
  - approved drug for this disease
  - drug in trials for this disease
  - approved drug, but for other diseases (a repurposing signal)
  - drug in trials, only for other diseases
  - no clinical-stage drug (with a count of those scoring 0.6 or more)
  - linked to another target of the set, or not (using the edge layers you have on)
- Filters: **Show** (by the categories above), **Connections** (all, only connected, only unconnected), **Evidence** (genetic, animal model, known drug, or weak evidence only) and **Evidence strength** (0.2, 0.4 or 0.6).
- A disease counts a drug as approved or in trials for it when the drug's indication is that disease, one of its subtypes, or a broader non-therapeutic-area parent in the disease ontology.
- Click a target while a disease is open: its edges and details show on top, and the disease stays. Esc clears the focus first, and a second Esc clears the disease.

### Hit list (triage)
Ranks the open disease's targets with a transparent score. Six sliders set the weights: disease score, genetic evidence, tractability, clinical precedent (a **negative** weight favours undrugged targets), a penalty for common essential genes and a penalty for safety liabilities.
- Presets: **Balanced**, **Repurposing** (also filters to approved-elsewhere) and **First-in-class** (prefers undrugged, tractable targets with genetic evidence).
- Moving a slider switches the dots to **hit score** sizing, and the table gets a **Hit score** column.
- Formula: hit score = (w1 disease + w2 genetic + w3 tractability + w4 precedent − w5 essential − w6 safety) ÷ (w1 + w2 + w3 + |w4|), scaled 0 to 100 across the open disease's selected targets.
- It is a ranking aid for triage, not a prediction of success.

### Bridge proteins
Needs an open disease or two or more selected nodes. Lists proteins outside your set that link to two or more of its members through the edge layers you have on, with how many links and whether they activate, inhibit or are regulated by those members. **Add the top 10 to the view** draws the wiring. This shows how a set of targets is connected.

### Perturbation preview
Select one node and press **Inhibit it** or **Activate it** in the details panel (the button matching the protein's clinical drugs is marked "like its drugs"). The map shows which nodes could rise or fall, following the signs of directed regulation edges for 1 to 3 steps: inhibiting a node that inhibits another predicts that one rises.
- Orange ↑ rises, blue ↓ falls, purple ↕ conflicting. Everything else dims.
- The panel lists direct effects, further downstream effects, phenotypes reached, and which shown targets of an open disease are reached.
- **Draw the cascade (top 40)** selects the source plus the 40 strongest effects and draws the links among them. The table then shows the predicted change and the step.
- Mixed and unsigned links are not followed. This is a qualitative hypothesis from curated signs, not a simulation.

### Drug filters
Show only targets by clinical status, on the whole map and on disease sets.
- Presets: All, With a drug, Approved, Ongoing trials, No drug or trial. Tier checkboxes with counts give finer control (colours below).
- **Ongoing trials:** any, only with, only without.
- **Drug type:** a target must have a drug of a ticked type (small molecule, antibody, antibody drug conjugate, protein, enzyme, oligonucleotide and so on).
- **Tractability:** small molecule, antibody, or **tractable, undrugged**, which finds targets with no clinical-stage drug that look druggable.

### Edge layers, Sign and Nodes
- **Edge layers:** choose which layers can be drawn. **Sign:** stimulating, inhibiting, mixed, unsigned.
- **Min evidence rows** and **Min resources** hide weakly supported edges. **Scope:** selected node(s) only, also edges among neighbours, or only edges between selected nodes.
- **Nodes:** show or hide node kinds. Drug-like chemicals are hidden by default, and so are noise hubs such as calcium, sodium and hydron (tick **Show noise hubs** to include them).

### Compartments and Display
- **Compartments:** show or hide the cell outlines, filter by compartment, and fly to one with the ⌖ button.
- **Display:** node size, **Colour by** (compartment, node type, drug status or tractability), labels on selection, **colour-blind-safe edge colours** (blue and orange) and a **high-contrast interface**. The two accessibility options are remembered in your browser.

## The details panel

Select one node and the right panel shows:
- Name, type, compartment and IDs, and the UniProt function text with disease links (to OMIM).
- **Druggability** (Pharos level) and **drugs**: type, action, stage, ongoing trials, approval sources, and the diseases each drug is approved or in trials for.
- **Upstream handles:** for a target with a late-stage drug, upstream regulators that could act in the same direction. They are ranked by druggability and by how few other proteins they regulate. It is a rough guide, not selectivity data, and it appears only when the direction of the drugs is clear.
- **Tractability** for four modalities, **essentiality** (DepMap), reported safety liabilities and tool compounds with links.
- With a disease open: the target's score and category, and **evidence bars** showing why it scores.
- A summary of what it regulates and is regulated by, then every edge. Click an edge row for the evidence: resources, mechanisms and PubMed links.

## Workflows

**Find untried targets for a disease.** Pick the disease, press **First-in-class** in the Hit list panel, and open **Table**. Tick Drug filters → **tractable, undrugged**. Use **Evidence strength** to demand stronger genetic support. Download the CSV.

**Find repurposing candidates.** Pick the disease and press **Repurposing**. These are shown targets with an approved drug for other diseases. Click a target to see which drug and which diseases. Use **Drug search** on the drug to see its other targets and approved uses.

**Judge a single target.** Search it and read the details panel: drugs and their stages, tractability, whether it is a common essential gene, safety liabilities, tool compounds and, if a disease is open, the evidence bars.

**See how a disease's targets are wired together.** Open the disease, set **Scope** to "only edges between selected nodes" (done automatically), then open **Bridge proteins** and add the top 10.

**Look for a way to move an approved target indirectly.** Select a target with a clinical drug and open **Upstream handles** in the details panel. Candidates with narrower regulation and better druggability are listed first. Then try **Perturbation preview** on a candidate to see what else moves.

**Ask what a drug does across the network.** Use **Drug search**, then **Select its targets**, and set scope to "only edges between selected nodes".

**Share or export.** **Copy link** stores the disease, filters, selection, weights and preview in the address. **Download CSV** and **Copy for spreadsheet** (paste into a spreadsheet; useful if a download manager blocks downloads) in the table.

## Colours, levels and scores explained

**Colour by drug status**

| Colour | Meaning |
| --- | --- |
| cyan | approved drug with an FDA, EMA, PMDA or drug-label source |
| mint | approval listed only by curated lists (ATC, TTD), or Pharos level Tclin |
| blue | drug in an ongoing trial, none approved |
| darker blue | drug reached trials, none ongoing or approved |
| orange | Tchem: small molecules bind it strongly, no drug |
| gray | Tbio: biology annotated, no qualifying compounds |
| dark gray | Tdark: little known |
| very dark | not a protein, or no data |

**Ongoing trials** means a trial whose status is recruiting, active not recruiting, not yet recruiting or enrolling by invitation.

**Pharos levels (Tclin, Tchem, Tbio, Tdark)** classify how much is known about a protein as a drug target. Tclin targets have approved drugs acting on them.

**Tractability (per modality, 0 to 5):** 0 no evidence, 1 predicted tractable (lower confidence), 2 predicted tractable (high confidence), 3 phase 1 precedent, 4 advanced clinical precedent, 5 approved-drug precedent. A high-confidence antibody prediction is counted only if the protein is on the cell surface or secreted according to this map's compartments.

**Evidence strength:** each evidence type has a level from 0 to 15 (the score times 15). The strength choices 0.2, 0.4 and 0.6 correspond to levels 3, 6 and 9. The types are genetic associations, genetic literature, somatic mutations, known drugs, affected pathways, RNA expression, literature and animal models.

**Essential (DepMap):** "common essential" means the gene is flagged essential across cancer cell-line knockout screens. The percentage is the share of screened cell lines with a gene effect of −0.5 or lower.

**Safety liabilities** are reported adverse effects and potential effects of modulating the target, from curated sources. Assay readouts from ToxCast and ecotoxicology outcomes from AOP-Wiki are left out. This is not a clinical safety verdict.

## Keyboard and accessibility

- **Tab** moves through the controls with a visible focus ring. **Enter** in a search box picks the first result.
- **Esc** clears the focus, then the disease or selection. **Backspace** (outside text boxes) goes back.
- Line symbols (arrow versus T-bar) carry the sign as well as the colour, and the colour-blind-safe option changes the palette. Muted text meets WCAG-AA contrast on the dark background, and there is a high-contrast interface option.
- The camera fly-to is skipped if your system asks for reduced motion.
- The 3D view itself is visual only. The table, search, filters and details panel give access to the same information without the mouse.

## Limitations

- **Coverage:** the map holds about 9,700 proteins that appear in the SIGNOR and OmniPath signalling databases. A disease may have targets outside this set, and they are not shown. Signalling databases are biased toward well-studied pathways, and transporters, structural proteins and many enzymes are under-represented.
- **"Unconnected" or "no link"** means no signalling edge in this map under the layers you have on. It does not mean no biology.
- **A line means "reported in a database", not "experimentally verified".** Signs come from curated records. Where OmniPath records both stimulation and inhibition, only the sign its resources agree on is kept when they point one way. Some edges stay mixed.
- **Compartments are heuristic** (Human Protein Atlas, UniProt and Human-GEM, plus manual overrides for the mTOR pathway proteins), and positions are schematic. Edges between compartments are long by design.
- **Metabolite regulation is sparse** in signalling databases. For example, AMP activation of AMPK is not curated in SIGNOR, and the map cannot show it. Only human and unspecified-species evidence is used.
- **Drug links** come from ChEMBL through Open Targets. The links of sirolimus, everolimus and temsirolimus to MTOR (through FKBP12) were added by hand and are labelled as such. "Approved" can rest on curated lists (ATC, TTD) and not a regulator. These cases are marked in the drug list. Drugs withdrawn in some products are flagged.
- **An approved drug means clinical precedent** for modulating the target, in some indication and at some dose. It does not mean the target is safe or effective for your disease. Oncology drugs accept toxicity that a chronic-disease drug cannot.
- **Hit list, handles and perturbation preview** are triage aids. The weights are opinions, the "fewer regulated proteins" ranking is not selectivity data, and the perturbation preview ignores dose, timing, feedback and cell type.
- **It is a snapshot.** The data reflect the releases listed below and are not updated automatically.

## Data sources and licences

| Source | Used for | Licence |
| --- | --- | --- |
| [SIGNOR](https://signor.uniroma2.it) | signed, directed interactions, complexes and families | CC BY 4.0 |
| [OmniPath](https://omnipathdb.org) | additional interactions | an aggregator: resources have different licences, and its default export includes entries with academic-use restrictions. Treat the merged interaction data as academic or non-commercial unless you rebuild with a stricter licence filter. |
| [Open Targets Platform](https://platform.opentargets.org) (release 26.09) | drug-target links, drug types, clinical stages, trials, indications, disease associations, evidence types, tractability, safety events, DepMap essentiality, chemical probes | CC0 1.0 |
| [UniProt](https://www.uniprot.org) | protein names, function and disease text, subcellular location | CC BY 4.0 |
| [Pharos / TCRD](https://pharos.nih.gov) (NCATS) | target development levels | see the Pharos site |
| [Human Protein Atlas](https://www.proteinatlas.org) and [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) | compartments and metabolite classes | check each project's terms before reuse |

An earlier version also used DrugCentral (CC BY-SA 4.0). It is no longer part of the data shipped with the site.

Please cite the sources if you use the results. For Open Targets: *Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery* (2025).

## Run it locally and rebuild the data

The viewer is a static site built with [three.js](https://threejs.org) and [Vite](https://vitejs.dev). To run it locally:

    cd viewer
    npm install
    npm run dev

`npm run build` writes the site into `docs/`, which GitHub Pages serves (Settings → Pages → branch `main`, folder `/docs`). Add an empty `docs/.nojekyll` file after each build.

The raw data (`data/raw/`) and intermediate data (`data/processed/`, `data/export/`) are not in the repository. To rebuild them you need Python 3.12 with polars, pyarrow, pandas, httpx, python-igraph, numpy and omnipath (managed with [uv](https://docs.astral.sh/uv/)). Run the scripts from the project root:

1. `scripts/01_fetch.py` downloads most raw files. The SIGNOR export was downloaded by hand from the SIGNOR site.
2. Three more downloads, with curl:

        curl -fSL "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,gene_primary,cc_subcellular_location&query=(organism_id:9606)%20AND%20(reviewed:true)" -o data/raw/uniprot/human_subcellular.tsv
        curl -fSL "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,gene_primary,protein_name,cc_function,cc_disease&query=(organism_id:9606)%20AND%20(reviewed:true)" -o data/raw/uniprot/human_function_disease.tsv
        curl -fSL https://opendata.ncats.nih.gov/public/pharos/pharos400_tdls_canonicals.csv -o data/raw/pharos/pharos400_tdls_canonicals.csv

3. `scripts/11_fetch_opentargets.py` downloads the Open Targets datasets (resumable; add `--extras`-style dataset names with `--datasets`, and see `--help`).
4. Build in order: `04_build_graph.py`, `05_compartments.py`, `06_layout.py`, `07_export.py`, `09_node_info.py`, `13_build_ot_drugs.py`, `14_build_ot_diseases.py`, `15_build_indications.py`, `16_build_target_risk.py`, `17_build_disease_evidence.py`, `18_build_drug_names.py`. (`10_drugs.py` is the older DrugCentral build, not used by the site.)
5. `scripts/08_pack_viewer.py` copies the export into `viewer/public/data/` (edges are shipped as text so download managers do not block them).

The scripts write reports to `data/*.txt` that show what each step produced.
