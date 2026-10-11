import './style.css';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// raw hex colours end to end (custom shader + basic materials)
THREE.ColorManagement.enabled = false;

const BASE = import.meta.env.BASE_URL;
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
async function get(f, type = 'json') {
  let last;
  for (let attempt = 0; attempt < 3; attempt++) {
    try {
      // retries bypass the HTTP cache: some browsers return an empty body on the first request
      const r = await fetch(BASE + 'data/' + f, attempt ? { cache: 'reload' } : undefined);
      if (!r.ok) throw new Error(`${f}: HTTP ${r.status}`);
      const buf = await r.arrayBuffer();
      if (!buf.byteLength) throw new Error(`${f}: empty response`);
      return type === 'json' ? JSON.parse(new TextDecoder().decode(buf)) : buf;
    } catch (err) {
      last = err;
      await new Promise((res) => setTimeout(res, 300 * (attempt + 1)));
    }
  }
  throw last;
}

const KIND_COLORS = {
  protein: '#4dabf7', complex: '#b197fc', family: '#9775fa', metabolite: '#ffd43b',
  chemical: '#ff8787', phenotype: '#69db7c', mirna: '#f783ac', stimulus: '#ffa94d',
  antibody: '#ced4da', fusion_protein: '#66d9e8', ncrna: '#e599f7',
};
const MAXE = 15000;     // max edges drawn at once
const ROWS = 300;       // max rows in the detail list
const CAP = 6;          // arrow / T-bar size in world units

const DRUG_COLORS = { approved: '#22d3ee', approvedc: '#67e8b0', Tclin: '#67e8b0', trial: '#5b8def', earlier: '#3b5fa8', Tchem: '#ffa94d', Tbio: '#8a93a3', Tdark: '#4a5468', none: '#363d4a' };
const REG_SRC = new Set(['FDA', 'DailyMed', 'EMA', 'EMA Human Drugs', 'PMDA']);
function drugClass(n, drugs) {
  if (n.t !== 'protein') return 'none';
  const r = drugs[n.k.slice(2)];
  if (!r) return 'none';
  if (r.d && r.d.length) {
    if (r.d.some((x) => x[3] === 'APPROVAL' && String(x[6]).split(',').some((s) => REG_SRC.has(s)))) return 'approved';
    if (r.s === 'APPROVAL') return 'approvedc';
    if (r.o) return 'trial';
    if (r.s) return 'earlier';
  }
  if (r.t === 'Tclin') return 'approvedc';
  return r.t || 'none';
}

// each sidebar section (an h2 plus what follows it) becomes an open/close panel
const OPEN_PANELS = new Set(['Edge layers', 'Disease targets']);
function makePanels() {
  const side = $('side');
  let cur = null;
  for (const el of [...side.children]) {
    if (el.tagName === 'H2') {
      cur = document.createElement('details');
      cur.className = 'panel';
      cur.open = OPEN_PANELS.has(el.textContent.trim());
      const sm = document.createElement('summary');
      sm.textContent = el.textContent;
      cur.append(sm);
      side.insertBefore(cur, el);
      el.remove();
    } else if (cur) cur.append(el);
  }
}

async function getEdges() {
  const bin = atob(await get('edges_b64.json'));
  const u = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  return u.buffer;
}

async function main() {
  makePanels();
  const [meta, nodes, regions, buf, drugs] = await Promise.all([
    get('meta.json'), get('nodes.json'), get('regions.json'), getEdges(), get('node_drugs_ot.json').catch(() => ({}))]);
  const N = nodes.length, E = meta.n_edges;
  if (buf.byteLength !== E * 20) throw new Error(`edges.bin is ${buf.byteLength} bytes, expected ${E * 20}`);

  // ---------------------------------------------------------------- edge accessors (20-byte records)
  const U32 = new Uint32Array(buf), U8 = new Uint8Array(buf), U16 = new Uint16Array(buf);
  const eSrc = (e) => U32[e * 5], eTgt = (e) => U32[e * 5 + 1];
  const eCls = (e) => U8[e * 20 + 8], eNet = (e) => U8[e * 20 + 9], eDir = (e) => U8[e * 20 + 10];
  const eEv = (e) => U16[e * 10 + 6], eRes = (e) => U16[e * 10 + 7];

  // adjacency (CSR): incident edge ids per node
  const adjOff = new Uint32Array(N + 1);
  for (let e = 0; e < E; e++) { adjOff[eSrc(e) + 1]++; adjOff[eTgt(e) + 1]++; }
  for (let i = 0; i < N; i++) adjOff[i + 1] += adjOff[i];
  const fillp = adjOff.slice(0, N), inc = new Uint32Array(adjOff[N]);
  for (let e = 0; e < E; e++) { inc[fillp[eSrc(e)]++] = e; inc[fillp[eTgt(e)]++] = e; }

  // ---------------------------------------------------------------- node data
  const pos = new Float32Array(N * 3), baseSize = new Float32Array(N);
  const colReg = new Float32Array(N * 3), colKind = new Float32Array(N * 3), colDrug = new Float32Array(N * 3);
  const lower = new Array(N), tmp = new THREE.Color();
  const kindCount = {}, regCount = {};
  nodes.forEach((n, i) => {
    pos[i * 3] = n.x; pos[i * 3 + 1] = n.y; pos[i * 3 + 2] = n.z;
    baseSize[i] = Math.min(16, 2.2 + 0.8 * Math.sqrt(Math.max(n.d, 1)));
    tmp.set(meta.region_colors[n.r] || '#adb5bd'); colReg.set([tmp.r, tmp.g, tmp.b], i * 3);
    tmp.set(KIND_COLORS[n.t] || '#adb5bd'); colKind.set([tmp.r, tmp.g, tmp.b], i * 3);
    tmp.set(DRUG_COLORS[drugClass(n, drugs)]); colDrug.set([tmp.r, tmp.g, tmp.b], i * 3);
    lower[i] = n.l.toLowerCase();
    kindCount[n.t] = (kindCount[n.t] || 0) + 1;
    regCount[n.r] = (regCount[n.r] || 0) + 1;
  });

  // ---------------------------------------------------------------- state
  const S = {
    cls: new Set(['regulation']), net: new Set(meta.net), minEv: 1, minRes: 1, scope: 'sel',
    kinds: new Set(Object.keys(kindCount).filter((k) => k !== 'chemical')),
    regions: new Set(Object.keys(meta.region_colors)), noise: false,
    selected: new Set(), forced: new Set(), highlight: null,
    colorby: 'region', sizeMul: 1, labels: true, shells: true, shown: 0,
  };
  const vis = (i) => {
    if (S.forced.has(i)) return true;
    const n = nodes[i];
    if (n.h === 2 && !S.noise) return false;
    return S.kinds.has(n.t) && S.regions.has(n.r);
  };

  // ---------------------------------------------------------------- three.js scene
  const view = $('view'), canvas = $('c');
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
  renderer.setClearColor(0x0b0e14);
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(50, 1, 1, 40000);
  const R = regions.cell_radius;
  camera.position.set(0, R * 0.6, R * 3.2);
  const controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true; controls.dampingFactor = 0.08;
  controls.minDistance = 5; controls.maxDistance = R * 12;
  function resize() {
    const w = view.clientWidth, h = view.clientHeight;
    renderer.setSize(w, h, false);
    camera.aspect = w / Math.max(h, 1); camera.updateProjectionMatrix();
  }
  new ResizeObserver(resize).observe(view); resize();

  // nodes: one Points object with a tiny shader (per-node size, colour, alpha, selection ring)
  const aCol = new Float32Array(N * 3), aSize = new Float32Array(N), aAlpha = new Float32Array(N), aSel = new Float32Array(N);
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  const bCol = new THREE.BufferAttribute(aCol, 3), bSize = new THREE.BufferAttribute(aSize, 1);
  const bAlpha = new THREE.BufferAttribute(aAlpha, 1), bSel = new THREE.BufferAttribute(aSel, 1);
  geom.setAttribute('aCol', bCol); geom.setAttribute('aSize', bSize);
  geom.setAttribute('aAlpha', bAlpha); geom.setAttribute('aSel', bSel);
  const mat = new THREE.ShaderMaterial({
    uniforms: { uScale: { value: 1 } },
    transparent: true, depthWrite: false,
    vertexShader: `
      attribute vec3 aCol; attribute float aSize; attribute float aAlpha; attribute float aSel;
      varying vec3 vCol; varying float vAlpha; varying float vSel;
      uniform float uScale;
      void main() {
        vCol = aCol; vAlpha = aAlpha; vSel = aSel;
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        gl_PointSize = aSize * uScale / -mv.z;
        gl_Position = projectionMatrix * mv;
      }`,
    fragmentShader: `
      varying vec3 vCol; varying float vAlpha; varying float vSel;
      void main() {
        if (vAlpha < 0.01) discard;
        float d = length(gl_PointCoord - 0.5);
        if (d > 0.5) discard;
        vec3 c = vCol;
        if (vSel > 0.5 && d > 0.36) c = vec3(1.0);
        gl_FragColor = vec4(c, vAlpha);
      }`,
  });
  const points = new THREE.Points(geom, mat);
  points.frustumCulled = false;
  scene.add(points);

  // compartment outlines + their name labels
  const shellGroup = new THREE.Group(); scene.add(shellGroup);
  const labelBox = $('labelbox'), regionLabels = [];
  for (const [name, g] of Object.entries(regions.regions)) {
    const col = new THREE.Color(meta.region_colors[name] || '#888');
    let at = null;
    if (g.shape === 'ball' || g.shape === 'sphere') {
      const geo = new THREE.SphereGeometry(g.radius, 40, 24);
      const solid = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: col, transparent: true,
        opacity: name === 'plasma_membrane' ? 0.03 : 0.06, depthWrite: false, side: THREE.BackSide }));
      const wire = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: col, wireframe: true,
        transparent: true, opacity: 0.09, depthWrite: false }));
      solid.renderOrder = -2; wire.renderOrder = -1;
      solid.position.set(...g.center); wire.position.set(...g.center);
      shellGroup.add(solid, wire);
      at = g.shape === 'ball' ? [g.center[0], g.center[1] + g.radius, g.center[2]] : [0, -g.radius, 0];
    } else if (name === 'cytosol') at = [0, 0.6 * R, 0];
    else if (name === 'extracellular') at = [0, 1.45 * R, 0];
    if (at) {
      const el = document.createElement('div'); el.className = 'rlabel';
      el.textContent = name.replace('_', ' '); labelBox.append(el);
      regionLabels.push({ el, at });
    }
  }

  // ---------------------------------------------------------------- edges (LineSegments, rebuilt on demand)
  const netRGB = meta.net.map((n) => new THREE.Color(meta.net_colors[n]));
  const dv = new THREE.Vector3(), pvec = new THREE.Vector3(), upv = new THREE.Vector3();
  let edgeObj = null;
  function drawEdges(list) {
    if (edgeObj) { scene.remove(edgeObj); edgeObj.geometry.dispose(); edgeObj = null; }
    list = list.slice(0, MAXE);
    if (!list.length) return;
    let segs = 0;
    for (const e of list) { segs += 1; if (eDir(e)) { const nt = eNet(e); segs += nt === 2 ? 1 : nt === 3 ? 3 : 2; } }
    const P = new Float32Array(segs * 6), C = new Float32Array(segs * 6); let o = 0;
    const seg = (ax, ay, az, bx, by, bz, c) => {
      P[o] = ax; P[o + 1] = ay; P[o + 2] = az; P[o + 3] = bx; P[o + 4] = by; P[o + 5] = bz;
      C[o] = c.r; C[o + 1] = c.g; C[o + 2] = c.b; C[o + 3] = c.r; C[o + 4] = c.g; C[o + 5] = c.b; o += 6;
    };
    for (const e of list) {
      const s = eSrc(e), t = eTgt(e), nt = eNet(e), c = netRGB[nt];
      const ax = pos[s * 3], ay = pos[s * 3 + 1], az = pos[s * 3 + 2];
      const bx = pos[t * 3], by = pos[t * 3 + 1], bz = pos[t * 3 + 2];
      if (!eDir(e)) { seg(ax, ay, az, bx, by, bz, c); continue; }
      dv.set(bx - ax, by - ay, bz - az);
      const len = dv.length() || 1; dv.divideScalar(len);
      const off = Math.min(len * 0.4, baseSize[t] * S.sizeMul * 0.5 + 0.5);
      const tx = bx - dv.x * off, ty = by - dv.y * off, tz = bz - dv.z * off;
      seg(ax, ay, az, tx, ty, tz, c);
      upv.set(0, 1, 0); if (Math.abs(dv.y) > 0.95) upv.set(1, 0, 0);
      pvec.crossVectors(dv, upv).normalize();
      const k = Math.min(CAP, len * 0.3);
      if (nt !== 2) {            // arrow head: stimulating, unsigned, mixed
        const hx = tx - dv.x * k, hy = ty - dv.y * k, hz = tz - dv.z * k;
        seg(tx, ty, tz, hx + pvec.x * k * 0.5, hy + pvec.y * k * 0.5, hz + pvec.z * k * 0.5, c);
        seg(tx, ty, tz, hx - pvec.x * k * 0.5, hy - pvec.y * k * 0.5, hz - pvec.z * k * 0.5, c);
      }
      if (nt === 2 || nt === 3) { // T-bar: inhibiting, mixed
        seg(tx + pvec.x * k * 0.6, ty + pvec.y * k * 0.6, tz + pvec.z * k * 0.6,
            tx - pvec.x * k * 0.6, ty - pvec.y * k * 0.6, tz - pvec.z * k * 0.6, c);
      }
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(P, 3));
    g.setAttribute('color', new THREE.BufferAttribute(C, 3));
    edgeObj = new THREE.LineSegments(g, new THREE.LineBasicMaterial({
      vertexColors: true, transparent: true, opacity: 0.9, depthWrite: false }));
    edgeObj.frustumCulled = false; edgeObj.renderOrder = 2;
    scene.add(edgeObj);
  }

  const edgePass = (e) => S.cls.has(meta.cls[eCls(e)]) && S.net.has(meta.net[eNet(e)])
    && eEv(e) >= S.minEv && eRes(e) >= S.minRes && vis(eSrc(e)) && vis(eTgt(e));
  function collectEdges() {
    if (!S.selected.size) return [];
    const set = new Set();
    const add = (n) => { for (let q = adjOff[n]; q < adjOff[n + 1]; q++) { const e = inc[q]; if (edgePass(e)) set.add(e); } };
    if (S.scope === 'among' && S.selected.size > 1) {
      S.selected.forEach((n) => {
        for (let q = adjOff[n]; q < adjOff[n + 1]; q++) {
          const e = inc[q];
          if (S.selected.has(eSrc(e)) && S.selected.has(eTgt(e)) && edgePass(e)) set.add(e);
        }
      });
      return [...set].sort((a, b) => eRes(b) - eRes(a) || eEv(b) - eEv(a));
    }
    S.selected.forEach(add);
    if (S.scope === 'nbh') {
      const nb = new Set(S.selected);
      set.forEach((e) => { nb.add(eSrc(e)); nb.add(eTgt(e)); });
      nb.forEach((n) => {
        for (let q = adjOff[n]; q < adjOff[n + 1]; q++) {
          const e = inc[q];
          if (!set.has(e) && nb.has(eSrc(e)) && nb.has(eTgt(e)) && edgePass(e)) set.add(e);
        }
      });
    }
    return [...set].sort((a, b) => eRes(b) - eRes(a) || eEv(b) - eEv(a));
  }

  // ---------------------------------------------------------------- node attribute refresh
  function refreshNodes() {
    const hl = S.highlight, cols = S.colorby === 'kind' ? colKind : S.colorby === 'drug' ? colDrug : colReg;
    let shown = 0;
    for (let i = 0; i < N; i++) {
      const v = vis(i), s = S.selected.has(i);
      const gs = S.group ? S.group.score.get(i) : undefined;
      aSize[i] = v ? baseSize[i] * S.sizeMul * (s ? 1.7 : 1) * (gs === undefined ? 1 : 0.8 + 3.2 * gs) : 0;
      aAlpha[i] = !v ? 0 : (hl && !hl.has(i) ? 0.08 : 1);
      aSel[i] = s ? 1 : 0;
      if (v) shown++;
    }
    aCol.set(cols);
    bCol.needsUpdate = bSize.needsUpdate = bAlpha.needsUpdate = bSel.needsUpdate = true;
    S.shown = shown;
  }

  // ---------------------------------------------------------------- labels
  let labelNodes = [], labelEls = [];
  function setLabels(list) {
    labelBox.querySelectorAll('.nlabel').forEach((x) => x.remove());
    labelEls = []; labelNodes = [];
    if (!S.labels || !S.selected.size) return;
    const nb = new Set(); list.forEach((e) => { nb.add(eSrc(e)); nb.add(eTgt(e)); });
    const rest = [...nb].filter((i) => !S.selected.has(i)).sort((a, b) => nodes[b].d - nodes[a].d).slice(0, 30);
    const selList = [...S.selected];
    if (selList.length > 60) {
      const gsc = (i) => (S.group && S.group.score.get(i)) || 0;
      selList.sort((a, b) => gsc(b) - gsc(a) || nodes[b].d - nodes[a].d);
      selList.length = 40;
    }
    labelNodes = [...selList, ...rest];
    for (const i of labelNodes) {
      const el = document.createElement('div');
      el.className = 'nlabel' + (S.selected.has(i) ? ' sel' : '');
      el.textContent = nodes[i].l; labelBox.append(el); labelEls.push(el);
    }
  }
  const pv = new THREE.Vector3();
  function place(el, x, y, z) {
    pv.set(x, y, z).project(camera);
    if (pv.z > 1 || pv.z < -1 || Math.abs(pv.x) > 1.05 || Math.abs(pv.y) > 1.05) { el.style.display = 'none'; return; }
    el.style.display = 'block';
    el.style.transform = `translate(${(pv.x * 0.5 + 0.5) * view.clientWidth}px,${(-pv.y * 0.5 + 0.5) * view.clientHeight}px) translate(-50%,-140%)`;
  }

  // ---------------------------------------------------------------- detail panel
    const TDL_TXT = { Tclin: 'has approved drugs acting on it', Tchem: 'small molecules bind it strongly, but no approved drug',
    Tbio: 'biology is annotated, but no qualifying small molecules', Tdark: 'little is known about it' };
  const AGN = { F: 'FDA', E: 'EMA', P: 'PMDA' };
  function drugBoxLegacy(i) {
    const n = nodes[i];
    if (n.t !== 'protein') return '';
    const acc = n.k.slice(2), r = drugs[acc];
    if (!r) return '';
    let h = '<div class="desc">';
    if (r.t) h += `<div><b>Druggability: ${esc(r.t)}</b> <span class="m">${esc(TDL_TXT[r.t] || '')}${r.f ? ' · ' + esc(r.f) : ''}</span></div>`;
    if (r.d && r.d.length) {
      h += `<div style="margin-top:4px"><b>Drugs (${r.n})</b></div>`;
      h += r.d.map(([name, act, ag]) => `<div>${esc(name)} <span class="m">${esc(act || 'action n/a')} · ${ag === '-' ? 'not on FDA/EMA/PMDA lists' : ag.split('').map((c) => AGN[c]).join(', ')}</span></div>`).join('');
      if (r.n > r.d.length) h += `<div class="m">+${r.n - r.d.length} more, see Pharos</div>`;
    } else h += '<div class="m">No drug with a mechanism-of-action entry in DrugCentral.</div>';
    h += `<div class="m"><a href="https://pharos.nih.gov/targets/${esc(acc)}" target="_blank" rel="noopener">Pharos</a> · DrugCentral mechanism-of-action drugs only; clinical-trial drugs are not included yet</div></div>`;
    return h;
  }
  const STAGE_TXT = { APPROVAL: 'approved', PREAPPROVAL: 'pre-approval', PHASE_4: 'phase 4', PHASE_3: 'phase 3', PHASE_2_3: 'phase 2/3',
    PHASE_2: 'phase 2', PHASE_1_2: 'phase 1/2', PHASE_1: 'phase 1', EARLY_PHASE_1: 'early phase 1', IND: 'IND filed', PRECLINICAL: 'preclinical', UNKNOWN: 'stage unknown' };
  const STAGE_ORDER = ['APPROVAL', 'PREAPPROVAL', 'PHASE_4', 'PHASE_3', 'PHASE_2_3', 'PHASE_2', 'PHASE_1_2', 'PHASE_1', 'EARLY_PHASE_1', 'IND', 'PRECLINICAL', 'UNKNOWN'];
  const REG_NAME = { FDA: 'FDA', DailyMed: 'FDA label', EMA: 'EMA', 'EMA Human Drugs': 'EMA', PMDA: 'PMDA' };
  function drugRow(x) {
    const [name, type, act, stage, nOn, onStage, srcs, wd] = x;
    const all = String(srcs || '').split(',').filter(Boolean);
    const reg = [...new Set(all.filter((s) => REG_NAME[s]).map((s) => REG_NAME[s]))];
    const cur = all.filter((s) => !REG_NAME[s] && s !== 'manual');
    const parts = [type || 'type n/a', act || 'action n/a', STAGE_TXT[stage] || stage];
    if (nOn) parts.push(`${nOn} ongoing trial${nOn > 1 ? 's' : ''}${onStage && onStage !== stage ? ' (up to ' + (STAGE_TXT[onStage] || onStage) + ')' : ''}`);
    if (stage === 'APPROVAL') parts.push(reg.length ? reg.join(', ') : cur.length ? 'approval listed by ' + cur.join(', ') + ' only' : 'source n/a');
    if (all.includes('manual')) parts.push('link curated by hand: acts through FKBP12');
    if (wd) parts.push('withdrawn in some products');
    return `<div>${esc(name)} <span class="m">${esc(parts.join(' · '))}</span></div>`;
  }
  function drugBox(i) {
    const n = nodes[i];
    if (n.t !== 'protein') return '';
    const acc = n.k.slice(2), r = drugs[acc];
    if (!r) return '';
    let h = '<div class="desc">';
    if (r.t) h += `<div><b>Druggability: ${esc(r.t)}</b> <span class="m">${esc(TDL_TXT[r.t] || '')}${r.f ? ' · ' + esc(r.f) : ''}</span></div>`;
    if (r.d && r.d.length) {
      const by = {};
      r.d.forEach((x) => { by[x[3]] = (by[x[3]] || 0) + 1; });
      const sumTxt = STAGE_ORDER.filter((k) => by[k]).map((k) => by[k] + ' ' + STAGE_TXT[k]).join(', ');
      const on = r.d.filter((x) => x[4] > 0).length;
      h += `<div style="margin-top:4px"><b>${r.d.length} drug${r.d.length > 1 ? 's' : ''}</b> <span class="m">${esc(sumTxt)}</span></div>`;
      h += on ? `<div class="m">${on} with trials still running${r.o ? ' (furthest: ' + esc(STAGE_TXT[r.o] || r.o) + ')' : ''}</div>`
              : '<div class="m">No trials currently recruiting or active.</div>';
      h += r.d.slice(0, 10).map(drugRow).join('');
      const rest = r.d.slice(10);
      if (rest.length) h += `<details><summary class="m">show ${rest.length} more</summary>${rest.map(drugRow).join('')}</details>`;
    } else h += '<div class="m">No clinical-stage drug is linked to this protein in Open Targets.</div>';
    h += `<div class="m"><a href="https://pharos.nih.gov/targets/${esc(acc)}" target="_blank" rel="noopener">Pharos</a> · Open Targets (CC0) · ongoing = recruiting, active, not yet recruiting or enrolling by invitation</div></div>`;
    return h;
  }
  let infoP = null;
  async function fillInfo(i) {
    const n = nodes[i];
    if (n.t !== 'protein') return;
    const acc = n.k.slice(2);
    let all;
    try { infoP = infoP || get('node_info.json'); all = await infoP; } catch (err) { infoP = null; return; }
    if (!(S.selected.size === 1 && S.selected.has(i))) return;
    const box = $('info'), d = all[acc];
    if (!box || !d) return;
    let h = '';
    if (d.p) h += `<div><b>${esc(d.p)}</b></div>`;
    if (d.a && d.a.length) h += `<div class="m">also known as: ${esc(d.a.join(' · '))}</div>`;
    if (d.f) h += `<div>${esc(d.f)}</div>`;
    if (d.d && d.d.length) {
      const items = d.d.map(([name, mim]) => mim
        ? `<a href="https://omim.org/entry/${esc(mim)}" target="_blank" rel="noopener">${esc(name)}</a>` : esc(name));
      h += `<div><b>Linked to disease:</b> ${items.join('; ')}${d.dm ? ` (+${d.dm} more)` : ''}</div>`;
    } else if (d.n) h += `<div class="m">Disease note: ${esc(d.n)}</div>`;
    h += `<div class="m">UniProt <a href="https://www.uniprot.org/uniprotkb/${esc(acc)}" target="_blank" rel="noopener">${esc(acc)}</a> · informational only, not medical advice</div>`;
    box.className = 'desc';
    box.innerHTML = h;
  }
  const DIRSYM = ['⇢', '→', '⊣', '±'];
  function edgeSentence(e, d) {
    const s = nodes[eSrc(e)].l, t = nodes[eTgt(e)].l, cls = meta.cls[eCls(e)], net = meta.net[eNet(e)];
    const mech = (d.mech[e] || '').split(';').filter(Boolean).slice(0, 3).join(', ');
    const via = mech ? ` (${mech})` : '';
    if (cls === 'binding') return `${s} and ${t} bind each other${via}.`;
    if (cls === 'membership') return `${s} has ${t} as a member.`;
    if (cls === 'substrate_product') return `${s} and ${t} are linked as substrate and product${via}.`;
    if (cls === 'transport') return `${s} and ${t}: relocalisation${via}.`;
    const verb = { stim: 'activates', inhib: 'inhibits', mixed: 'activates or inhibits (reports conflict)',
      unsigned: 'regulates (direction not given)' }[net];
    return `${s} ${verb} ${t}${via}.`;
  }
  function nodeSummary(i) {
    const cnt = { out: [0, 0, 0], in: [0, 0, 0] };            // [activate, inhibit, mixed/unsigned]
    const top = { oS: [], oI: [], iS: [], iI: [] };
    for (let q = adjOff[i]; q < adjOff[i + 1]; q++) {
      const e = inc[q];
      if (meta.cls[eCls(e)] !== 'regulation' || !eDir(e)) continue;
      const nt = eNet(e), out = eSrc(e) === i, other = out ? eTgt(e) : eSrc(e);
      const k = nt === 1 ? 0 : nt === 2 ? 1 : 2;
      cnt[out ? 'out' : 'in'][k]++;
      if (k < 2) top[(out ? 'o' : 'i') + (k ? 'I' : 'S')].push([other, eRes(e), eEv(e)]);
    }
    const names = (a) => a.sort((x, y) => y[1] - x[1] || y[2] - x[2]).slice(0, 4)
      .map(([n]) => `<a data-node="${n}">${esc(nodes[n].l)}</a>`).join(', ');
    const part = (c, tS, tI, label, vS, vI) => {
      const tot = c[0] + c[1] + c[2];
      if (!tot) return '';
      let h = `<div><b>${label}</b> ${tot}: ${c[0]} activated, ${c[1]} inhibited, ${c[2]} mixed or unsigned.</div>`;
      const bits = [];
      if (tS.length) bits.push(`${vS} ${names(tS)}`);
      if (tI.length) bits.push(`${vI} ${names(tI)}`);
      return h + (bits.length ? `<div class="m">Strongest evidence: ${bits.join('; ')}.</div>` : '');
    };
    const h = part(cnt.out, top.oS, top.oI, 'Regulates', 'activates', 'inhibits')
      + part(cnt.in, top.iS, top.iI, 'Regulated by', 'activated by', 'inhibited by');
    return h ? `<div class="desc">${h}<div class="m">All regulation edges, ignoring the filters.</div></div>` : '';
  }
  let details = null;
  function renderDetail(list) {
    const el = $('detail');
    if (!S.selected.size) { el.hidden = true; return; }
    el.hidden = false;
    const sel = [...S.selected];
    let h = chipsHtml(sel);
    if (sel.length === 1) {
      const n = nodes[sel[0]];
      h += `<h2>${esc(n.l)}</h2><div class="sub">${esc(n.t)} · ${esc(n.r)}${n.u ? ' · UniProt ' + esc(n.u) : ''}${n.c ? ' · ' + esc(n.c) : ''}</div>`;
      h += `<div class="sub">${n.d} edges in total (${n.dr} regulation)${n.m ? ' · ' + n.m + ' members' : ''}</div>`;
      h += '<div id="info"></div>';
      h += drugBox(sel[0]);
      h += nodeSummary(sel[0]);
    } else {
      h += `<h2>${sel.length} nodes selected</h2><div class="sub">${sel.slice(0, 8).map((i) => esc(nodes[i].l)).join(', ')}${sel.length > 8 ? '…' : ''}</div>`;
    }
    h += `<div class="sub">${list.length} edges match the filters${list.length > MAXE ? ` (drawing the top ${MAXE})` : ''}${list.length > ROWS ? ` · listing the top ${ROWS} by resources` : ''}</div><div style="margin-top:8px">`;
    for (const e of list.slice(0, ROWS)) {
      const s = eSrc(e), t = eTgt(e), nt = eNet(e);
      h += `<div class="row" data-e="${e}"><div><a data-node="${s}">${esc(nodes[s].l)}</a> <b style="color:${meta.net_colors[meta.net[nt]]}">${eDir(e) ? DIRSYM[nt] : '—'}</b> <a data-node="${t}">${esc(nodes[t].l)}</a></div>`
        + `<div class="m">${meta.cls[eCls(e)]} · ${eEv(e)} ev · ${eRes(e)} res</div><div class="more" hidden></div></div>`;
    }
    el.innerHTML = h + '</div>';
    if (sel.length === 1) fillInfo(sel[0]);
  }
  $('detail').addEventListener('click', async (ev) => {
    const rg = ev.target.closest('button[data-rmgroup]');
    if (rg) { clearGroup(); return; }
    const rm = ev.target.closest('button[data-rm]');
    if (rm) { const i = +rm.dataset.rm; S.selected.delete(i); S.forced.delete(i); S.selected = new Set(S.selected); update(); return; }
    if (ev.target.closest('button[data-clear]')) { clearSel(); return; }
    const a = ev.target.closest('a[data-node]');
    if (a) { selectOnly(+a.dataset.node, true); return; }
    const row = ev.target.closest('.row');
    if (!row) return;
    const more = row.querySelector('.more');
    if (!more.hidden) { more.hidden = true; return; }
    more.hidden = false; more.textContent = 'loading…';
    try { details = details || await get('edge_details.json'); } catch (err) { more.textContent = String(err); return; }
    const e = +row.dataset.e, ns = details.ns[e];
    const pm = (details.pmids[e] || '').split(';').filter(Boolean)
      .map((p) => `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(p)}/" target="_blank" rel="noopener">${esc(p)}</a>`).join(', ');
    more.innerHTML = `<div><b>${esc(edgeSentence(e, details))}</b></div><div>evidence: stimulating ${ns[0]} · inhibiting ${ns[1]} · unsigned ${ns[2]}</div>`
      + `<div>resources: ${esc((details.res[e] || '').split(';').join(', ') || '—')}</div>`
      + `<div>mechanisms: ${esc((details.mech[e] || '').split(';').join(', ') || '—')}</div>`
      + `<div>PMIDs (${details.npm[e]} total): ${pm || '—'}</div>`;
  });

  // ---------------------------------------------------------------- central update
  let lastDrawn = 0;
  function update() {
    const list = collectEdges();
    drawEdges(list);
    lastDrawn = Math.min(list.length, MAXE);
    if (S.selected.size) {
      const hl = new Set(S.selected);
      for (const e of list.slice(0, MAXE)) { hl.add(eSrc(e)); hl.add(eTgt(e)); }
      S.highlight = hl;
    } else S.highlight = null;
    refreshNodes(); setLabels(list); renderDetail(list);
    $('status').textContent = `${S.shown.toLocaleString()} of ${N.toLocaleString()} nodes shown · ${lastDrawn.toLocaleString()} edges drawn`;
  }
  function selectOnly(i, fly) {
    remember();
    if (!vis(i)) S.forced.add(i);
    S.selected = new Set([i]); S.group = null;
    update();
    if (fly) flyTo(new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]), 120);
  }
  const hist = [];
  const backBtn = document.createElement('button');
  backBtn.id = 'backbtn'; backBtn.textContent = '↩ Back'; backBtn.title = 'return to the previous selection (Backspace)'; backBtn.hidden = true;
  view.append(backBtn);
  function sameSel(a) {
    return a.group === S.group && a.sel.size === S.selected.size && [...S.selected].every((x) => a.sel.has(x));
  }
  function remember() {
    if (!S.selected.size) return;
    const last = hist[hist.length - 1];
    if (last && sameSel(last)) return;
    hist.push({ sel: new Set(S.selected), forced: new Set(S.forced), group: S.group, scope: S.scope, d: curD });
    if (hist.length > 30) hist.shift();
    backBtn.hidden = false;
  }
  function goBack() {
    const h = hist.pop();
    backBtn.hidden = hist.length === 0;
    if (!h) return;
    S.selected = h.sel; S.forced = h.forced; S.group = h.group; S.scope = h.scope; curD = h.d;
    $('scope').value = h.scope;
    if (h.d) $('dq').value = h.d.name;
    update();
    if (S.selected.size) frameSelection();
  }
  backBtn.addEventListener('click', goBack);
  window.addEventListener('keydown', (ev) => {
    const t = document.activeElement && document.activeElement.tagName;
    if (ev.key === 'Backspace' && t !== 'INPUT' && t !== 'SELECT' && t !== 'TEXTAREA') { ev.preventDefault(); goBack(); }
  });
  function clearSel() { remember(); S.selected = new Set(); S.forced.clear(); S.group = null; update(); }
  function addToSel(i, fly) {
    remember();
    if (!vis(i)) S.forced.add(i);
    S.selected = new Set([...S.selected, i]);
    update();
    if (fly) frameSelection();
  }
  function frameSelection() {
    const ids = [...S.selected];
    if (!ids.length) return;
    const c = new THREE.Vector3();
    ids.forEach((i) => c.add(new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2])));
    c.divideScalar(ids.length);
    let rad = 0;
    ids.forEach((i) => { rad = Math.max(rad, c.distanceTo(new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]))); });
    flyTo(c, Math.max(120, rad * 2.6));
  }
  const chipsHtmlOld = (sel) => '<div class="chips">'
    + sel.map((i) => `<span class="chip">${esc(nodes[i].l)}<button data-rm="${i}" title="remove from selection">×</button></span>`).join('')
    + (sel.length > 1 ? '<button class="clearall" data-clear="1">clear all</button>' : '') + '</div>';

  // ---------------------------------------------------------------- camera tween
  let tween = null;
  function flyTo(target, dist) {
    const dir = camera.position.clone().sub(controls.target).normalize();
    tween = { t0: performance.now(), dur: 700, fromT: controls.target.clone(), toT: target.clone(),
      fromP: camera.position.clone(), toP: target.clone().addScaledVector(dir, dist) };
  }

  // ---------------------------------------------------------------- picking + pointer handling
  function pick(mx, my) {
    const W = view.clientWidth, H = view.clientHeight;
    const k = H / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov / 2)));
    camera.updateMatrixWorld();
    let best = -1, bs = 1e9;
    for (let i = 0; i < N; i++) {
      const sz = aSize[i]; if (sz === 0) continue;
      pv.set(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
      const dist = camera.position.distanceTo(pv);
      pv.project(camera);
      if (pv.z > 1 || pv.z < -1) continue;
      const dx = (pv.x * 0.5 + 0.5) * W - mx, dy = (-pv.y * 0.5 + 0.5) * H - my;
      const r = Math.max(0.5 * sz * k / dist, 3), d = Math.hypot(dx, dy);
      if (d <= r + 4 && d - r < bs) { bs = d - r; best = i; }
    }
    return best;
  }
  const rel = (ev) => { const r = canvas.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
  let down = null, hoverRaf = 0, hoverEv = null;
  canvas.addEventListener('pointerdown', (ev) => { down = { x: ev.clientX, y: ev.clientY }; $('tip').hidden = true; });
  canvas.addEventListener('pointerup', (ev) => {
    if (!down) return;
    const moved = Math.hypot(ev.clientX - down.x, ev.clientY - down.y); down = null;
    if (moved > 5 || ev.button !== 0) return;
    const [mx, my] = rel(ev), i = pick(mx, my);
    if (i >= 0) {
      if (ev.shiftKey) { S.selected.has(i) ? S.selected.delete(i) : S.selected.add(i); S.selected = new Set(S.selected); update(); }
      else selectOnly(i, false);
    } else if (!ev.shiftKey) clearSel();
  });
  canvas.addEventListener('pointermove', (ev) => {
    if (ev.buttons) return;
    hoverEv = ev;
    if (hoverRaf) return;
    hoverRaf = requestAnimationFrame(() => {
      hoverRaf = 0;
      const [mx, my] = rel(hoverEv), i = pick(mx, my), tip = $('tip');
      if (i < 0) { tip.hidden = true; canvas.style.cursor = 'default'; return; }
      const n = nodes[i];
      tip.textContent = `${n.l} · ${n.t} · ${n.r} · ${n.d} edges`;
      tip.style.left = mx + 'px'; tip.style.top = my + 'px'; tip.hidden = false; canvas.style.cursor = 'pointer';
    });
  });
  canvas.addEventListener('pointerleave', () => { $('tip').hidden = true; });
  window.addEventListener('keydown', (ev) => { if (ev.key === 'Escape') { document.activeElement?.blur(); clearSel(); } });

  // ---------------------------------------------------------------- sidebar
  function addCheck(parent, text, checked, onChange, color, count) {
    const l = document.createElement('label'); l.className = 'chk';
    const i = document.createElement('input'); i.type = 'checkbox'; i.checked = checked;
    i.addEventListener('change', () => onChange(i.checked));
    l.append(i);
    if (color) { const s = document.createElement('span'); s.className = 'sw'; s.style.background = color; l.append(s); }
    const t = document.createElement('span'); t.textContent = text; l.append(t);
    if (count != null) { const c = document.createElement('span'); c.className = 'cnt'; c.textContent = count.toLocaleString(); l.append(c); }
    parent.append(l); return l;
  }
  const clsN = {}, netN = {};
  for (let e = 0; e < E; e++) { const c = meta.cls[eCls(e)], n = meta.net[eNet(e)]; clsN[c] = (clsN[c] || 0) + 1; netN[n] = (netN[n] || 0) + 1; }
  meta.cls.forEach((c) => addCheck($('cls'), c.replace('_', ' '), S.cls.has(c),
    (on) => { on ? S.cls.add(c) : S.cls.delete(c); update(); }, null, clsN[c] || 0));
  const NETTXT = { stim: '→ stimulating', inhib: '⊣ inhibiting', mixed: '± mixed', unsigned: '— unsigned' };
  meta.net.slice().sort((a, b) => ['stim', 'inhib', 'mixed', 'unsigned'].indexOf(a) - ['stim', 'inhib', 'mixed', 'unsigned'].indexOf(b))
    .forEach((n) => addCheck($('net'), NETTXT[n] || n, true,
      (on) => { on ? S.net.add(n) : S.net.delete(n); update(); }, meta.net_colors[n], netN[n] || 0));
  $('minEv').addEventListener('input', (ev) => { S.minEv = Math.max(1, +ev.target.value || 1); update(); });
  $('minRes').addEventListener('input', (ev) => { S.minRes = Math.max(1, +ev.target.value || 1); update(); });
  $('scope').addEventListener('change', (ev) => { S.scope = ev.target.value; update(); });

  Object.keys(kindCount).sort((a, b) => kindCount[b] - kindCount[a]).forEach((k) => addCheck($('kinds'),
    k.replace('_', ' '), S.kinds.has(k), (on) => { on ? S.kinds.add(k) : S.kinds.delete(k); update(); },
    KIND_COLORS[k] || '#adb5bd', kindCount[k]));
  $('noise').addEventListener('change', (ev) => { S.noise = ev.target.checked; update(); });

  $('shells').addEventListener('change', (ev) => { S.shells = ev.target.checked; shellGroup.visible = S.shells; });
  Object.keys(meta.region_colors).forEach((r) => {
    const l = addCheck($('regions'), r.replace('_', ' '), true,
      (on) => { on ? S.regions.add(r) : S.regions.delete(r); update(); }, meta.region_colors[r], regCount[r] || 0);
    const b = document.createElement('button'); b.textContent = '⌖'; b.title = 'fly to ' + r;
    b.addEventListener('click', (ev) => {
      ev.preventDefault();
      const g = regions.regions[r]; if (!g) return;
      flyTo(new THREE.Vector3(...g.center), g.shape === 'ball' ? g.radius * 3.2 : R * 2.6);
    });
    l.append(b);
  });
  $('size').addEventListener('input', (ev) => { S.sizeMul = +ev.target.value; update(); });
  $('colorby').addEventListener('change', (ev) => { S.colorby = ev.target.value; refreshNodes(); setLegend(); });
  const legend = document.createElement('div'); legend.id = 'legend';
  $('colorby').parentElement.after(legend);
  const DRUG_TXT = [['approved', 'Target of an approved drug (Pharos Tclin)'], ['drug', 'Other drug in DrugCentral'],
    ['Tchem', 'Small molecules bind it, no drug (Tchem)'], ['Tbio', 'Biology known, no compounds (Tbio)'],
    ['Tdark', 'Little known (Tdark)'], ['none', 'Not a protein, or no data']];
  DRUG_TXT.length = 0;
  DRUG_TXT.push(['approved', 'Approved drug (FDA / EMA / PMDA / label)'], ['approvedc', 'Approval only on curated lists (ATC/TTD), or Pharos Tclin'],
    ['trial', 'Drug in an ongoing trial, none approved'], ['earlier', 'Drug reached trials, none ongoing or approved'],
    ['Tchem', 'Small molecules bind it, no drug (Tchem)'], ['Tbio', 'Biology known, no compounds (Tbio)'],
    ['Tdark', 'Little known (Tdark)'], ['none', 'Not a protein, or no data']);
  function setLegend() {
    legend.innerHTML = S.colorby === 'drug'
      ? DRUG_TXT.map(([k, t]) => `<div class="leg"><span class="sw" style="background:${DRUG_COLORS[k]}"></span><span>${t}</span></div>`).join('') : '';
  }
  $('labels').addEventListener('change', (ev) => { S.labels = ev.target.checked; update(); });

  // search
  $('q').addEventListener('input', () => {
    const q = $('q').value.trim().toLowerCase(), box = $('results');
    box.innerHTML = '';
    if (q.length < 2) return;
    const m = [];
    for (let i = 0; i < N; i++) { const p = lower[i].indexOf(q); if (p >= 0) m.push([i, p === 0 ? 0 : 1]); }
    m.sort((a, b) => a[1] - b[1] || nodes[b[0]].d - nodes[a[0]].d);
    m.slice(0, 12).forEach(([i]) => {
      const n = nodes[i], b = document.createElement('button');
      b.innerHTML = `${esc(n.l)} <small>${esc(n.t)} · ${esc(n.r)} · ${n.d}</small>`;
      b.addEventListener('click', () => { addToSel(i, true); $('q').value = ''; $('results').innerHTML = ''; });
      box.append(b);
    });
  });
  $('q').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') $('results').querySelector('button')?.click(); });

  // ---------------------------------------------------------------- render loop
  // ---------------------------------------------------------------- disease panel (Open Targets, CC0)
  S.group = null;
  function chipsHtml(sel) {
    const g = S.group;
    const manual = g ? sel.filter((i) => !g.set.has(i)) : sel;
    let h = '<div class="chips">';
    if (g) h += `<span class="chip group">${esc(g.name)} · ${g.set.size} targets<button data-rmgroup="1" title="remove these targets">×</button></span>`;
    h += manual.slice(0, 40).map((i) => `<span class="chip">${esc(nodes[i].l)}<button data-rm="${i}" title="remove from selection">×</button></span>`).join('');
    if (manual.length > 40) h += `<span class="m">+${manual.length - 40} more</span>`;
    if (sel.length > 1) h += '<button class="clearall" data-clear="1">clear all</button>';
    return h + '</div>';
  }
  function clearGroup() {
    remember();
    if (!S.group) return;
    S.group.set.forEach((i) => { S.selected.delete(i); S.forced.delete(i); });
    S.selected = new Set(S.selected);
    S.group = null;
    $('dinfo').textContent = 'Disease selection removed.';
    update();
  }
  const DZ = { index: null, blob: null, off: null, low: null, loading: null };
  function loadDiseases() {
    if (DZ.loading) return DZ.loading;
    $('dinfo').textContent = 'loading disease data…';
    DZ.loading = (async () => {
      const [idx, b64] = await Promise.all([get('disease_index.json'), get('disease_targets_b64.json')]);
      const bin = atob(b64), u = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
      const off = new Array(idx.length);
      let o = 0;
      idx.forEach((r, k) => { off[k] = o; o += r[2]; });
      DZ.index = idx; DZ.blob = u; DZ.off = off;
      DZ.low = idx.map((r) => (r[1] + ' ' + r[3]).toLowerCase());
      $('dinfo').textContent = idx.length.toLocaleString() + ' diseases ready. Type a name above.';
    })().catch((err) => { DZ.loading = null; $('dinfo').textContent = 'Could not load disease data: ' + err.message; throw err; });
    return DZ.loading;
  }
  let curD = null;
  function applyDisease() {
    if (!curD) return;
    const minS = +$('dscore').value, cap = Math.max(1, Math.floor(+$('dcap').value) || 1);
    if (S.group) S.group.set.forEach((i) => { S.selected.delete(i); S.forced.delete(i); });
    const set = new Set(), score = new Map();
    for (let j = 0; j < curD.idxs.length && set.size < cap; j++) {
      if (curD.sc[j] < minS) break;
      set.add(curD.idxs[j]); score.set(curD.idxs[j], curD.sc[j]);
    }
    S.group = { name: curD.name, set, score };
    set.forEach((i) => { if (!vis(i)) S.forced.add(i); S.selected.add(i); });
    S.selected = new Set(S.selected);
    S.scope = 'among'; $('scope').value = 'among';
    $('dinfo').textContent = `${set.size} of ${curD.idxs.length} scored targets shown (score ≥ ${minS.toFixed(2)}; top score ${curD.sc[0].toFixed(2)}). Only targets present in this map are listed, and scores below 0.1 are not stored.`;
    update();
  }
  function pickDisease(k) {
    remember();
    const r = DZ.index[k], o = DZ.off[k], n = r[2];
    const idxs = new Array(n), sc = new Array(n);
    for (let j = 0; j < n; j++) {
      const b = (o + j) * 3;
      idxs[j] = DZ.blob[b] | (DZ.blob[b + 1] << 8);
      sc[j] = DZ.blob[b + 2] / 255;
    }
    curD = { id: r[0], name: r[1], idxs, sc };
    $('dres').innerHTML = ''; $('dq').value = r[1];
    $('dscore').max = String(Math.max(0.2, Math.ceil(sc[0] * 100) / 100));
    $('dscore').value = '0.1'; $('dscoreV').textContent = '0.10';
    $('dcap').max = String(n); $('dcap').value = String(Math.min(150, n));
    applyDisease();
    frameSelection();
  }
  $('dq').addEventListener('focus', () => { loadDiseases().catch(() => {}); });
  $('dq').addEventListener('input', async () => {
    const box = $('dres'), q = $('dq').value.trim().toLowerCase();
    box.innerHTML = '';
    if (q.length < 2) return;
    try { await loadDiseases(); } catch (e) { return; }
    if ($('dq').value.trim().toLowerCase() !== q) return;
    const m = [];
    for (let k = 0; k < DZ.index.length; k++) {
      const p = DZ.index[k][1].toLowerCase().indexOf(q);
      if (p === 0) m.push([k, 0]);
      else if (p > 0) m.push([k, 1]);
      else if (DZ.low[k].indexOf(q) >= 0) m.push([k, 2]);
    }
    m.sort((a, b) => a[1] - b[1] || DZ.index[b[0]][2] - DZ.index[a[0]][2]);
    if (!m.length) { box.textContent = 'No disease with scored targets in this map matches.'; return; }
    m.slice(0, 12).forEach(([k]) => {
      const b = document.createElement('button');
      b.innerHTML = `${esc(DZ.index[k][1])} <small>${DZ.index[k][2]} targets</small>`;
      b.addEventListener('click', () => pickDisease(k));
      box.append(b);
    });
  });
  $('dq').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') $('dres').querySelector('button')?.click(); });
  $('dscore').addEventListener('input', () => { $('dscoreV').textContent = (+$('dscore').value).toFixed(2); applyDisease(); });
  $('dcap').addEventListener('input', applyDisease);
  $('dclear').addEventListener('click', clearGroup);

  const dbs = new THREE.Vector2();
  function frame(now) {
    requestAnimationFrame(frame);
    if (tween) {
      const k = Math.min(1, (now - tween.t0) / tween.dur), e = k * k * (3 - 2 * k);
      controls.target.lerpVectors(tween.fromT, tween.toT, e);
      camera.position.lerpVectors(tween.fromP, tween.toP, e);
      if (k >= 1) tween = null;
    }
    controls.update();
    renderer.getDrawingBufferSize(dbs);
    mat.uniforms.uScale.value = dbs.y / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov / 2)));
    labelNodes.forEach((i, j) => place(labelEls[j], pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]));
    regionLabels.forEach((r) => {
      if (!S.shells) { r.el.style.display = 'none'; return; }
      place(r.el, r.at[0], r.at[1], r.at[2]);
    });
    renderer.render(scene, camera);
  }
  update();
  requestAnimationFrame(frame);
}

main().catch((e) => {
  console.error(e);
  document.body.insertAdjacentHTML('beforeend', `<pre style="position:fixed;top:0;left:0;right:0;margin:0;padding:12px;background:#300;color:#fbb;z-index:9">${esc(e.stack || e)}</pre>`);
});
