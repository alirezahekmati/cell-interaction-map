import './style.css';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// raw hex colours end to end (custom shader + basic materials)
THREE.ColorManagement.enabled = false;

const BASE = import.meta.env.BASE_URL;
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
async function get(f, type = 'json') {
  const r = await fetch(BASE + 'data/' + f);
  if (!r.ok) throw new Error(`${f}: HTTP ${r.status}`);
  return type === 'json' ? r.json() : r.arrayBuffer();
}

const KIND_COLORS = {
  protein: '#4dabf7', complex: '#b197fc', family: '#9775fa', metabolite: '#ffd43b',
  chemical: '#ff8787', phenotype: '#69db7c', mirna: '#f783ac', stimulus: '#ffa94d',
  antibody: '#ced4da', fusion_protein: '#66d9e8', ncrna: '#e599f7',
};
const MAXE = 15000;     // max edges drawn at once
const ROWS = 300;       // max rows in the detail list
const CAP = 6;          // arrow / T-bar size in world units

async function main() {
  const [meta, nodes, regions, buf] = await Promise.all([
    get('meta.json'), get('nodes.json'), get('regions.json'), get('edges.bin', 'buf')]);
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
  const colReg = new Float32Array(N * 3), colKind = new Float32Array(N * 3);
  const lower = new Array(N), tmp = new THREE.Color();
  const kindCount = {}, regCount = {};
  nodes.forEach((n, i) => {
    pos[i * 3] = n.x; pos[i * 3 + 1] = n.y; pos[i * 3 + 2] = n.z;
    baseSize[i] = Math.min(16, 2.2 + 0.8 * Math.sqrt(Math.max(n.d, 1)));
    tmp.set(meta.region_colors[n.r] || '#adb5bd'); colReg.set([tmp.r, tmp.g, tmp.b], i * 3);
    tmp.set(KIND_COLORS[n.t] || '#adb5bd'); colKind.set([tmp.r, tmp.g, tmp.b], i * 3);
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
    const hl = S.highlight, cols = S.colorby === 'kind' ? colKind : colReg;
    let shown = 0;
    for (let i = 0; i < N; i++) {
      const v = vis(i), s = S.selected.has(i);
      aSize[i] = v ? baseSize[i] * S.sizeMul * (s ? 1.7 : 1) : 0;
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
    labelNodes = [...S.selected, ...rest];
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
  const DIRSYM = ['⇢', '→', '⊣', '±'];
  let details = null;
  function renderDetail(list) {
    const el = $('detail');
    if (!S.selected.size) { el.hidden = true; return; }
    el.hidden = false;
    const sel = [...S.selected];
    let h = '';
    if (sel.length === 1) {
      const n = nodes[sel[0]];
      h += `<h2>${esc(n.l)}</h2><div class="sub">${esc(n.t)} · ${esc(n.r)}${n.u ? ' · UniProt ' + esc(n.u) : ''}${n.c ? ' · ' + esc(n.c) : ''}</div>`;
      h += `<div class="sub">${n.d} edges in total (${n.dr} regulation)${n.m ? ' · ' + n.m + ' members' : ''}</div>`;
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
  }
  $('detail').addEventListener('click', async (ev) => {
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
    more.innerHTML = `<div>evidence: stimulating ${ns[0]} · inhibiting ${ns[1]} · unsigned ${ns[2]}</div>`
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
    if (!vis(i)) S.forced.add(i);
    S.selected = new Set([i]);
    update();
    if (fly) flyTo(new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]), 120);
  }
  function clearSel() { S.selected = new Set(); S.forced.clear(); update(); }

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
  $('colorby').addEventListener('change', (ev) => { S.colorby = ev.target.value; refreshNodes(); });
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
      b.addEventListener('click', () => { selectOnly(i, true); });
      box.append(b);
    });
  });
  $('q').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') $('results').querySelector('button')?.click(); });

  // ---------------------------------------------------------------- render loop
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
