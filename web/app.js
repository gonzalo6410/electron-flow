// Electron Flow viewer: one Lewis structure drawn in 3D. Renders what the core computes; nothing
// chemical is decided here. Ledger items 13, 14, 24, 34.
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { CASES, runCase } from './core/cases.js';
import { COMPOUNDS } from './core/compounds.js';
import { nodeFromSmiles } from './core/smiles.js';
import { buildNode, enumerate } from './core/lewis.js';
import { relax } from './core/geometry.js';
import { piFrame } from './core/pi.js';
import { structureFamilies } from './core/symmetry.js';
import { ELEMENTS } from './core/elements.js';

// ---------- timing ----------
const DWELL = 4000;   // ms a structure stays before the next equivalent one
const FADE = 1200;    // ms for the old picture to go and the new one to come, one after the other

// ---------- formatting ----------
const MINUS = '−';
// a formal charge as written on a structure: 2+, -, -1/3
function chargeMark(q) {
  if (Math.abs(q) < 1e-9) return '';
  const sign = q < 0 ? MINUS : '+';
  const a = Math.abs(q);
  if (Math.abs(a - 1) < 1e-9) return sign;
  if (Number.isInteger(a)) return `${a}${sign}`;
  for (let d = 2; d <= 12; d++) {
    const num = Math.round(a * d);
    if (Math.abs(a * d - num) < 1e-6) return `${sign}${num}/${d}`;
  }
  return sign + a.toFixed(2);
}

// ---------- scene ----------
const view = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff);
const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 100);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xffffff, 0x999999, 1.0));
const sun = new THREE.DirectionalLight(0xffffff, 1.1);
sun.position.set(3, 5, 6);
scene.add(sun);

function resize() {
  const w = view.clientWidth, h = view.clientHeight;
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);
resize();

// ---------- shapes ----------
const sphereGeo = new THREE.SphereGeometry(1, 40, 28);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 18, 1);
const atomRadius = el => 0.16 + 0.15 * ELEMENTS[el].radius;
const isLight = hex => (0.299 * ((hex >> 16) & 255) + 0.587 * ((hex >> 8) & 255) + 0.114 * (hex & 255)) / 255 > 0.6;
const solid = color => new THREE.MeshLambertMaterial({ color, transparent: true, opacity: 1 });

function cylinder(a, b, radius, material) {
  const d = b.clone().sub(a);
  const m = new THREE.Mesh(cylGeo, material);
  m.position.copy(a).add(b).multiplyScalar(0.5);
  m.scale.set(radius, d.length(), radius);
  m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.normalize());
  return m;
}

function dashed(a, b, radius, material, group) {
  const d = b.clone().sub(a);
  const L = d.length();
  const seg = 0.15, gap = 0.11;
  const n = Math.max(1, Math.floor((L + gap) / (seg + gap)));
  const start = (L - (n * seg + (n - 1) * gap)) / 2;
  const u = d.normalize();
  for (let k = 0; k < n; k++) {
    const s = start + k * (seg + gap);
    group.add(cylinder(a.clone().addScaledVector(u, s), a.clone().addScaledVector(u, s + seg), radius, material));
  }
}

// ---------- sprites (always face the camera) ----------
function canvasSprite(draw, width, height, opacity = 1, depthTest = true) {
  const c = document.createElement('canvas');
  c.width = 256; c.height = 128;
  draw(c.getContext('2d'));
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, opacity, depthTest, depthWrite: false }));
  sp.scale.set(width, height, 1);
  sp.userData.baseOpacity = opacity;
  return sp;
}

// element symbol
function symbolSprite(text, width, color) {
  return canvasSprite(g => {
    g.font = '500 72px "Helvetica Neue", Helvetica, Arial, sans-serif';
    g.textAlign = 'center'; g.textBaseline = 'middle';
    g.fillStyle = color;
    g.fillText(text, 128, 66);
  }, width, width / 2);
}

// formal charge in a small white badge, legible over anything
function chargeSprite(text) {
  const sp = canvasSprite(g => {
    g.font = '600 60px "Helvetica Neue", Helvetica, Arial, sans-serif';
    g.textAlign = 'center'; g.textBaseline = 'middle';
    const w = Math.max(56, g.measureText(text).width + 30), h = 64, r = 32;
    const x0 = 128 - w / 2, y0 = 64 - h / 2;
    g.beginPath();
    g.moveTo(x0 + r, y0); g.lineTo(x0 + w - r, y0); g.arc(x0 + w - r, y0 + r, r, -Math.PI / 2, Math.PI / 2);
    g.lineTo(x0 + r, y0 + h); g.arc(x0 + r, y0 + r, r, Math.PI / 2, 3 * Math.PI / 2); g.closePath();
    g.fillStyle = '#fff'; g.fill();
    g.lineWidth = 5; g.strokeStyle = '#111'; g.stroke();
    g.fillStyle = '#111';
    g.fillText(text, 128, 66);
  }, 0.56, 0.28, 1, false);
  sp.renderOrder = 10;
  return sp;
}

// a lone pair: two dots
function pairSprite(opacity) {
  return canvasSprite(g => {
    g.fillStyle = '#111';
    for (const x of [96, 160]) { g.beginPath(); g.arc(x, 64, 18, 0, Math.PI * 2); g.fill(); }
  }, 0.30, 0.15, opacity);
}

// ---------- camera framing ----------
// principal axes of a point cloud (largest, middle, smallest spread), by Jacobi rotations
function principalAxes(points) {
  const n = Math.max(1, points.length);
  const c = [0, 1, 2].map(a => points.reduce((s, p) => s + p[a], 0) / n);
  const m = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of points) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) m[a][b] += (p[a] - c[a]) * (p[b] - c[b]);
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 30; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) off += m[p][q] * m[p][q];
    if (off < 1e-12) break;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(m[p][q]) < 1e-12) continue;
      const th = 0.5 * Math.atan2(2 * m[p][q], m[q][q] - m[p][p]);
      const cs = Math.cos(th), sn = Math.sin(th);
      for (let k = 0; k < 3; k++) { const a = m[k][p], b = m[k][q]; m[k][p] = cs * a - sn * b; m[k][q] = sn * a + cs * b; }
      for (let k = 0; k < 3; k++) { const a = m[p][k], b = m[q][k]; m[p][k] = cs * a - sn * b; m[q][k] = sn * a + cs * b; }
      for (let k = 0; k < 3; k++) { const a = v[k][p], b = v[k][q]; v[k][p] = cs * a - sn * b; v[k][q] = sn * a + cs * b; }
    }
  }
  const order = [0, 1, 2].sort((a, b) => m[b][b] - m[a][a]);
  const col = k => new THREE.Vector3(v[0][k], v[1][k], v[2][k]).normalize();
  const major = col(order[0]), middle = col(order[1]);
  const minor = major.clone().cross(middle).normalize();
  if (minor.z < 0) minor.negate();
  const flat = m[order[2]][order[2]] < 0.5 * m[order[0]][order[0]];
  return { major, middle, minor, flat };
}

function frameCamera(node, geo) {
  const heavy = geo.pos.filter((p, i) => node.atoms[i].el !== 'H');
  const extent = Math.max(1.5, ...geo.pos.map(p => Math.hypot(...p)));
  const dist = extent / Math.tan((camera.fov / 2) * Math.PI / 180) * 1.5 + 1;
  const axes = principalAxes(heavy.length ? heavy : geo.pos);
  let dir;
  if (axes.flat) {
    // a flat molecule: look at its face, tilted a little so depth still shows
    dir = axes.minor.clone().multiplyScalar(Math.cos(0.35)).addScaledVector(axes.major, Math.sin(0.35) * 0.6).addScaledVector(axes.middle, Math.sin(0.35) * 0.8).normalize();
    camera.up.copy(axes.middle);
  } else {
    dir = new THREE.Vector3(0.38, 0.42, 0.82).normalize(); // a round molecule: off every axis
    camera.up.set(0, 1, 0);
  }
  camera.position.copy(dir.multiplyScalar(dist));
  controls.target.set(0, 0, 0);
  controls.update();
}

// ---------- state ----------
let current = null;          // { node, r, geo, families, main }
let showResonance = false;
let familyIdx = 0, memberIdx = 0;
let lastSwitch = 0;
let fade = null;             // { from, to, start }
const atomGroup = new THREE.Group();      // spheres and symbols, fixed per molecule
let structureGroup = new THREE.Group();   // bonds, dots, charges, per structure
scene.add(atomGroup, structureGroup);
const tracked = [];          // sprites placed relative to the camera every frame

function load(entry) {
  let node, r;
  if (entry.kind === 'case') ({ node, r } = runCase(entry.c));
  else { node = buildNode(nodeFromSmiles(entry.smiles)); r = enumerate(node); }
  const geo = relax(node, r);
  const { families } = structureFamilies(node, r);
  // the main structure: least total formal charge, first one on a tie (all weights are equal; this
  // only picks what to show when resonance is switched off)
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  let main = 0;
  r.structures.forEach((st, q) => { if (spread(st) < spread(r.structures[main])) main = q; });
  current = { node, r, geo, families, main };
  familyIdx = 0; memberIdx = 0; fade = null;
  lastSwitch = performance.now();
  frameCamera(node, geo);
  buildAtoms();
  swapStructure(buildStructure(), false);
  updatePanel();
}

function shownStructure() {
  const { r, families, main } = current;
  return showResonance ? r.structures[families[familyIdx][memberIdx]] : r.structures[main];
}

function buildAtoms() {
  atomGroup.clear();
  tracked.length = 0;
  const { node, geo } = current;
  node.atoms.forEach((a, i) => {
    const rad = atomRadius(a.el);
    const m = new THREE.Mesh(sphereGeo, solid(ELEMENTS[a.el].color));
    m.position.set(...geo.pos[i]);
    m.scale.setScalar(rad);
    atomGroup.add(m);
    const sp = symbolSprite(a.el, Math.max(0.4, 1.5 * rad), isLight(ELEMENTS[a.el].color) ? '#111' : '#fff');
    atomGroup.add(sp);
    tracked.push({ sprite: sp, index: i, kind: 'symbol' });
  });
}

// the Lewis drawing of one structure: lines for bonds (a dashed one for a fraction), dots for lone
// pairs, badges for charges, a grey ghost where an adjacent pair carries no bond in this structure
function buildStructure() {
  const g = new THREE.Group();
  const { node, r, geo } = current;
  const st = shownStructure();
  const pos = geo.pos.map(p => new THREE.Vector3(...p));
  const sym = node.atoms.map(a => a.el);
  const frame = piFrame(node, r, geo.pos, geo.groups, st);
  const bondMat = solid(0x3a3a3a), ghostMat = solid(0xd0d0d0), ionMat = solid(0xc8c8c8);

  r.edges.forEach(([i, j], k) => {
    const o = st.orders[k];
    const a = pos[i], b = pos[j];
    const bd = b.clone().sub(a).normalize();
    const a2 = a.clone().addScaledVector(bd, atomRadius(sym[i]) * 0.9);
    const b2 = b.clone().addScaledVector(bd, -atomRadius(sym[j]) * 0.9);
    if (o <= 1e-9) {
      if (node.fixed[i] || node.fixed[j]) dashed(a2, b2, 0.03, ionMat, g);
      else g.add(cylinder(a2, b2, 0.04, ghostMat));
      return;
    }
    // lines side by side in the plane of the molecule (across the p orbitals), extra lines toward
    // the neighbours so they sit inside a ring
    const full = Math.floor(o + 1e-9), fr = o - full;
    const lines = full + (fr > 1e-9 ? 1 : 0);
    let side;
    if (frame.normals[k].length) side = bd.clone().cross(new THREE.Vector3(...frame.normals[k][0])).normalize();
    else {
      const ref = Math.abs(bd.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
      side = ref.addScaledVector(bd, -ref.dot(bd)).normalize();
    }
    if (lines > 1) {
      const mid = a.clone().add(b).multiplyScalar(0.5);
      const toward = new THREE.Vector3();
      for (const [p, q] of r.edges) {
        const other = p === i ? q : q === i ? p : p === j ? q : q === j ? p : -1;
        if (other >= 0 && other !== i && other !== j) toward.add(pos[other].clone().sub(mid));
      }
      if (toward.dot(side) < 0) side.negate();
    }
    for (let q = 0; q < lines; q++) {
      const shift = side.clone().multiplyScalar(0.17 * (q - (lines - 1) / 2));
      const pa = a2.clone().add(shift), pb = b2.clone().add(shift);
      if (q === lines - 1 && fr > 1e-9) dashed(pa, pb, 0.055, bondMat, g);
      else g.add(cylinder(pa, pb, 0.07, bondMat));
    }
  });

  // lone pairs: laid out around the atom in the plane of the screen, every frame (placeSprites)
  sym.forEach((el, i) => {
    if (node.fixed[i]) return;
    const pairs = st.lone[i] / 2;
    const whole = Math.floor(pairs + 1e-9);
    const partners = frame.lonePairs[i].filter(lp => lp.kind === 'partner').map(lp => lp.partner);
    const count = frame.lonePairs[i].length;
    for (let q = 0; q < count; q++) {
      const sp = pairSprite(q < whole ? 1 : Math.max(0.2, pairs - whole));
      g.add(sp);
      tracked.push({ sprite: sp, index: i, kind: 'pair', q, count, partners, group: g });
    }
  });

  // charges
  sym.forEach((el, i) => {
    if (Math.abs(st.charge[i]) < 1e-9) return;
    const sp = chargeSprite(chargeMark(st.charge[i]));
    g.add(sp);
    tracked.push({ sprite: sp, index: i, kind: 'charge', group: g });
  });
  return g;
}

function setOpacity(group, f) {
  group.traverse(o => { if (o.material) o.material.opacity = f * (o.userData.baseOpacity ?? 1); });
}

function dropTracked(group) {
  for (let k = tracked.length - 1; k >= 0; k--) if (tracked[k].group === group) tracked.splice(k, 1);
}

function swapStructure(next, animate) {
  if (fade) { // a fade still running is finished on the spot
    scene.remove(fade.from); dropTracked(fade.from);
    structureGroup = fade.to; setOpacity(structureGroup, 1); fade = null;
  }
  if (!animate) {
    scene.remove(structureGroup); dropTracked(structureGroup);
    structureGroup = next; scene.add(next);
    return;
  }
  setOpacity(next, 0);
  scene.add(next);
  fade = { from: structureGroup, to: next, start: performance.now() };
}

function showStructure(animate) {
  swapStructure(buildStructure(), animate);
  updatePanel();
}

// ---------- per-frame sprite placement ----------
const camRight = new THREE.Vector3(), camUp = new THREE.Vector3(), camBack = new THREE.Vector3();
const angleCache = new Map();

// screen-plane angles for an atom's dot pairs: pairs that belong to a missing or half bond point at
// their partner, the rest go into the widest gaps between the bonds as seen from the camera
function pairAngles(i, count, partners) {
  const key = `${i}:${count}:${partners.join(',')}`;
  if (angleCache.has(key)) return angleCache.get(key);
  const pos = current.geo.pos;
  const p = new THREE.Vector3(...pos[i]);
  const angleOf = j => { const v = new THREE.Vector3(...pos[j]).sub(p); const x = v.dot(camRight), y = v.dot(camUp); return Math.hypot(x, y) > 0.25 ? Math.atan2(y, x) : null; };
  const taken = [];
  for (const nb of current.geo.groups[i].nbrs) { const a = angleOf(nb.j); if (a !== null) taken.push(a); }
  for (const j of current.geo.groups[i].assoc) { const a = angleOf(j); if (a !== null) taken.push(a); }
  const out = partners.map(j => angleOf(j) ?? 0);
  while (out.length < count) {
    const all = [...taken, ...out].sort((a, b) => a - b);
    if (!all.length) { out.push(Math.PI / 2); continue; }
    let best = { mid: 0, size: -1 };
    for (let q = 0; q < all.length; q++) {
      const a = all[q], b = q + 1 < all.length ? all[q + 1] : all[0] + 2 * Math.PI;
      if (b - a > best.size) best = { mid: (a + b) / 2, size: b - a };
    }
    out.push(best.mid);
  }
  angleCache.set(key, out);
  return out;
}

function placeSprites() {
  camera.matrixWorld.extractBasis(camRight, camUp, camBack);
  angleCache.clear();
  const pos = current.geo.pos;
  const toCam = new THREE.Vector3();
  for (const t of tracked) {
    const p = new THREE.Vector3(...pos[t.index]);
    const rad = atomRadius(current.node.atoms[t.index].el);
    toCam.copy(camera.position).sub(p).normalize();
    if (t.kind === 'symbol') {
      t.sprite.position.copy(p).addScaledVector(toCam, rad * 1.25);
    } else if (t.kind === 'pair') {
      const a = pairAngles(t.index, t.count, t.partners)[t.q];
      t.sprite.position.copy(p).addScaledVector(toCam, rad * 1.1).addScaledVector(camRight, Math.cos(a) * (rad + 0.2)).addScaledVector(camUp, Math.sin(a) * (rad + 0.2));
    } else {
      t.sprite.position.copy(p).addScaledVector(toCam, rad).addScaledVector(camRight, rad * 0.95 + 0.08).addScaledVector(camUp, rad * 0.85 + 0.06);
    }
  }
}

// ---------- panel ----------
const sel = document.getElementById('case');
const ENTRIES = [
  ...COMPOUNDS.map(c => ({ kind: 'smiles', name: c.name, smiles: c.smiles, group: 'compounds' })),
  ...CASES.filter(c => c.before).map(c => ({ kind: 'case', name: c.name, c, group: 'transitions' })),
];
for (const groupName of ['compounds', 'transitions']) {
  const og = document.createElement('optgroup');
  og.label = groupName;
  ENTRIES.forEach((e, i) => { if (e.group !== groupName) return; const o = document.createElement('option'); o.value = i; o.textContent = e.name; og.appendChild(o); });
  sel.appendChild(og);
}
const typedOption = document.createElement('option');
typedOption.value = 'typed'; typedOption.hidden = true;
sel.prepend(typedOption);
const errBox = document.getElementById('error');
function tryLoad(entry) {
  try { load(entry); errBox.textContent = ''; }
  catch (e) { errBox.textContent = e.message; }
}
sel.addEventListener('change', () => { if (sel.value !== 'typed') tryLoad(ENTRIES[Number(sel.value)]); });
const smilesBox = document.getElementById('smiles');
document.getElementById('smiles-form').addEventListener('submit', e => {
  e.preventDefault();
  const text = smilesBox.value.trim();
  if (!text) return;
  tryLoad({ kind: 'smiles', name: text, smiles: text });
  typedOption.textContent = `typed: ${text}`;
  typedOption.hidden = false;
  sel.value = 'typed';
  smilesBox.blur();
});

const resonanceBox = document.getElementById('resonance');
resonanceBox.addEventListener('change', () => {
  showResonance = resonanceBox.checked;
  familyIdx = 0; memberIdx = 0;
  lastSwitch = performance.now();
  showStructure(true);
});

function stepFamily(d) {
  if (!current || !showResonance) return;
  const nf = current.families.length;
  if (nf < 2) return;
  familyIdx = (familyIdx + d + nf) % nf;
  memberIdx = 0;
  lastSwitch = performance.now();
  showStructure(true);
}
document.getElementById('prev').addEventListener('click', () => stepFamily(-1));
document.getElementById('next').addEventListener('click', () => stepFamily(1));
window.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
  if (e.key === 'ArrowLeft') stepFamily(-1);
  if (e.key === 'ArrowRight') stepFamily(1);
});

function updatePanel() {
  const { r, families } = current;
  const n = r.structures.length;
  document.getElementById('steps').style.visibility = (showResonance && families.length > 1) ? 'visible' : 'hidden';
  let line = `<b>${r.pool}</b> valence electrons · <b>${n}</b> contributing structure${n === 1 ? '' : 's'}`;
  if (r.transitionEdges) line += '<br>transition: before and after each weighted 1/2';
  if (showResonance && n > 1) {
    const size = families[familyIdx].length;
    if (families.length > 1) line += `<br>structure <b>${familyIdx + 1}</b> of ${families.length}`;
    if (size > 1) line += `<br><b>${size}</b> equivalent by symmetry, alternating`;
  } else if (n > 1) line += '<br>main structure shown';
  document.getElementById('status').innerHTML = line;
}

// ---------- loop ----------
function animate(t) {
  requestAnimationFrame(animate);
  if (current) {
    if (fade) {
      const f = Math.min(1, (t - fade.start) / FADE);
      const outPart = Math.min(1, f * 2), inPart = Math.max(0, f * 2 - 1); // old goes, then new comes
      const ease = x => x * x * (3 - 2 * x);
      setOpacity(fade.from, 1 - ease(outPart));
      setOpacity(fade.to, ease(inPart));
      if (f >= 1) { scene.remove(fade.from); dropTracked(fade.from); structureGroup = fade.to; fade = null; }
    }
    const family = current.families[familyIdx];
    if (showResonance && family.length > 1 && !fade && t - lastSwitch > DWELL) {
      memberIdx = (memberIdx + 1) % family.length;
      lastSwitch = t;
      showStructure(true);
    }
    placeSprites();
  }
  controls.update();
  renderer.render(scene, camera);
}

const first = ENTRIES.findIndex(e => e.name === 'SF6');
sel.value = first;
tryLoad(ENTRIES[first]);
requestAnimationFrame(animate);
