// Electron Flow viewer. Renders what the core computes; nothing chemical is decided here.
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { CASES, runCase } from './core/cases.js';
import { COMPOUNDS } from './core/compounds.js';
import { nodeFromSmiles } from './core/smiles.js';
import { buildNode, enumerate } from './core/lewis.js';
import { relax } from './core/geometry.js';
import { piFrame, repulsion } from './core/pi.js';
import { ELEMENTS } from './core/elements.js';
import { structureFamilies } from './core/symmetry.js';

// ---------- formatting ----------
const MINUS = '−';
function frac(x, signed = true) {
  if (Math.abs(x) < 1e-9) return '0';
  const sign = x < 0 ? MINUS : (signed ? '+' : '');
  const a = Math.abs(x);
  for (let d = 1; d <= 12; d++) {
    const nmr = Math.round(a * d);
    if (Math.abs(a * d - nmr) < 1e-6) return sign + (d === 1 ? String(nmr) : `${nmr}/${d}`);
  }
  return sign + a.toFixed(2);
}
// charge as chemists write it on a structure: 2+, −, −1/3
function chargeMark(q) {
  if (Math.abs(q) < 1e-9) return '';
  const s = q < 0 ? MINUS : '+';
  const a = Math.abs(q);
  if (Math.abs(a - 1) < 1e-9) return s;
  if (Number.isInteger(a)) return `${a}${s}`;
  return frac(q);
}

// ---------- scene ----------
const view = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff);
const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 100);
camera.position.set(0, 2, 9);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xffffff, 0x999999, 1.0));
const dir = new THREE.DirectionalLight(0xffffff, 1.1);
dir.position.set(3, 5, 6);
scene.add(dir);

function resize() {
  const w = view.clientWidth, h = view.clientHeight;
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);
resize();

const sphereGeo = new THREE.SphereGeometry(1, 40, 28);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 18, 1);
const atomRadius = el => 0.16 + 0.15 * ELEMENTS[el].radius;
const isLight = hex => (0.299 * ((hex >> 16) & 255) + 0.587 * ((hex >> 8) & 255) + 0.114 * (hex & 255)) / 255 > 0.6;

function solid(color) { return new THREE.MeshLambertMaterial({ color, transparent: true, opacity: 1 }); }

function cylinder(a, b, radius, material) {
  const va = new THREE.Vector3(...a), vb = new THREE.Vector3(...b);
  const d = vb.clone().sub(va);
  const m = new THREE.Mesh(cylGeo, material);
  m.position.copy(va).add(vb).multiplyScalar(0.5);
  m.scale.set(radius, d.length(), radius);
  m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.clone().normalize());
  return m;
}

function dashed(a, b, radius, material, group) {
  const va = new THREE.Vector3(...a), vb = new THREE.Vector3(...b);
  const d = vb.clone().sub(va);
  const L = d.length();
  const seg = 0.15, gap = 0.11;
  const n = Math.max(1, Math.floor((L + gap) / (seg + gap)));
  const total = n * seg + (n - 1) * gap;
  const start = (L - total) / 2;
  const u = d.clone().normalize();
  for (let k = 0; k < n; k++) {
    const s = start + k * (seg + gap);
    group.add(cylinder(va.clone().add(u.clone().multiplyScalar(s)).toArray(), va.clone().add(u.clone().multiplyScalar(s + seg)).toArray(), radius, material));
  }
}

// text sprite; width in scene units, canvas aspect 2:1
function textSprite(text, width, color = '#111', weight = '500', opacity = 1) {
  const c = document.createElement('canvas');
  c.width = 256; c.height = 128;
  const g = c.getContext('2d');
  g.font = `${weight} 72px "Helvetica Neue", Helvetica, Arial, sans-serif`;
  g.textAlign = 'center'; g.textBaseline = 'middle';
  g.fillStyle = color;
  g.fillText(text, 128, 66);
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, opacity, depthTest: true, depthWrite: false }));
  sp.scale.set(width, width / 2, 1);
  return sp;
}

// a formal charge in a small white badge, so it stays legible over bonds and atoms
function chargeSprite(text, opacity = 1) {
  const c = document.createElement('canvas');
  c.width = 256; c.height = 128;
  const g = c.getContext('2d');
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
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, opacity, depthTest: false, depthWrite: false }));
  sp.scale.set(0.56, 0.28, 1);
  sp.renderOrder = 10;
  return sp;
}

// a lone pair: two dots side by side, always facing the camera
function pairSprite(opacity = 1) {
  const c = document.createElement('canvas');
  c.width = 128; c.height = 64;
  const g = c.getContext('2d');
  g.fillStyle = '#111';
  for (const x of [40, 88]) { g.beginPath(); g.arc(x, 32, 15, 0, Math.PI * 2); g.fill(); }
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, opacity, depthTest: true, depthWrite: false }));
  sp.scale.set(0.30, 0.15, 1);
  return sp;
}

// ---------- state ----------
let current = null;   // {c, node, r, geo, lpDirs}
let mode = 'cycle';
let drawStyle = 'lewis'; // 'lewis' or 'orbital'
let familyIdx = 0;   // which family of equivalent structures is shown (stepped by hand)
let memberIdx = 0;   // which member of the family (cycles on its own)
let lastSwitch = 0;
let fade = null;      // {from, to, start}
const atomGroup = new THREE.Group();   // spheres and symbols, fixed per case
let structureGroup = new THREE.Group(); // bonds, charges, lone pairs, per structure
scene.add(atomGroup);
scene.add(structureGroup);
const faceCamera = []; // {sprite, index, kind} sprites repositioned every frame

function load(entry) {
  let node, r, c;
  if (entry.kind === 'case') {
    c = entry.c;
    ({ node, r } = runCase(c));
  } else {
    c = { name: entry.name };
    node = buildNode(nodeFromSmiles(entry.smiles));
    r = enumerate(node);
  }
  const geo = relax(node, r);
  const fam = structureFamilies(node, r);
  current = { c, node, r, geo, families: fam.families };
  familyIdx = 0; memberIdx = 0;
  lastSwitch = performance.now();
  fade = null;
  document.getElementById('rules').textContent =
    'Rules: the octet is never exceeded; contributing structures are the ones with the fewest octet holes; sigma bonds break only where the octet forces it; equal weights; geometry from electron groups and springs.';
  // frame the molecule: look at its flattest face, tilted a little so depth still shows
  const extent = Math.max(1.5, ...geo.pos.map(p => Math.hypot(...p)));
  const dist = extent / Math.tan((camera.fov / 2) * Math.PI / 180) * 1.5 + 1;
  const axes = principalAxes(geo.pos.filter((p, i) => node.atoms[i].el !== 'H'));
  const view = axes.minor.clone().multiplyScalar(Math.cos(0.35)).add(axes.major.clone().multiplyScalar(Math.sin(0.35) * 0.6)).add(axes.middle.clone().multiplyScalar(Math.sin(0.35) * 0.8)).normalize();
  camera.position.copy(view.multiplyScalar(dist));
  camera.up.copy(axes.middle);
  controls.target.set(0, 0, 0);
  controls.update();
  buildAtoms();
  swapStructure(buildStructure(picture()), false);
  document.getElementById('steps').style.visibility = (mode === 'cycle' && fam.families.length > 1) ? 'visible' : 'hidden';
  readout();
}

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
      for (let k = 0; k < 3; k++) {
        const a = m[k][p], b = m[k][q];
        m[k][p] = cs * a - sn * b; m[k][q] = sn * a + cs * b;
      }
      for (let k = 0; k < 3; k++) {
        const a = m[p][k], b = m[q][k];
        m[p][k] = cs * a - sn * b; m[q][k] = sn * a + cs * b;
      }
      for (let k = 0; k < 3; k++) {
        const a = v[k][p], b = v[k][q];
        v[k][p] = cs * a - sn * b; v[k][q] = sn * a + cs * b;
      }
    }
  }
  const order = [0, 1, 2].sort((a, b) => m[b][b] - m[a][a]);
  const col = k => new THREE.Vector3(v[0][k], v[1][k], v[2][k]).normalize();
  const major = col(order[0]), middle = col(order[1]);
  const minor = major.clone().cross(middle).normalize();
  if (minor.z < 0) minor.negate();
  return { major, middle, minor };
}

function picture() {
  const { r } = current;
  if (mode === 'cycle') {
    const st = r.structures[current.families[familyIdx][memberIdx]];
    return { orders: st.orders, charge: st.charge, lone: st.lone, count: st.count, single: true };
  }
  return { orders: r.avg.orders, charge: r.avg.charge, lone: r.avg.lone, count: r.avg.count, single: false };
}

function buildAtoms() {
  atomGroup.clear();
  faceCamera.length = 0;
  const { node, geo } = current;
  node.atoms.forEach((a, i) => {
    const rad = atomRadius(a.el);
    const m = new THREE.Mesh(sphereGeo, solid(ELEMENTS[a.el].color));
    m.position.set(...geo.pos[i]);
    m.scale.setScalar(rad);
    atomGroup.add(m);
    const sp = textSprite(a.el, Math.max(0.4, 1.5 * rad), isLight(ELEMENTS[a.el].color) ? '#111' : '#fff', '500');
    atomGroup.add(sp);
    faceCamera.push({ sprite: sp, index: i, kind: 'symbol' });
  });
}

function buildStructure(pic) {
  const style = drawStyle;
  const g = new THREE.Group();
  const { node, r, geo } = current;
  const pos = geo.pos;
  const sym = node.atoms.map(a => a.el);
  const frame = piFrame(node, r, pos, geo.groups, pic);
  const bondMat = solid(0x3a3a3a);
  const ghostMat = solid(0xd0d0d0);
  const ionMat = solid(0xc8c8c8);

  // bonds
  r.edges.forEach(([i, j], k) => {
    const o = pic.orders[k];
    const loose = node.fixed[i] || node.fixed[j];
    const a = new THREE.Vector3(...pos[i]), b = new THREE.Vector3(...pos[j]);
    const bd = b.clone().sub(a).normalize();
    const ra = atomRadius(sym[i]) * 0.9, rb = atomRadius(sym[j]) * 0.9;
    const a2 = a.clone().add(bd.clone().multiplyScalar(ra)), b2 = b.clone().sub(bd.clone().multiplyScalar(rb));
    if (o <= 1e-9) {
      // a pair that is adjacent but carries no bond in this picture: ion pair dashed, missing bond ghosted
      if (loose) dashed(a2.toArray(), b2.toArray(), 0.03, ionMat, g);
      else g.add(cylinder(a2.toArray(), b2.toArray(), 0.04, ghostMat));
      return;
    }
    if (style === 'orbital') {
      // sigma part a rod, solid if whole, dashed if partial
      if (o >= 1 - 1e-9) g.add(cylinder(a2.toArray(), b2.toArray(), 0.07, bondMat));
      else dashed(a2.toArray(), b2.toArray(), 0.06, bondMat, g);
      // each pi part a translucent lobe above and below, along the p orbitals; fainter when partial
      frame.normals[k].forEach((nrm, q) => {
        const amount = Math.min(1, o - 1 - q);
        if (amount <= 1e-9) return;
        const mid = a2.clone().add(b2).multiplyScalar(0.5);
        const L = a2.distanceTo(b2);
        const nv = new THREE.Vector3(...nrm);
        const mat = new THREE.MeshLambertMaterial({ color: 0x5577aa, transparent: true, opacity: 0.22 + 0.3 * amount, depthWrite: false });
        for (const side of [1, -1]) {
          const lobe = new THREE.Mesh(sphereGeo, mat);
          lobe.userData.baseOpacity = 0.22 + 0.3 * amount;
          lobe.renderOrder = 1;
          lobe.position.copy(mid).add(nv.clone().multiplyScalar(0.24 * side));
          lobe.scale.set(0.17, L * 0.46, 0.15);
          const m = new THREE.Matrix4().makeBasis(bd.clone().cross(nv).normalize(), bd.clone(), nv.clone()); // right-handed: x = y cross z
          lobe.quaternion.setFromRotationMatrix(m);
          g.add(lobe);
        }
      });
      return;
    }
    // Lewis style: one line per whole order plus a dashed line for a fraction, side by side in the
    // plane of the molecule (the in-plane perpendicular to the p orbitals), extra lines toward the
    // neighbours (inside a ring)
    const full = Math.floor(o + 1e-9), fr = o - full;
    const lines = full + (fr > 1e-9 ? 1 : 0);
    let side;
    if (frame.normals[k].length) side = bd.clone().cross(new THREE.Vector3(...frame.normals[k][0])).normalize();
    else {
      const ref = Math.abs(bd.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
      side = ref.sub(bd.clone().multiplyScalar(ref.dot(bd))).normalize();
    }
    if (lines > 1) {
      const mid = a.clone().add(b).multiplyScalar(0.5);
      let toward = new THREE.Vector3();
      for (const [p, q] of r.edges) {
        const other = p === i ? q : q === i ? p : p === j ? q : q === j ? p : -1;
        if (other < 0 || other === i || other === j) continue;
        toward.add(new THREE.Vector3(...pos[other]).sub(mid));
      }
      if (toward.dot(side) < 0) side.negate();
    }
    const off = side.multiplyScalar(0.17);
    for (let q = 0; q < lines; q++) {
      const shift = off.clone().multiplyScalar(q - (lines - 1) / 2);
      const pa = a2.clone().add(shift).toArray(), pb = b2.clone().add(shift).toArray();
      const partial = q === lines - 1 && fr > 1e-9;
      if (partial) dashed(pa, pb, 0.055, bondMat, g);
      else g.add(cylinder(pa, pb, 0.07, bondMat));
    }
  });

  // lone pairs: whole pairs full, a fractional pair faded. Orbital style places them in their
  // orbitals; Lewis style keeps the in-plane and partner placements but spreads the p-orbital pairs
  // around the atom by repulsion, as dots on a drawing would be
  sym.forEach((el, i) => {
    if (node.fixed[i]) return;
    const rad = atomRadius(el);
    const pairs = pic.lone[i] / 2;
    const whole = Math.floor(pairs + 1e-9);
    let lps = frame.lonePairs[i];
    if (style === 'lewis' && lps.some(lp => lp.kind === 'p' || lp.kind === 'free')) {
      const kept = lps.filter(lp => lp.kind !== 'p' && lp.kind !== 'free');
      const need = lps.length - kept.length;
      const bonds = geo.groups[i].nbrs.map(nb => new THREE.Vector3(...pos[nb.j]).sub(new THREE.Vector3(...pos[i])).normalize());
      let extra = [];
      const pO = frame.pOrbital[i];
      if (pO && bonds.length >= 2) {
        // in the plane, in the widest gaps between the bonds, as dots on a drawing
        const nrm = new THREE.Vector3(...pO);
        const ex = bonds[0].clone().sub(nrm.clone().multiplyScalar(bonds[0].dot(nrm))).normalize();
        const ey = nrm.clone().cross(ex);
        const angles = [...bonds, ...kept.map(lp => new THREE.Vector3(...lp.dir))].map(v => Math.atan2(v.dot(ey), v.dot(ex))).sort((a, b) => a - b);
        const gaps = angles.map((a, q) => { const b = q + 1 < angles.length ? angles[q + 1] : angles[0] + 2 * Math.PI; return { mid: (a + b) / 2, size: b - a }; }).sort((a, b) => b.size - a.size);
        for (let q = 0; q < need && q < gaps.length; q++) {
          const v = ex.clone().multiplyScalar(Math.cos(gaps[q].mid)).add(ey.clone().multiplyScalar(Math.sin(gaps[q].mid)));
          extra.push({ dir: v.toArray(), kind: 'free', far: true });
        }
      }
      if (extra.length < need) extra.push(...repulsion(bonds.map(b => b.toArray()), [...kept, ...extra].map(lp => lp.dir), need - extra.length).map(dir => ({ dir, kind: 'free' })));
      lps = [...kept, ...extra];
    }
    lps.forEach((lp, q) => {
      const opacity = q < whole ? 1 : Math.max(0.2, pairs - whole);
      const sp = pairSprite(opacity);
      sp.userData.baseOpacity = opacity;
      const d = lp.dir;
      const out = rad + (lp.far ? 0.3 : 0.17);
      sp.position.set(pos[i][0] + d[0] * out, pos[i][1] + d[1] * out, pos[i][2] + d[2] * out);
      g.add(sp);
    });
  });

  // charges, placed beside each atom every frame
  sym.forEach((el, i) => {
    const q = pic.charge[i];
    if (Math.abs(q) < 1e-9) return;
    const sp = chargeSprite(chargeMark(q));
    sp.userData.baseOpacity = 1;
    g.add(sp);
    faceCamera.push({ sprite: sp, index: i, kind: 'charge', group: g });
  });
  return g;
}

function setOpacity(group, f) {
  group.traverse(o => {
    if (o.material) o.material.opacity = f * (o.userData.baseOpacity ?? 1);
  });
}

function swapStructure(next, animate = true) {
  // a fade still running is finished on the spot
  if (fade) {
    scene.remove(fade.from);
    for (let k = faceCamera.length - 1; k >= 0; k--) if (faceCamera[k].group === fade.from) faceCamera.splice(k, 1);
    structureGroup = fade.to;
    setOpacity(structureGroup, 1);
    fade = null;
  }
  // drop face-camera entries that belong to no live group
  for (let k = faceCamera.length - 1; k >= 0; k--) if (faceCamera[k].group && faceCamera[k].group !== next && faceCamera[k].group !== structureGroup) faceCamera.splice(k, 1);
  if (!animate) {
    scene.remove(structureGroup);
    for (let k = faceCamera.length - 1; k >= 0; k--) if (faceCamera[k].group === structureGroup) faceCamera.splice(k, 1);
    structureGroup = next;
    scene.add(next);
    fade = null;
    return;
  }
  setOpacity(next, 0);
  scene.add(next);
  fade = { from: structureGroup, to: next, start: performance.now() };
}

const FADE = 600;   // ms crossfade
const DWELL = 2200; // ms per structure, equal weights so equal time

function familyLine() {
  const fams = current.families;
  const size = fams[familyIdx].length;
  const a = fams.length > 1 ? `<br>structure <b>${familyIdx + 1}</b> of ${fams.length} (step with the arrows)` : '';
  const b = size > 1 ? `<br><b>${size}</b> equivalent by symmetry, cycling` : '';
  return a + b;
}

function stepFamily(d) {
  if (!current || mode !== 'cycle') return;
  const nf = current.families.length;
  if (nf < 2) return;
  familyIdx = (familyIdx + d + nf) % nf;
  memberIdx = 0;
  lastSwitch = performance.now();
  swapStructure(buildStructure(picture()), true);
  readout();
}

function readout() {
  const { node, r } = current;
  const pic = picture();
  const sym = node.atoms.map(a => a.el);
  const n = r.structures.length;
  const st = document.getElementById('status');
  st.innerHTML = `<b>${r.pool}</b> valence electrons · <b>${n}</b> contributing structure${n === 1 ? '' : 's'} · <b>${r.holes / 2}</b> octet hole${r.holes === 2 ? '' : 's'}` +
    (r.transitionEdges ? '<br>transition: before and after each weighted 1/2' : '') +
    (mode === 'cycle' && n > 1 ? familyLine() : '') +
    (mode === 'average' ? '<br>weighted average of all contributing structures' : '');
  const label = i => `${sym[i]}<sub>${i}</sub>`;
  document.getElementById('atoms').innerHTML =
    '<tr><th></th><th>charge</th><th>count</th><th>lone e</th></tr>' +
    sym.map((el, i) => node.fixed[i]
      ? `<tr><td class="sym">${label(i)}</td><td>${frac(pic.charge[i])}</td><td class="dim">ion</td><td class="dim"></td></tr>`
      : `<tr><td class="sym">${label(i)}</td><td>${frac(pic.charge[i])}</td><td>${frac(pic.count[i], false)}</td><td>${frac(pic.lone[i], false)}</td></tr>`).join('');
  document.getElementById('bonds').innerHTML =
    '<tr><th></th><th>order</th></tr>' +
    r.edges.map(([i, j], k) => {
      const o = pic.orders[k];
      const loose = node.fixed[i] || node.fixed[j];
      return `<tr><td class="sym">${label(i)} – ${label(j)}</td><td${loose ? ' class="dim"' : ''}>${loose && o <= 1e-9 ? 'ion pair' : frac(o, false)}</td></tr>`;
    }).join('');
}

// ---------- controls ----------
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
const errBox = document.getElementById('error');
function tryLoad(entry) {
  try { load(entry); errBox.textContent = ''; }
  catch (e) { errBox.textContent = e.message; }
}
sel.addEventListener('change', () => tryLoad(ENTRIES[Number(sel.value)]));
const smilesBox = document.getElementById('smiles');
document.getElementById('smiles-form').addEventListener('submit', e => {
  e.preventDefault();
  const text = smilesBox.value.trim();
  if (!text) return;
  tryLoad({ kind: 'smiles', name: text, smiles: text });
  smilesBox.blur();
});
const bCycle = document.getElementById('btn-cycle'), bAvg = document.getElementById('btn-avg');
function setMode(m) {
  mode = m;
  document.getElementById('steps').style.visibility = (m === 'cycle' && current.families.length > 1) ? 'visible' : 'hidden';
  bCycle.classList.toggle('on', m === 'cycle');
  bAvg.classList.toggle('on', m === 'average');
  lastSwitch = performance.now();
  swapStructure(buildStructure(picture()), true);
  readout();
}
bCycle.addEventListener('click', () => setMode('cycle'));
bAvg.addEventListener('click', () => setMode('average'));
const bLewis = document.getElementById('btn-lewis'), bOrb = document.getElementById('btn-orbital');
function setStyle(st) {
  drawStyle = st;
  bLewis.classList.toggle('on', st === 'lewis');
  bOrb.classList.toggle('on', st === 'orbital');
  if (current) { swapStructure(buildStructure(picture()), true); }
}
bLewis.addEventListener('click', () => setStyle('lewis'));
bOrb.addEventListener('click', () => setStyle('orbital'));
document.getElementById('prev').addEventListener('click', () => stepFamily(-1));
document.getElementById('next').addEventListener('click', () => stepFamily(1));
window.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
  if (e.key === 'ArrowLeft') stepFamily(-1);
  if (e.key === 'ArrowRight') stepFamily(1);
});

const camRight = new THREE.Vector3(), camUp = new THREE.Vector3(), toCam = new THREE.Vector3();
function placeSprites() {
  camera.matrixWorld.extractBasis(camRight, camUp, toCam);
  const pos = current.geo.pos;
  for (const f of faceCamera) {
    const p = pos[f.index];
    const rad = atomRadius(current.node.atoms[f.index].el);
    toCam.copy(camera.position).sub(new THREE.Vector3(...p)).normalize();
    if (f.kind === 'symbol') {
      f.sprite.position.set(p[0] + toCam.x * rad * 1.25, p[1] + toCam.y * rad * 1.25, p[2] + toCam.z * rad * 1.25);
    } else {
      const o = new THREE.Vector3(...p).add(toCam.clone().multiplyScalar(rad)).add(camRight.clone().multiplyScalar(rad * 0.95 + 0.08)).add(camUp.clone().multiplyScalar(rad * 0.85 + 0.06));
      f.sprite.position.copy(o);
    }
  }
}

function animate(t) {
  requestAnimationFrame(animate);
  if (current) {
    if (fade) {
      const f = Math.min(1, (t - fade.start) / FADE);
      const s = f * f * (3 - 2 * f);
      setOpacity(fade.from, 1 - s);
      setOpacity(fade.to, s);
      if (f >= 1) {
        scene.remove(fade.from);
        for (let k = faceCamera.length - 1; k >= 0; k--) if (faceCamera[k].group === fade.from) faceCamera.splice(k, 1);
        structureGroup = fade.to;
        fade = null;
      }
    }
    const family = current.families[familyIdx];
    if (mode === 'cycle' && family.length > 1 && !fade && t - lastSwitch > DWELL) {
      memberIdx = (memberIdx + 1) % family.length;
      lastSwitch = t;
      swapStructure(buildStructure(picture()), true);
      readout();
    }
    placeSprites();
  }
  controls.update();
  renderer.render(scene, camera);
}
tryLoad(ENTRIES.find(e => e.name === 'SF6') || ENTRIES[0]);
sel.value = ENTRIES.findIndex(e => e.name === 'SF6');
requestAnimationFrame(animate);
