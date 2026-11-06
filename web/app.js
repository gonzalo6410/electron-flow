// Electron Flow viewer: one Lewis structure drawn in 3D. Renders what the core computes; nothing
// chemical is decided here. Ledger items 13, 14, 24, 34.
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { CASES, runCase } from './core/cases.js';
import { COMPOUNDS } from './core/compounds.js';
import { nodeFromSmiles, toSmiles } from './core/smiles.js';
import { buildNode, enumerate } from './core/lewis.js';
import { relax } from './core/geometry.js';
import { piFrame } from './core/pi.js';
import { structureFamilies } from './core/symmetry.js';
import { ELEMENTS } from './core/elements.js';
import { attachEditor } from './edit.js';

// ---------- timing ----------
const DWELL = 4000;   // ms a structure stays before the next equivalent one
const FADE = 1200;    // ms for the old picture to go and the new one to come, one after the other
const MOVE = 800;     // ms for the camera to travel to a newly loaded molecule

// ---------- sizes, all in angstrom ----------
const LINE = 0.055;       // radius of every bond line, single or one of a double
const LINE_GAP = 0.21;    // centre to centre between the lines of a multiple bond
const HASH = 0.075;       // radius of the thin discs of a hashed line
const DOT_OUT = 0.12;     // a lone pair sits this far outside the ball's rim
const BADGE_OUT = 0.12;   // a charge badge sits this far outside the rim
// a loose cation is drawn the size of a chlorine ball, not of its metal atom, whose covalent radius
// would make it the biggest thing on screen
const atomRadius = el => 0.16 + 0.15 * (ELEMENTS[el].ionic ? Math.min(ELEMENTS[el].radius, 1.0) : ELEMENTS[el].radius);

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
// a bond order in words
function orderWords(o) {
  if (Math.abs(o - 1) < 1e-9) return 'single bond';
  if (Math.abs(o - 2) < 1e-9) return 'double bond';
  if (Math.abs(o - 3) < 1e-9) return 'triple bond';
  for (let d = 2; d <= 12; d++) {
    const num = Math.round(o * d);
    if (Math.abs(o * d - num) < 1e-6) return `bond of order ${num}/${d}`;
  }
  return `bond of order ${o.toFixed(2)}`;
}

// ---------- scene ----------
const view = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ antialias: true });
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff);
const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 100);
// the mouse turns the molecule about the camera's up axis; OrbitControls reads that axis once when
// it is made, so a fresh one is made whenever a molecule sets a new up
let controls = null;
function newControls(target) {
  if (controls) controls.dispose();
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.target.copy(target);
  controls.update();
}
newControls(new THREE.Vector3());
// the light travels with the camera: from the front, a little above and to the right of the viewer,
// so every molecule is lit the same way whichever way it was placed; the even light keeps the dark
// side of a ball from going black and the sun is weak enough that a white ball never burns out
scene.add(new THREE.AmbientLight(0xffffff, 1.1));
const sun = new THREE.DirectionalLight(0xffffff, 1.6);
sun.position.set(1.2, 2.0, 4.0);
sun.target.position.set(0, 0, 0);
camera.add(sun, sun.target);
scene.add(camera);

// size and sharpness follow the box the canvas sits in and the screen it is shown on
function resize() {
  const w = view.clientWidth, h = view.clientHeight;
  if (!w || !h) return;
  renderer.setPixelRatio(window.devicePixelRatio || 1);
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
new ResizeObserver(resize).observe(view);
// a window moved to a screen of another density changes devicePixelRatio without a resize
function watchPixelRatio() {
  const mq = window.matchMedia(`(resolution: ${window.devicePixelRatio}dppx)`);
  mq.addEventListener('change', () => { resize(); watchPixelRatio(); }, { once: true });
}
watchPixelRatio();
resize();

// ---------- shapes ----------
const sphereGeo = new THREE.SphereGeometry(1, 48, 32);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 24, 1);
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

// a hashed line: a faint thread with thin discs across it, for a pair of atoms adjacent but
// unbonded in this structure; it must read as "no bond here", quieter than any bond
function hashed(a, b, radius, material, group) {
  const d = b.clone().sub(a);
  const L = d.length();
  const step = 0.12;
  const n = Math.max(2, Math.floor(L / step));
  const start = (L - (n - 1) * step) / 2;
  const u = d.normalize();
  group.add(cylinder(a, b, 0.012, material));
  for (let k = 0; k < n; k++) {
    const p = a.clone().addScaledVector(u, start + k * step);
    group.add(cylinder(p.clone().addScaledVector(u, -0.008), p.clone().addScaledVector(u, 0.008), radius, material));
  }
}

// ---------- sprites (always face the camera) ----------
// Pictures drawn once on a canvas and kept: one per element symbol, one per charge mark, one for a
// pair of dots. A sprite gets its own material (opacity is per material) but shares the picture.
const textures = new Map();
function texture(key, draw) {
  if (textures.has(key)) return textures.get(key);
  const c = document.createElement('canvas');
  c.width = 512; c.height = 256;
  draw(c.getContext('2d'));
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;   // the canvas colours are screen colours, not light values
  tex.anisotropy = renderer.capabilities.getMaxAnisotropy();
  textures.set(key, tex);
  return tex;
}

function sprite(tex, width, height, opacity = 1, depthTest = true) {
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, opacity, depthTest, depthWrite: false }));
  sp.scale.set(width, height, 1);
  sp.userData.baseOpacity = opacity;
  return sp;
}

const FONT = '"Helvetica Neue", Helvetica, Arial, sans-serif';

// black on a light ball, white on a dark one; the split is at the lightness where the two read equally
function symbolColor(hex) {
  const lin = v => { const c = v / 255; return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4; };
  const Y = 0.2126 * lin((hex >> 16) & 255) + 0.7152 * lin((hex >> 8) & 255) + 0.0722 * lin(hex & 255);
  const lightness = 116 * Math.cbrt(Y) - 16;
  return lightness > 58 ? '#111' : '#fff';
}

// element symbol; the sprite is sized so the letters stand about as tall as the ball's radius
function symbolSprite(el, rad) {
  const tex = texture(`symbol:${el}`, g => {
    g.font = `500 190px ${FONT}`;
    g.textAlign = 'center'; g.textBaseline = 'middle';
    g.fillStyle = symbolColor(ELEMENTS[el].color);
    g.fillText(el, 256, 136);
  });
  const h = 1.5 * rad;
  return sprite(tex, 2 * h, h);
}

// formal charge in a small white badge, legible over anything; the signs are drawn as strokes so
// that + and - come out the same size whatever the font does with them
function chargeSprite(text) {
  const tex = texture(`charge:${text}`, g => {
    const digits = text.replace(/[+−]/g, '');
    const sign = text.includes(MINUS) ? '-' : '+';
    g.font = `600 112px ${FONT}`;
    g.textAlign = 'left'; g.textBaseline = 'middle';
    const signW = 64, pad = 22;
    const dw = digits ? g.measureText(digits).width + 10 : 0;
    const h = 150, r = h / 2;
    const w = Math.max(h, dw + signW + 2 * pad);
    const x0 = 256 - w / 2, y0 = 128 - h / 2;
    g.beginPath();
    g.moveTo(x0 + r, y0); g.lineTo(x0 + w - r, y0); g.arc(x0 + w - r, y0 + r, r, -Math.PI / 2, Math.PI / 2);
    g.lineTo(x0 + r, y0 + h); g.arc(x0 + r, y0 + r, r, Math.PI / 2, 3 * Math.PI / 2); g.closePath();
    g.fillStyle = '#fff'; g.fill();
    g.lineWidth = 9; g.strokeStyle = '#111'; g.stroke();
    g.fillStyle = '#111';
    const left = 256 - (dw + signW) / 2;
    if (digits) g.fillText(digits, left, 134);
    const cx = left + dw + signW / 2, cy = 128;
    g.lineWidth = 13; g.lineCap = 'round';
    g.beginPath(); g.moveTo(cx - 26, cy); g.lineTo(cx + 26, cy); g.stroke();
    if (sign === '+') { g.beginPath(); g.moveTo(cx, cy - 26); g.lineTo(cx, cy + 26); g.stroke(); }
  });
  const sp = sprite(tex, 0.56, 0.28, 1, false);
  sp.renderOrder = 10;
  return sp;
}

// a lone pair: two dots
function pairSprite(opacity) {
  const tex = texture('pair', g => {
    g.fillStyle = '#111';
    for (const x of [184, 328]) { g.beginPath(); g.arc(x, 128, 44, 0, Math.PI * 2); g.fill(); }
  });
  return sprite(tex, 0.32, 0.16, opacity);
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
  const spread = order.map(k => m[k][k]);
  const flat = spread[2] < 0.5 * spread[0];
  const linear = spread[1] < 0.02 * spread[0];
  return { major, middle, minor, flat, linear };
}

// Where the camera should stand for a molecule: looking at the face of a flat one, along the lines of
// a straight one, obliquely at a round one; far enough that the whole thing fits with a margin.
function frameFor(node, geo, frame) {
  const pos = geo.pos.map(p => new THREE.Vector3(...p));
  const box = new THREE.Box3();
  pos.forEach((p, i) => { const r = atomRadius(node.atoms[i].el); box.expandByPoint(p.clone().addScalar(r)); box.expandByPoint(p.clone().addScalar(-r)); });
  const target = box.getCenter(new THREE.Vector3());
  // the ball that holds every atom with its dots and badge; a small molecule is not blown up to fill the screen
  const radius = Math.max(2.4, ...pos.map((p, i) => p.distanceTo(target) + atomRadius(node.atoms[i].el) + 0.35));
  const vfov = (camera.fov / 2) * Math.PI / 180;
  const hfov = Math.atan(Math.tan(vfov) * camera.aspect);
  const dist = radius / Math.sin(Math.min(vfov, hfov)) * 1.12;

  const heavy = geo.pos.filter((p, i) => node.atoms[i].el !== 'H');
  const axes = principalAxes(heavy.length >= 3 ? heavy : geo.pos);
  let dir, up;
  if (axes.linear) {
    // a straight molecule: the lines of a multiple bond lie across the direction its pi part points
    // in, so looking along that direction shows them side by side; with two such directions at right
    // angles (a cumulene) the view splits the difference; the long axis runs across the screen
    const look = new THREE.Vector3();
    frame.normals.forEach(ns => {
      if (!ns.length) return;
      const v = new THREE.Vector3(...ns[0]);
      look.add(look.dot(v) < 0 ? v.negate() : v);
    });
    if (look.lengthSq() < 1e-6) look.copy(axes.middle);
    look.addScaledVector(axes.major, -look.dot(axes.major));
    dir = look.normalize();
    up = dir.clone().cross(axes.major).normalize();
    if (up.y < 0) up.negate();
  } else if (axes.flat) {
    // a flat molecule: look at its face, tilted a little so depth still shows
    dir = axes.minor.clone().multiplyScalar(Math.cos(0.35)).addScaledVector(axes.major, Math.sin(0.35) * 0.6).addScaledVector(axes.middle, Math.sin(0.35) * 0.8).normalize();
    up = axes.middle.clone();
  } else {
    dir = new THREE.Vector3(0.38, 0.42, 0.82).normalize(); // a round molecule: off every axis
    up = new THREE.Vector3(0, 1, 0);
  }
  up.addScaledVector(dir, -up.dot(dir)).normalize();
  return { position: target.clone().addScaledVector(dir, dist), up, target };
}

// The camera travels from where it is to the new frame: the viewing direction turns at a steady rate,
// the distance and the point looked at change smoothly, nothing jumps.
let move = null;   // { from, to, start }
function moveCamera(to, animate) {
  const from = { position: camera.position.clone(), up: camera.up.clone(), target: controls.target.clone() };
  if (!animate) { placeCamera(to); return; }
  move = { from, to, start: performance.now() };
  controls.enabled = false;
}
function placeCamera(f) {
  camera.position.copy(f.position);
  camera.up.copy(f.up);
  camera.lookAt(f.target);
  newControls(f.target);
}
function stepMove(t) {
  const f = Math.min(1, (t - move.start) / MOVE);
  const e = f * f * (3 - 2 * f);
  const { from, to } = move;
  const target = from.target.clone().lerp(to.target, e);
  const d0 = from.position.clone().sub(from.target), d1 = to.position.clone().sub(to.target);
  const dist = d0.length() + (d1.length() - d0.length()) * e;
  const q = new THREE.Quaternion().setFromUnitVectors(d0.clone().normalize(), d1.clone().normalize());
  const dir = d0.clone().normalize().applyQuaternion(new THREE.Quaternion().slerp(q, e));
  const up = from.up.clone().lerp(to.up, e).normalize();
  camera.position.copy(target).addScaledVector(dir, dist);
  camera.up.copy(up);
  camera.lookAt(target);
  if (f >= 1) {
    placeCamera(to);
    move = null;
  }
}

// ---------- state ----------
let current = null;          // { node, r, geo, families, main }
let showResonance = false;
let paused = false;
let familyIdx = 0, memberIdx = 0;
let lastSwitch = 0;
let fade = null;             // { from, to, start }
const atomGroup = new THREE.Group();      // spheres and symbols, fixed per molecule
let structureGroup = new THREE.Group();   // bonds, dots, charges, per structure
scene.add(atomGroup, structureGroup);
const tracked = [];          // sprites placed relative to the camera every frame
let editor = null;           // set up at the end, once the page exists

function load(entry, animate) {
  if (entry.kind === 'case') {
    const { node, r } = runCase(entry.c);
    show(null, node, r, null, animate);
  } else {
    loadSpec(nodeFromSmiles(entry.smiles), null, animate);
  }
}

// a molecule given as atoms, adjacency and charge (from a SMILES or from editing); `start` holds
// positions to begin the relaxation from, used when editing so the picture does not jump
function loadSpec(spec, start, animate = true) {
  if (!spec.atoms.length) { showEmpty(); return; }
  const node = buildNode(spec);
  const r = enumerate(node);
  show(spec, node, r, start, animate);
}

function showEmpty() {
  current = null;
  fade = null;
  disposeGroup(atomGroup); atomGroup.clear(); tracked.length = 0;
  retire(structureGroup); structureGroup = new THREE.Group(); scene.add(structureGroup);
  document.getElementById('steps').style.visibility = 'hidden';
  document.getElementById('status').innerHTML = '<div>nothing here yet</div>';
  hideTip();
  if (editor) editor.onLoaded();
}

// true when some atom sits outside the picture
function outOfView(geo) {
  camera.updateMatrixWorld();
  const v = new THREE.Vector3();
  return geo.pos.some(p => { v.set(...p).project(camera); return Math.abs(v.x) > 0.92 || Math.abs(v.y) > 0.92 || v.z > 1; });
}

function show(spec, node, r, start, animate) {
  const geo = relax(node, r, start ? { start } : {});
  // the main structure: least total formal charge, first one on a tie (all weights are equal; this
  // only picks what to show when resonance is switched off, and which family comes first)
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  let main = 0;
  r.structures.forEach((st, q) => { if (spread(st) < spread(r.structures[main])) main = q; });
  const families = structureFamilies(node, r).families
    .sort((a, b) => (spread(r.structures[a[0]]) - spread(r.structures[b[0]])) || (b.length - a.length) || (a[0] - b[0]));
  current = { spec, node, r, geo, families, main };
  familyIdx = 0; memberIdx = 0; fade = null;
  lastSwitch = performance.now();
  if (!start || outOfView(geo)) moveCamera(frameFor(node, geo, piFrame(node, r, geo.pos, geo.groups, r.structures[main])), animate);
  buildAtoms();
  swapStructure(buildStructure(), false);
  hideTip();
  updatePanel();
  if (editor && !start) editor.onLoaded();
  if (editor) editor.refresh();
}

function shownStructure() {
  const { r, families, main } = current;
  return showResonance ? r.structures[families[familyIdx][memberIdx]] : r.structures[main];
}

// free the materials of a group that leaves the scene (geometries and pictures are shared and kept)
function disposeGroup(group) {
  group.traverse(o => { if (o.material) o.material.dispose(); });
}

function buildAtoms() {
  disposeGroup(atomGroup);
  atomGroup.clear();
  tracked.length = 0;
  const { node, geo } = current;
  node.atoms.forEach((a, i) => {
    const rad = atomRadius(a.el);
    const m = new THREE.Mesh(sphereGeo, solid(ELEMENTS[a.el].color));
    m.position.set(...geo.pos[i]);
    m.scale.setScalar(rad);
    m.userData.atom = i;
    atomGroup.add(m);
    const sp = symbolSprite(a.el, rad);
    atomGroup.add(sp);
    tracked.push({ sprite: sp, index: i, kind: 'symbol' });
    // a loose anion (a free chloride) keeps its full shell: four pairs around it, as on paper
    const fixed = node.fixed[i];
    if (fixed && fixed.lone > 0) {
      const pairs = fixed.lone / 2;
      for (let q = 0; q < pairs; q++) {
        const dot = pairSprite(1);
        atomGroup.add(dot);
        tracked.push({ sprite: dot, index: i, kind: 'pair', dir: null, angle: Math.PI / 4 + (q * 2 * Math.PI) / pairs, group: atomGroup });
      }
    }
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
  const bondMat = solid(0x3a3a3a), ghostMat = solid(0xb4b4b4), ionMat = solid(0xc8c8c8);

  // every part gets a key; parts with the same key in two structures are the same picture
  const part = (key, fn) => { const sub = new THREE.Group(); sub.userData.key = key; fn(sub); g.add(sub); };
  // what the pointer finds when it rests on a piece of a bond
  const tag = (sub, k, info) => sub.traverse(o => { if (o.isMesh) { o.userData.edge = k; o.userData.info = info; } });

  r.edges.forEach(([i, j], k) => {
    const o = st.orders[k];
    const a = pos[i], b = pos[j];
    const bd = b.clone().sub(a).normalize();
    const a2 = a.clone().addScaledVector(bd, atomRadius(sym[i]) * 0.9);
    const b2 = b.clone().addScaledVector(bd, -atomRadius(sym[j]) * 0.9);
    const pair = `${sym[i]}–${sym[j]}`;
    if (o <= 1e-9) {
      part(`edge:${k}:0`, sub => {
        if (node.fixed[i] || node.fixed[j]) { dashed(a2, b2, 0.035, ionMat, sub); tag(sub, k, `${pair}: an ion kept nearby, no bond`); }
        else { hashed(a2, b2, HASH, ghostMat, sub); tag(sub, k, `${pair}: no bond in this structure`); }
      });
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
    part(`edge:${k}:${o.toFixed(3)}:${side.x.toFixed(2)},${side.y.toFixed(2)},${side.z.toFixed(2)}`, sub => {
      for (let q = 0; q < lines; q++) {
        const shift = side.clone().multiplyScalar(LINE_GAP * (q - (lines - 1) / 2));
        const pa = a2.clone().add(shift), pb = b2.clone().add(shift);
        if (q === lines - 1 && fr > 1e-9) dashed(pa, pb, LINE * 0.8, bondMat, sub);
        else sub.add(cylinder(pa, pb, LINE, bondMat));
      }
      tag(sub, k, `${pair}: ${orderWords(o)}`);
    });
  });

  // lone pairs: each pair gets a fixed direction in space (the in-plane slots and the partner
  // directions from the pi frame; pairs that would sit above or below the plane are put in the
  // widest gaps between the bonds in the plane instead, as dots on a drawing would be); every frame
  // that direction is projected onto the screen, so the dots move smoothly as the molecule turns
  sym.forEach((el, i) => {
    if (node.fixed[i]) return;
    const pairs = st.lone[i] / 2;
    const whole = Math.floor(pairs + 1e-9);
    const lps = frame.lonePairs[i];
    if (!lps.length) return;
    const dirs = drawingDirections(i, lps, frame, pos);
    const partners = lps.filter(lp => lp.kind === 'partner').map(lp => lp.partner);
    part(`pairs:${i}:${st.lone[i]}:${partners.join(',')}`, sub => {
      dirs.forEach((dir, q) => {
        const sp = pairSprite(q < whole ? 1 : Math.max(0.2, pairs - whole));
        sub.add(sp);
        tracked.push({ sprite: sp, index: i, kind: 'pair', dir, angle: null, group: g });
      });
    });
  });

  // charges
  sym.forEach((el, i) => {
    if (Math.abs(st.charge[i]) < 1e-9) return;
    part(`charge:${i}:${st.charge[i]}`, sub => {
      const sp = chargeSprite(chargeMark(st.charge[i]));
      sub.add(sp);
      tracked.push({ sprite: sp, index: i, kind: 'charge', angle: null, goal: null, group: g });
    });
  });
  return g;
}

// directions in space for an atom's dot pairs
function drawingDirections(i, lps, frame, pos) {
  const kept = lps.filter(lp => lp.kind !== 'p' && lp.kind !== 'free').map(lp => new THREE.Vector3(...lp.dir));
  const need = lps.length - kept.length;
  if (!need) return kept;
  const bonds = current.geo.groups[i].nbrs.map(nb => pos[nb.j].clone().sub(pos[i]).normalize());
  const extra = [];
  const pO = frame.pOrbital[i];
  if (pO && bonds.length >= 2) {
    // in the plane of the atom, in the widest gaps between its bonds and the pairs already placed
    const nrm = new THREE.Vector3(...pO);
    const ex = bonds[0].clone().addScaledVector(nrm, -bonds[0].dot(nrm)).normalize();
    const ey = nrm.clone().cross(ex);
    const placed = [...bonds, ...kept];
    for (let q = 0; q < need; q++) {
      const angles = placed.map(v => Math.atan2(v.dot(ey), v.dot(ex))).sort((a, b) => a - b);
      let best = { mid: 0, size: -1 };
      for (let k = 0; k < angles.length; k++) {
        const a = angles[k], b = k + 1 < angles.length ? angles[k + 1] : angles[0] + 2 * Math.PI;
        if (b - a > best.size) best = { mid: (a + b) / 2, size: b - a };
      }
      const v = ex.clone().multiplyScalar(Math.cos(best.mid)).addScaledVector(ey, Math.sin(best.mid));
      extra.push(v); placed.push(v);
    }
  } else {
    // nothing flat to go by: spread them away from the bonds
    let away = bonds.reduce((s, v) => s.add(v), new THREE.Vector3()).negate();
    if (away.length() < 1e-6) away.set(0, 1, 0);
    away.normalize();
    const p1 = Math.abs(away.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
    p1.addScaledVector(away, -p1.dot(away)).normalize();
    const p2 = away.clone().cross(p1);
    for (let q = 0; q < need; q++) {
      const ang = (2 * Math.PI * q) / need;
      extra.push(away.clone().addScaledVector(p1, 0.7 * Math.cos(ang)).addScaledVector(p2, 0.7 * Math.sin(ang)).normalize());
    }
  }
  return [...kept, ...extra];
}

function setOpacity(group, f) {
  group.traverse(o => { if (o.material) o.material.opacity = f * (o.userData.baseOpacity ?? 1); });
}

function dropTracked(group) {
  for (let k = tracked.length - 1; k >= 0; k--) if (tracked[k].group === group) tracked.splice(k, 1);
}

function retire(group) {
  scene.remove(group);
  dropTracked(group);
  disposeGroup(group);
}

function finishFade() {
  if (!fade) return;
  retire(fade.from);
  fade.to.children.forEach(c => { c.visible = true; });
  setOpacity(fade.to, 1);
  structureGroup = fade.to;
  fade = null;
}

// Replace the picture. With animate, parts that are the same in both structures stay as they are;
// a part that goes fades out while the part that replaces it fades in, in place.
function swapStructure(next, animate) {
  finishFade();
  if (!animate) {
    retire(structureGroup);
    structureGroup = next; scene.add(next);
    return;
  }
  const oldKeys = new Set(structureGroup.children.map(c => c.userData.key));
  const newKeys = new Set(next.children.map(c => c.userData.key));
  const going = structureGroup.children.filter(c => !newKeys.has(c.userData.key));
  const coming = next.children.filter(c => !oldKeys.has(c.userData.key));
  // parts common to both: the old copy goes now, the new copy shows at once
  structureGroup.children.filter(c => newKeys.has(c.userData.key)).forEach(c => { c.visible = false; });
  coming.forEach(c => { setOpacity(c, 0); c.visible = true; });
  scene.add(next);
  fade = { from: structureGroup, to: next, going, coming, start: performance.now() };
}

function showStructure(animate) {
  swapStructure(buildStructure(), animate);
  updatePanel();
}

// ---------- per-frame sprite placement ----------
const camRight = new THREE.Vector3(), camUp = new THREE.Vector3(), camBack = new THREE.Vector3();
const TWO_PI = 2 * Math.PI;
const wrap = a => ((a % TWO_PI) + TWO_PI) % TWO_PI;            // into [0, 2 pi)
const turn = (a, b) => { const d = wrap(b - a); return d > Math.PI ? d - TWO_PI : d; };   // shortest signed turn from a to b
const gap = (a, b) => Math.abs(turn(a, b));                   // angle between two directions

// screen angle (right = 0, up = pi/2) of a direction or of an atom seen from another
const screenAngle = v => Math.atan2(v.dot(camUp), v.dot(camRight));

// Dots of one atom that project to nearly the same spot are nudged apart along the rim, so two pairs
// never read as one cluster of four; the nudge is the same every frame, so the dots move smoothly.
function spreadAngles(angles, minSep) {
  const order = angles.map((a, q) => q).sort((p, q) => wrap(angles[p]) - wrap(angles[q]));
  for (let it = 0; it < 12; it++) {
    let moved = false;
    for (let k = 0; k < order.length; k++) {
      const p = order[k], q = order[(k + 1) % order.length];
      const d = order.length === 1 ? TWO_PI : wrap(angles[q] - angles[p]);
      if (d < minSep - 1e-6) { const push = (minSep - d) / 2; angles[p] -= push; angles[q] += push; moved = true; }
    }
    if (!moved) break;
  }
  return angles;
}

// A charge badge prefers the upper right of its atom, as on paper, and moves to the clearest spot
// around the rim only when a bond or a pair comes to lie under it; the move is eased, never a jump.
function badgeAngle(t, obstacles, dt) {
  const clearance = a => obstacles.reduce((m, o) => Math.min(m, gap(a, o)), Math.PI);
  const preferred = Math.PI / 4;
  if (t.goal === null || clearance(t.goal) < 0.4) {
    let best = preferred, bestScore = -Infinity;
    for (let k = 0; k < 24; k++) {
      const a = preferred + (k * TWO_PI) / 24;
      const score = Math.min(clearance(a), 0.9) - 0.12 * gap(a, preferred);
      if (score > bestScore) { bestScore = score; best = a; }
    }
    t.goal = best;
  }
  if (t.angle === null) t.angle = t.goal;
  const d = turn(t.angle, t.goal);
  const step = Math.min(Math.abs(d), 5 * dt);   // radians per second
  t.angle = wrap(t.angle + Math.sign(d) * step);
  return t.angle;
}

function placeSprites(dt) {
  camera.matrixWorld.extractBasis(camRight, camUp, camBack);
  const { node, r, geo } = current;
  const pos = geo.pos.map(p => new THREE.Vector3(...p));
  const toCam = new THREE.Vector3();
  // per atom: the screen angles of its bonds and of its dots, for the badges to keep clear of
  const bondAngles = pos.map(() => []);
  r.edges.forEach(([i, j]) => {
    const v = pos[j].clone().sub(pos[i]);
    bondAngles[i].push(screenAngle(v)); bondAngles[j].push(screenAngle(v.negate()));
  });
  const dotAngles = pos.map(() => []);
  const pairs = tracked.filter(t => t.kind === 'pair');
  for (const t of pairs) {
    if (t.dir) {
      const x = t.dir.dot(camRight), y = t.dir.dot(camUp);
      if (Math.hypot(x, y) > 0.2 || t.angle === null) t.angle = Math.atan2(y, x);
    }
    t.shown = t.angle;
  }
  for (let i = 0; i < pos.length; i++) {
    const mine = pairs.filter(t => t.index === i);
    if (mine.length < 2) continue;
    const rad = atomRadius(node.atoms[i].el);
    const spread = spreadAngles(mine.map(t => t.angle), 0.24 / (rad + DOT_OUT));
    mine.forEach((t, q) => { t.shown = spread[q]; });
  }
  for (const t of tracked) {
    const p = pos[t.index];
    const rad = atomRadius(node.atoms[t.index].el);
    toCam.copy(camera.position).sub(p).normalize();
    if (t.kind === 'symbol') {
      t.sprite.position.copy(p).addScaledVector(toCam, rad * 1.25);
    } else if (t.kind === 'pair') {
      const d = rad + DOT_OUT;
      t.sprite.position.copy(p).addScaledVector(toCam, rad * 1.1).addScaledVector(camRight, Math.cos(t.shown) * d).addScaledVector(camUp, Math.sin(t.shown) * d);
      t.sprite.material.rotation = t.shown + Math.PI / 2;   // the two dots lie along the rim
      dotAngles[t.index].push(t.shown);
    } else {
      const a = badgeAngle(t, [...bondAngles[t.index], ...dotAngles[t.index]], dt);
      const d = rad + BADGE_OUT;
      t.sprite.position.copy(p).addScaledVector(toCam, rad).addScaledVector(camRight, Math.cos(a) * d).addScaledVector(camUp, Math.sin(a) * d);
    }
  }
}

// ---------- hover: a small plain note on the atom or bond under the pointer ----------
const tip = document.getElementById('tip');
const pointer = { x: 0, y: 0, inside: false, down: false };
const raycaster = new THREE.Raycaster();
let tipShown = '';

renderer.domElement.addEventListener('pointermove', e => {
  const rect = renderer.domElement.getBoundingClientRect();
  pointer.x = e.clientX - rect.left; pointer.y = e.clientY - rect.top;
  pointer.inside = true;
});
renderer.domElement.addEventListener('pointerleave', () => { pointer.inside = false; hideTip(); });
renderer.domElement.addEventListener('pointerdown', () => { pointer.down = true; hideTip(); });
window.addEventListener('pointerup', () => { pointer.down = false; });

function hideTip() { tip.style.display = 'none'; tipShown = ''; }

function atomNote(i) {
  const { node } = current;
  const st = shownStructure();
  const el = node.atoms[i].el;
  const q = st.charge[i], lone = st.lone[i];
  const charge = Math.abs(q) < 1e-9 ? 'no charge' : `charge ${chargeMark(q)}`;
  if (node.fixed[i]) return `<b>${el}</b> ${charge}<br>a loose ion, no bonds${lone ? `, ${lone} lone electrons` : ''}`;
  const ceiling = ELEMENTS[el].ceiling;
  const count = st.count[i];
  const octet = count >= ceiling ? `counts ${count} electrons, full` : `counts ${count} electrons, ${ceiling - count} short`;
  const lones = lone === 0 ? 'no lone electrons' : lone === 1 ? '1 lone electron' : `${lone} lone electrons`;
  return `<b>${el}</b> ${charge}<br>${octet}, ${lones}`;
}

function updateTip() {
  if (!current || !pointer.inside || pointer.down || move) { if (tipShown) hideTip(); return; }
  const rect = renderer.domElement.getBoundingClientRect();
  raycaster.setFromCamera(new THREE.Vector2((pointer.x / rect.width) * 2 - 1, -(pointer.y / rect.height) * 2 + 1), camera);
  const meshes = [];
  atomGroup.children.forEach(o => { if (o.isMesh) meshes.push(o); });
  structureGroup.traverse(o => { if (o.isMesh && o.visible && o.material.opacity > 0.5) meshes.push(o); });
  const hit = raycaster.intersectObjects(meshes, false)[0];
  let text = '';
  if (hit) text = hit.object.userData.atom !== undefined ? atomNote(hit.object.userData.atom) : hit.object.userData.info;
  if (!text) { if (tipShown) hideTip(); return; }
  if (text !== tipShown) { tip.innerHTML = text; tip.style.display = 'block'; tipShown = text; }
  const x = Math.min(pointer.x + 14, rect.width - tip.offsetWidth - 6);
  const y = Math.min(pointer.y + 16, rect.height - tip.offsetHeight - 6);
  tip.style.left = `${x}px`; tip.style.top = `${y}px`;
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
function tryLoad(entry, animate = true) {
  try { load(entry, animate); errBox.textContent = ''; }
  catch (e) { errBox.textContent = e.message; }
}
sel.addEventListener('change', () => {
  if (sel.value === 'typed') return;
  const entry = ENTRIES[Number(sel.value)];
  tryLoad(entry);
  smilesBox.value = entry.kind === 'smiles' ? entry.smiles : '';
});
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
  resonanceBox.blur();   // so the keys go back to the page
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
function togglePause() {
  paused = !paused;
  lastSwitch = performance.now();
  updatePanel();
}
for (const [id, d] of [['prev', -1], ['next', 1]]) {
  const b = document.getElementById(id);
  b.addEventListener('click', () => { stepFamily(d); b.blur(); });
}
// keys: left and right step the structures, space pauses the alternation; nothing while typing
window.addEventListener('keydown', e => {
  if (['INPUT', 'SELECT', 'TEXTAREA', 'BUTTON'].includes(e.target.tagName) || e.metaKey || e.ctrlKey || e.altKey) return;
  if (e.key === 'ArrowLeft') { stepFamily(-1); e.preventDefault(); }
  else if (e.key === 'ArrowRight') { stepFamily(1); e.preventDefault(); }
  else if (e.key === ' ') { togglePause(); e.preventDefault(); }
});

// the status: four short lines, always the same height so nothing below it moves
function updatePanel() {
  const { r, families } = current;
  const n = r.structures.length;
  const stepping = showResonance && families.length > 1;
  document.getElementById('steps').style.visibility = stepping ? 'visible' : 'hidden';
  const lines = [`<b>${r.pool}</b> valence electrons` + (current.spec && current.spec.charge ? ` · charge <b>${(current.spec.charge > 0 ? '+' : MINUS) + Math.abs(current.spec.charge)}</b>` : ''), `<b>${n}</b> contributing structure${n === 1 ? '' : 's'}`];
  if (r.transitionEdges) lines.push('before and after, equal weight');
  if (showResonance && n > 1) {
    const size = families[familyIdx].length;
    if (families.length > 1) lines.push(`structure <b>${familyIdx + 1}</b> of ${families.length}`);
    if (size > 1) lines.push(`<b>${size}</b> equivalent, ${paused ? 'paused' : 'alternating'}`);
  } else if (n > 1) lines.push('main structure shown');
  document.getElementById('status').innerHTML = lines.map(l => `<div>${l}</div>`).join('');
}

// ---------- loop ----------
let lastT = 0;
function animate(t) {
  requestAnimationFrame(animate);
  const dt = Math.min(0.1, (t - lastT) / 1000);
  lastT = t;
  if (move) stepMove(t);
  if (current) {
    if (fade) {
      const f = Math.min(1, (t - fade.start) / FADE);
      const e = f * f * (3 - 2 * f);
      fade.going.forEach(c => setOpacity(c, 1 - e));
      fade.coming.forEach(c => setOpacity(c, e));
      if (f >= 1) finishFade();
    }
    const family = current.families[familyIdx];
    if (showResonance && !paused && family.length > 1 && !fade && t - lastSwitch > DWELL) {
      memberIdx = (memberIdx + 1) % family.length;
      lastSwitch = t;
      showStructure(true);
    }
    placeSprites(dt);
    updateTip();
  }
  if (!move) controls.update();
  renderer.render(scene, camera);
}

// ---------- editing ----------
editor = attachEditor({
  scene, camera, domElement: renderer.domElement,
  getSpec: () => (current ? current.spec : { atoms: [], adj: [], charge: 0 }),
  getPositions: () => (current ? current.geo.pos : []),
  atomMeshes: () => atomGroup.children.filter(o => o.userData.atom !== undefined),
  isEditable: () => !current || !!current.spec,
  loadSpec: (spec, start) => {
    loadSpec(spec, start);
    typedOption.textContent = 'edited molecule'; typedOption.hidden = false; sel.value = 'typed';
    // the box shows the SMILES of what is on screen now, in its main structure, ready to copy
    smilesBox.value = current ? toSmiles(current.node, current.r.structures[current.main]) : '';
    errBox.textContent = '';
  },
  setNote: text => { document.getElementById('note').textContent = text; },
});
document.getElementById('fill').addEventListener('click', e => { editor.fill(); e.target.blur(); });
document.getElementById('clear').addEventListener('click', e => { editor.clear(); e.target.blur(); });
// a periodic table of the elements the engine knows, in their places (main group only for now)
const TABLE = [
  ['H', '', '', '', '', '', '', 'He'],
  ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
  ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
  ['K', 'Ca', '', '', '', 'Se', 'Br', 'Kr'],
  ['', '', '', '', '', '', 'I', 'Xe'],
];
const tableBox = document.getElementById('table');
for (const row of TABLE) for (const el of row) {
  const b = document.createElement('button');
  if (!el) { b.className = 'gap'; b.tabIndex = -1; tableBox.appendChild(b); continue; }
  b.textContent = el;
  b.draggable = true;
  b.addEventListener('click', () => { editor.addElement(el, editor.selectedAtom()); b.blur(); });
  b.addEventListener('dragstart', e => { e.dataTransfer.setData('text/element', el); e.dataTransfer.effectAllowed = 'copy'; });
  tableBox.appendChild(b);
}

const first = ENTRIES.findIndex(e => e.name === 'SF6');
sel.value = first;
tryLoad(ENTRIES[first], false);
requestAnimationFrame(animate);
