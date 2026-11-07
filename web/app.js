// Electron Flow viewer: one Lewis structure drawn in 3D. Renders what the core computes; nothing
// chemical is decided here. Ledger items 13, 14, 24, 34.
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { CASES, runCase } from './core/cases.js';
import { nodeFromSmiles, toSmiles } from './core/smiles.js';
import { buildNode, enumerate } from './core/lewis.js';
import { relax } from './core/geometry.js';
import { piFrame } from './core/pi.js';
import { structureFamilies } from './core/symmetry.js';
import { ELEMENTS } from './core/elements.js';
import { attachEditor } from './edit.js';
import { attachArrows } from './arrows.js';
import { commit, filmSteps, combine } from './core/moves.js';
import { REACTIONS, resolveArrows, pieceFormulas } from './core/reactions.js';
import { makeIndex, search, namesFor, reasonText, balanceSmiles } from './core/names.js';

// ---------- timing ----------
const DWELL = 3000;   // ms a structure stays before the next equivalent one (4 s until 7 Nov 2026, his call: a bit faster)
const FADE = 1200;    // ms for the old picture to go and the new one to come, one after the other
const MOVE = 800;     // ms for the camera to travel to a newly loaded molecule
const MOVIE = 3200;   // ms for a move to play from the reagents to the products

// ---------- sizes, all in angstrom ----------
const LINE = 0.055;       // radius of every bond line, single or one of a double
const LINE_GAP = 0.21;    // centre to centre between the lines of a multiple bond
const HASH = 0.075;       // radius of the thin discs of a hashed line
const DOT_OUT = 0.12;     // a lone pair sits this far outside the ball's rim
const BADGE_OUT = 0.12;   // a charge badge sits this far outside the rim
// a loose cation is drawn the size of a chlorine ball, not of its metal atom, whose covalent radius
// would make it the biggest thing on screen
const atomRadius = el => 0.16 + 0.15 * (ELEMENTS[el].ionic ? Math.min(ELEMENTS[el].radius, 1.0) : ELEMENTS[el].radius);

// ---------- angles ----------
const TWO_PI = 2 * Math.PI, DEG = Math.PI / 180;
const wrap = a => ((a % TWO_PI) + TWO_PI) % TWO_PI;            // into [0, 2 pi)
const turn = (a, b) => { const d = wrap(b - a); return d > Math.PI ? d - TWO_PI : d; };   // shortest signed turn from a to b
const gap = (a, b) => Math.abs(turn(a, b));                   // angle between two directions

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
function pairSprite() {
  const tex = texture('pair', g => {
    g.fillStyle = '#111';
    for (const x of [184, 328]) { g.beginPath(); g.arc(x, 128, 44, 0, Math.PI * 2); g.fill(); }
  });
  return sprite(tex, 0.32, 0.16);
}

// a lone electron: one dot, the size of one of a pair's
function dotSprite() {
  const tex = texture('dot', g => {
    g.fillStyle = '#111';
    g.beginPath(); g.arc(256, 128, 44, 0, Math.PI * 2); g.fill();
  });
  return sprite(tex, 0.32, 0.16);
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
// a straight one, obliquely at a round one; far enough that the whole thing fits with a margin. With
// `keep`, the camera's present direction and up are kept and only the distance and the point looked
// at change, so a molecule being edited does not tumble when it outgrows the view.
function frameFor(node, geo, frame, keep) {
  const pos = geo.pos.map(p => new THREE.Vector3(...p));
  const box = new THREE.Box3();
  pos.forEach((p, i) => { const r = atomRadius(node.atoms[i].el); box.expandByPoint(p.clone().addScalar(r)); box.expandByPoint(p.clone().addScalar(-r)); });
  const target = box.getCenter(new THREE.Vector3());
  let dir, up;
  if (keep) { dir = keep.dir.clone(); up = keep.up.clone(); } else { [dir, up] = viewFor(node, geo, frame); }
  up.addScaledVector(dir, -up.dot(dir)).normalize();

  // far enough that every atom, with its dots and badge, falls inside both view angles as seen from
  // this direction, with a margin; never closer than a 2.4 ball, so a water molecule is not blown up
  const right = up.clone().cross(dir);
  const vfov = (camera.fov / 2) * Math.PI / 180;
  const hfov = Math.atan(Math.tan(vfov) * camera.aspect);
  let dist = 2.4 / Math.sin(Math.min(vfov, hfov));
  pos.forEach((p, i) => {
    const q = p.clone().sub(target);
    const reach = atomRadius(node.atoms[i].el) + 0.35;
    dist = Math.max(dist, q.dot(dir) + (Math.abs(q.dot(right)) + reach) / Math.tan(hfov), q.dot(dir) + (Math.abs(q.dot(up)) + reach) / Math.tan(vfov));
  });
  dist *= 1.12;
  return { position: target.clone().addScaledVector(dir, dist), up, target };
}

// A unit direction and a ring of twelve directions tilted from it by `tilt`.
function around(dir, tilt) {
  const p1 = new THREE.Vector3(Math.abs(dir.x) < 0.9 ? 1 : 0, Math.abs(dir.x) < 0.9 ? 0 : 1, 0);
  p1.addScaledVector(dir, -p1.dot(dir)).normalize();
  const p2 = dir.clone().cross(p1);
  const ring = [dir.clone()];
  for (let k = 0; k < 12; k++) {
    const a = (k * TWO_PI) / 12;
    ring.push(dir.clone().multiplyScalar(Math.cos(tilt)).addScaledVector(p1, Math.sin(tilt) * Math.cos(a)).addScaledVector(p2, Math.sin(tilt) * Math.sin(a)));
  }
  return ring;
}

// Of several viewing directions, the one from which the two closest atoms on the screen are
// furthest apart.
function clearest(dirs, pos) {
  let best = dirs[0], bestApart = -1;
  for (const d of dirs) {
    let apart = Infinity;
    for (let i = 0; i < pos.length; i++) for (let j = i + 1; j < pos.length; j++) {
      const v = pos[j].clone().sub(pos[i]);
      apart = Math.min(apart, v.addScaledVector(d, -v.dot(d)).length());
    }
    if (apart > bestApart) { bestApart = apart; best = d; }
  }
  return best;
}

// The direction to look from, and which way is up, for a molecule seen for the first time.
function viewFor(node, geo, frame) {
  const pos = geo.pos.map(p => new THREE.Vector3(...p));
  // the molecule's own atoms decide the view (its heavy atoms, or all of them when there are fewer
  // than three), not the loose ions beside it, so a CH3+ with a Cl- nearby is still seen face-on
  const ownIndex = geo.pos.map((p, i) => i).filter(i => !node.fixed[i]);
  const own = ownIndex.map(i => geo.pos[i]);
  const heavy = ownIndex.filter(i => node.atoms[i].el !== 'H').map(i => geo.pos[i]);
  const axes = principalAxes(heavy.length >= 3 ? heavy : own.length ? own : geo.pos);
  // the directions of the tethers to loose ions, seen from the molecule's own centre
  const ownCentre = own.reduce((s, p) => s.add(new THREE.Vector3(...p)), new THREE.Vector3()).divideScalar(Math.max(1, own.length));
  const ties = geo.pos.filter((p, i) => node.fixed[i]).map(p => new THREE.Vector3(...p).sub(ownCentre).normalize());
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
    if (look.lengthSq() < 1e-6 && ties.length) look.copy(axes.major).cross(ties[0]); // a loose ion beside the line stays beside it
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
    // with loose ions beside it, look across their tethers so no ion hides behind the molecule: with
    // one tether the view loses its part along it; with more, the view goes along the direction the
    // tethers span least
    let tilt = 0.35;
    if (ties.length) {
      const t = principalAxes([...ties, ...ties.map(v => v.clone().negate())].map(v => v.toArray()));
      if (t.linear) dir.addScaledVector(t.major, -dir.dot(t.major));
      else { dir.copy(t.minor).multiplyScalar(Math.sign(dir.dot(t.minor)) || 1); tilt = 0.5; }
      if (dir.lengthSq() < 1e-6) dir.copy(t.minor);
      dir.normalize();
    }
    // then, of that direction and a ring of tilts around it, the one that keeps the atoms furthest
    // apart on the screen: straight along a tetrahedron's axis an atom hides behind the centre, and
    // straight along the normal of two tethers the atoms mirrored across their plane overlap
    dir.copy(clearest(around(dir, tilt), pos));
  }
  return [dir, up];
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
let arrows = null;           // the curved arrows, set up at the end too
let movie = null;            // a move playing: { node, steps, start, out, shown }

let cutNote = '';             // after an edit: the bonds the octet left no room for
function load(entry, animate) {
  movie = null;
  cutNote = '';
  if (entry.kind === 'case') {
    const { node, r } = runCase(entry.c);
    show(null, node, r, null, animate);
  } else {
    loadSpec(nodeFromSmiles(entry.smiles), null, animate);
  }
  if (current) current.reaction = entry.kind === 'reaction' ? entry.rx : null;
  showHint();
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

// true when some atom sits outside the picture; given the atoms, when some atom with its dots
// around it (its radius and a little more) would be cut by the edge (the bromide of SN1 step 1 left
// its dots below the picture with its centre inside)
function outOfView(geo, node = null) {
  camera.updateMatrixWorld();
  const v = new THREE.Vector3();
  const beyond = (q, edge) => { v.copy(q).project(camera); return Math.abs(v.x) > edge || Math.abs(v.y) > edge || v.z > 1; };
  if (!node) return geo.pos.some(p => beyond(new THREE.Vector3(...p), 0.92));
  camera.matrixWorld.extractBasis(camRight, camUp, camBack);
  return geo.pos.some((p, i) => {
    const c = new THREE.Vector3(...p), reach = atomRadius(node.atoms[i].el) + 0.35;
    return [camRight, camUp].some(d => beyond(c.clone().addScaledVector(d, reach), 0.98) || beyond(c.clone().addScaledVector(d, -reach), 0.98));
  });
}

// the main structure and the families in the order they are shown, for a molecule just loaded
function mainAndFamilies(node, r) {
  // the main structure: least total formal charge; on a tie, the one whose loose ions sit by atoms of
  // the opposite charge (a sodium of carbonate by an O-, not by the C=O); then the one with its
  // negative charges on the more electronegative atoms (the O- of pyridine N-oxide rather than a ring
  // C-); then the first written. All weights stay equal; this only picks what to show when resonance
  // is switched off, and which family comes first.
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  const ionsUnplaced = st => r.edges.reduce((s, [i, j]) => {
    const ion = node.fixed[i] ? i : node.fixed[j] ? j : -1;
    if (ion < 0) return s;
    return s + (st.charge[ion === i ? j : i] * node.fixed[ion].charge < 0 ? 0 : 1);
  }, 0);
  const chargeOnElectronegative = st => st.charge.reduce((s, q, i) => s + q * ELEMENTS[node.atoms[i].el].en, 0);
  const rank = st => [spread(st), ionsUnplaced(st), chargeOnElectronegative(st)];
  const before = (a, b) => { const ra = rank(a), rb = rank(b); return (ra[0] - rb[0]) || (ra[1] - rb[1]) || (ra[2] - rb[2]); };
  let main = 0;
  r.structures.forEach((st, q) => { if (before(st, r.structures[main]) < 0) main = q; });
  const families = structureFamilies(node, r).families
    .sort((a, b) => before(r.structures[a[0]], r.structures[b[0]]) || (b.length - a.length) || (a[0] - b[0]));
  return { main, families };
}

function show(spec, node, r, start, animate) {
  const geo = relax(node, r, start ? { start } : {});
  const { main, families } = mainAndFamilies(node, r);
  current = { spec, node, r, geo, families, main };
  familyIdx = 0; memberIdx = 0; fade = null;
  if (arrows) arrows.dropIfStale();
  if (!start) paused = false; // a new molecule alternates; an edit keeps the pause
  lastSwitch = performance.now();
  if (!start) moveCamera(frameFor(node, geo, piFrame(node, r, geo.pos, geo.groups, r.structures[main])), animate);
  else if (outOfView(geo)) moveCamera(frameFor(node, geo, null, { dir: camera.position.clone().sub(controls.target).normalize(), up: camera.up }), true);
  buildAtoms();
  swapStructure(buildStructure(), false);
  hideTip();
  updatePanel();
  if (editor && !start) editor.onLoaded();
  if (editor) editor.refresh();
}

function shownStructure() {
  if (movie) return movie.shown;
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
        const dot = pairSprite();
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
  const pickMat = new THREE.MeshBasicMaterial({ transparent: true, opacity: 0, depthWrite: false });

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
    // an unseen thicker line along every pair, so a drag can start from or land on a bond even
    // where its two drawn lines leave a gap between them
    part(`pick:${k}`, sub => {
      const m = cylinder(a2, b2, 0.16, pickMat);
      m.userData.edge = k; m.userData.pick = true; m.userData.baseOpacity = 0;
      sub.add(m);
    });
    // while a move plays: a bond being made or broken is drawn as faint as it is weak and not at all
    // once gone; an ion's old tie fades out as its new one fades in (core/moves.js, filmRoles)
    const role = movie ? movie.roles[k] : null;
    const changing = !!role && role.changing;
    if (role && (role.tieBefore || role.tieAfter)) {
      const f = (role.tieBefore ? 1 - movieMix() : 0) + (role.tieAfter ? movieMix() : 0);
      if (f > 0.02) {
        part(`tie:${k}:${f.toFixed(2)}`, sub => {
          const mat = f < 1 ? solid(0xc8c8c8) : ionMat;
          if (f < 1) mat.opacity = f;
          dashed(a2, b2, 0.035, mat, sub);
          sub.traverse(m => { if (m.isMesh) m.userData.baseOpacity = f; });
          tag(sub, k, `${pair}: an ion kept nearby, no bond`);
        });
      }
      if (o <= 1e-9) return;
    }
    if (o <= 1e-9 && changing) return; // a bond already broken, or not yet made: nothing to draw
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
    // while a move plays, a bond being made or broken is as faint as it is weak, up to half a bond
    const faint = changing && o < 0.5 ? Math.max(0.08, o / 0.5) : 1;
    const mat = faint < 1 ? solid(0x3a3a3a) : bondMat;
    if (faint < 1) mat.opacity = faint;
    part(`edge:${k}:${o.toFixed(3)}:${side.x.toFixed(2)},${side.y.toFixed(2)},${side.z.toFixed(2)}`, sub => {
      for (let q = 0; q < lines; q++) {
        const shift = side.clone().multiplyScalar(LINE_GAP * (q - (lines - 1) / 2));
        const pa = a2.clone().add(shift), pb = b2.clone().add(shift);
        if (q === lines - 1 && fr > 1e-9) dashed(pa, pb, LINE * 0.8, mat, sub);
        else sub.add(cylinder(pa, pb, LINE, mat));
      }
      if (faint < 1) sub.traverse(m => { if (m.isMesh) m.userData.baseOpacity = faint; });
      tag(sub, k, `${pair}: ${changing ? 'a bond being made or broken' : orderWords(o)}`);
    });
  });

  // lone pairs: each pair gets a fixed direction in space (see drawingDirections); every frame that
  // direction is projected onto the screen, so the dots move smoothly as the molecule turns. An odd
  // lone electron, the one of a radical, is a single dot.
  const heavy = pos.filter((p, i) => sym[i] !== 'H');
  const centre = (heavy.length ? heavy : pos).reduce((s, p) => s.add(p), new THREE.Vector3()).divideScalar(Math.max(1, heavy.length || pos.length));
  sym.forEach((el, i) => {
    if (node.fixed[i]) return;
    const lps = frame.lonePairs[i];
    if (!lps.length) return;
    const dirs = drawingDirections(i, lps, frame, pos, centre);
    const partners = lps.filter(lp => lp.kind === 'partner').map(lp => lp.partner);
    part(`pairs:${i}:${st.lone[i]}:${partners.join(',')}`, sub => {
      dirs.forEach((dir, q) => {
        const single = st.lone[i] % 2 === 1 && q === dirs.length - 1;
        const sp = single ? dotSprite() : pairSprite();
        // while a move plays a pair that is coming or going is drawn as faint as it is partial
        const pairs = st.lone[i] / 2, whole = Math.floor(pairs + 1e-9);
        if (movie && q >= whole && pairs - whole > 1e-9) { sp.material.opacity = Math.max(0.15, pairs - whole); sp.userData.baseOpacity = sp.material.opacity; }
        sub.add(sp);
        tracked.push({ sprite: sp, index: i, kind: 'pair', dir, angle: null, group: g });
      });
    });
  });

  // charges
  // While a move plays the charges are partial: a part of a charge shows as delta plus or delta
  // minus, fainter the smaller it is, with no number on it (the hover note has the number); within a
  // tenth of a whole charge it shows as the whole charge.
  sym.forEach((el, i) => {
    const q = st.charge[i];
    if (Math.abs(q) < (movie ? 0.08 : 1e-9)) return;
    const near = Math.round(q);
    const whole = !movie || (near !== 0 && Math.abs(q - near) < 0.1);
    const text = whole ? chargeMark(movie ? near : q) : `δ${q < 0 ? MINUS : '+'}`;
    const strength = whole ? 1 : 0.35 + 0.65 * Math.min(1, Math.abs(q) / 0.9);
    part(`charge:${i}:${text}:${strength.toFixed(2)}`, sub => {
      const sp = chargeSprite(text);
      sp.material.opacity = strength; sp.userData.baseOpacity = strength;
      sub.add(sp);
      tracked.push({ sprite: sp, index: i, kind: 'charge', angle: null, goal: null, group: g });
    });
  });
  return g;
}

// Directions in space for an atom's dot pairs. A pair pointing at a partner keeps that direction, and
// so does a pair in a slot that stands out of the atom's plane (the two pairs of water). Every other
// pair of a flat atom is drawn in its plane, as dots on paper are: the pairs are spread evenly over
// the gaps between the bonds, a wider gap taking more, so a terminal chlorine shows its three pairs
// at the top, the far side and the bottom, and the pair of a pyrrole or amide nitrogen, whose three
// gaps are equal, goes in the gap that faces away from the rest of the molecule. An atom with no
// plane to go by spreads its remaining pairs away from its bonds.
function drawingDirections(i, lps, frame, pos, centre) {
  const dirs = lps.map(lp => new THREE.Vector3(...lp.dir));
  const bonds = current.geo.groups[i].nbrs.map(nb => pos[nb.j].clone().sub(pos[i]).normalize());
  const pO = frame.pOrbital[i];
  if (!pO || !bonds.length) {
    const kept = dirs.filter((d, q) => lps[q].kind === 'slot' || lps[q].kind === 'partner');
    const need = dirs.length - kept.length;
    if (!need) return kept;
    let away = bonds.reduce((s, v) => s.add(v), new THREE.Vector3()).negate();
    if (away.length() < 1e-6) away.set(0, 1, 0);
    away.normalize();
    const p1 = Math.abs(away.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
    p1.addScaledVector(away, -p1.dot(away)).normalize();
    const p2 = away.clone().cross(p1);
    const extra = [];
    for (let q = 0; q < need; q++) {
      const ang = (2 * Math.PI * q) / need;
      extra.push(away.clone().addScaledVector(p1, 0.7 * Math.cos(ang)).addScaledVector(p2, 0.7 * Math.sin(ang)).normalize());
    }
    return [...kept, ...extra];
  }
  const nrm = new THREE.Vector3(...pO);
  const inPlane = v => Math.abs(v.dot(nrm)) < 0.5;
  const ex = bonds[0].clone().addScaledVector(nrm, -bonds[0].dot(nrm)).normalize();
  const ey = nrm.clone().cross(ex);
  const angleOf = v => Math.atan2(v.dot(ey), v.dot(ex));
  const held = dirs.filter((d, q) => lps[q].kind === 'partner' || (lps[q].kind === 'slot' && !inPlane(d)));
  const loose = dirs.length - held.length;
  // the gaps between the bonds and the held pairs in the plane, each with how far it faces outward
  const outward = angleOf(pos[i].clone().sub(centre));
  const taken = [...bonds, ...held.filter(inPlane)].map(angleOf).sort((a, b) => a - b);
  const gaps = taken.map((a, k) => {
    const width = (k + 1 < taken.length ? taken[k + 1] : taken[0] + TWO_PI) - a;
    return { start: a, width, count: 0, out: Math.cos(a + width / 2 - outward) };
  });
  for (let q = 0; q < loose; q++) {
    const room = g => g.width / (g.count + 1);
    const most = Math.max(...gaps.map(room));
    const best = gaps.filter(g => room(g) > most - 10 * DEG).sort((a, b) => b.out - a.out)[0];
    best.count++;
  }
  const placed = [];
  for (const g of gaps) for (let q = 0; q < g.count; q++) {
    const a = g.start + (g.width * (q + 1)) / (g.count + 1);
    placed.push(ex.clone().multiplyScalar(Math.cos(a)).addScaledVector(ey, Math.sin(a)));
  }
  return [...held, ...placed];
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
  const name = `${el} ${i + 1}`;   // the notes about arrows name atoms this way
  if (node.fixed[i]) return `<b>${name}</b> ${charge}<br>a loose ion, no bonds${lone ? `, ${lone} lone electrons` : ''}`;
  const ceiling = ELEMENTS[el].ceiling;
  const count = st.count[i];
  const num = v => (Math.abs(v - Math.round(v)) < 1e-9 ? String(Math.round(v)) : v.toFixed(2));   // partial during a move
  const octet = count >= ceiling - 1e-9 ? `counts ${num(count)} electrons, full` : `counts ${num(count)} electrons, ${num(ceiling - count)} short`;
  const lones = lone === 0 ? 'no lone electrons' : lone === 1 ? '1 lone electron' : `${num(lone)} lone electrons`;
  return `<b>${name}</b> ${charge}<br>${octet}, ${lones}`;
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
// The reaction scenes (and the SN2 transition state) in a short list; compounds come from the search.
const sel = document.getElementById('case');
const ENTRIES = [
  ...REACTIONS.map(rx => ({ kind: 'reaction', name: rx.name, smiles: rx.smiles, rx, group: 'reactions' })),
  ...CASES.filter(c => c.before).map(c => ({ kind: 'case', name: c.name, c, group: 'transitions' })),
];
{
  const o = document.createElement('option'); o.value = ''; o.textContent = 'reaction scenes…'; sel.appendChild(o);
  for (const groupName of ['reactions', 'transitions']) {
    const og = document.createElement('optgroup');
    og.label = groupName;
    ENTRIES.forEach((e, i) => { if (e.group !== groupName) return; const op = document.createElement('option'); op.value = i; op.textContent = e.name; og.appendChild(op); });
    sel.appendChild(og);
  }
}
const errBox = document.getElementById('error');
const subBox = document.getElementById('sub');
function tryLoad(entry, animate = true) {
  try { load(entry, animate); errBox.textContent = ''; return true; }
  catch (e) { errBox.textContent = e.message; return false; }
}
sel.addEventListener('change', () => {
  if (sel.value === '') return;
  const entry = ENTRIES[Number(sel.value)];
  sel.value = ''; sel.blur();
  if (tryLoad(entry)) { onScreen(entry.name); subBox.textContent = entry.kind === 'reaction' ? 'reagents placed together' : 'before and after at equal weight'; }
});

// ---------- the search: a name or a SMILES ----------
// The names (web/names.json.gz, built by tools/build-names.mjs from PubChem's names) arrive after
// the page; until then, and for anything not among them, the box reads SMILES.
const searchBox = document.getElementById('search');
const hitsBox = document.getElementById('hits');
let names = null;            // the index, once loaded
let namesNote = 'the names are still loading; a SMILES works now';
let shownName = '';          // what the box says when it is not being typed in
let hits = [], hitOn = -1;
(async () => {
  try {
    // names.txt is web/names.json.gz written out as base64 by build.sh (the page host serves no .gz)
    const res = await fetch('./names.txt');
    if (!res.ok) throw new Error(String(res.status));
    const raw = (await res.text()).trim();
    let text = raw;
    if (raw[0] !== '{') {
      const bin = atob(raw);
      const bytes = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      text = await new Response(stream).text();
    }
    names = makeIndex(JSON.parse(text));
    namesNote = '';
    if (document.activeElement === searchBox) suggest();
    if (current && !shownName) nameOnScreen();
  } catch (e) {
    namesNote = 'the names could not be loaded; a SMILES still works';
  }
})();

// the box shows the name of what is on screen, or a grey line when it has none
function onScreen(name, placeholder = 'a name or a SMILES') {
  shownName = name || '';
  searchBox.value = shownName;
  searchBox.placeholder = placeholder;
}
// after an edit or a move: the name of the molecule now on screen if the list has it
function nameOnScreen(fallback) {
  const found = current && current.spec ? namesFor(names, current.spec) : null;
  if (found) onScreen(found.join(' + '));
  else onScreen('', fallback || 'a name or a SMILES');
}

function suggest() {
  const text = searchBox.value.trim();
  hits = names && text ? search(names, text, 8) : [];
  hitOn = -1;
  drawHits(text);
}
function drawHits(text) {
  hitsBox.innerHTML = '';
  if (!text) { hitsBox.classList.remove('open'); return; }
  if (!hits.length) {
    const d = document.createElement('div');
    d.className = 'none';
    d.textContent = namesNote || 'no name matches; enter reads it as a SMILES';
    hitsBox.appendChild(d);
  }
  hits.forEach((h, q) => {
    const row = document.createElement('div');
    row.className = 'hit' + (q === hitOn ? ' on' : '') + (h.entry.smiles ? '' : ' off');
    const nm = document.createElement('div');
    nm.className = 'nm';
    nm.textContent = h.entry.name;
    if (h.matched !== h.entry.name || !h.entry.smiles) {
      const why = document.createElement('span');
      why.className = 'why';
      why.textContent = h.entry.smiles ? `also called ${h.matched}` : `not shown: ${reasonText(h.entry.key)}`;
      nm.appendChild(why);
    }
    row.appendChild(nm);
    if (h.entry.smiles) {
      const add = document.createElement('span');
      add.className = 'add';
      add.textContent = 'add';
      add.title = 'add it beside what is on screen';
      add.addEventListener('mousedown', e => { e.preventDefault(); e.stopPropagation(); pickHit(q, true); });
      row.appendChild(add);
    }
    row.addEventListener('mousedown', e => { e.preventDefault(); pickHit(q, false); });
    hitsBox.appendChild(row);
  });
  hitsBox.classList.add('open');
}
function closeHits() { hits = []; hitOn = -1; hitsBox.classList.remove('open'); hitsBox.innerHTML = ''; }

// load a compound, or add it beside what is on screen
function loadCompound(name, smiles, beside) {
  const balanced = balanceSmiles(smiles);
  const ok = beside ? addBeside(balanced.smiles) : tryLoad({ kind: 'smiles', name, smiles: balanced.smiles });
  if (!ok) return false;
  if (beside) nameOnScreen('two compounds on screen');
  else onScreen(name);
  const added = balanced.added ? `, with ${balanced.added.count > 1 ? balanced.added.count + ' ' : ''}${balanced.added.ion} beside it to balance the charge` : '';
  subBox.textContent = beside ? `${name} added beside${added}` : added ? added.slice(2) : '';
  return true;
}
function pickHit(q, beside) {
  const h = hits[q];
  if (!h) return;
  closeHits();
  if (!h.entry.smiles) { errBox.textContent = `${h.entry.name} cannot be shown yet: ${reasonText(h.entry.key)}`; searchBox.select(); return; }
  if (loadCompound(h.entry.name, h.entry.smiles, beside)) searchBox.blur();
}
// Enter: the highlighted suggestion; else a name equal to what was typed; else the text as a SMILES;
// else the first suggestion
function enter(beside) {
  const text = searchBox.value.trim();
  if (!text) return;
  if (hitOn >= 0) { pickHit(hitOn, beside); return; }
  if (names) {
    const found = search(names, text, 8);
    if (found.length && found[0].exact) { hits = found; pickHit(0, beside); return; }
    hits = found;
  }
  let smilesOk = true;
  try { nodeFromSmiles(text); } catch (e) { smilesOk = false; }
  if (smilesOk) {
    closeHits();
    if (loadCompound(text, text, beside)) {
      if (!beside) nameOnScreen('typed as a SMILES');
      searchBox.blur();
    } else searchBox.select();
    return;
  }
  if (hits.length) { pickHit(0, beside); return; }
  closeHits();
  errBox.textContent = names ? `no compound called "${text}" in the list, and it does not read as a SMILES` : `"${text}" does not read as a SMILES, and ${namesNote}`;
  searchBox.select();
}
searchBox.addEventListener('focus', () => { searchBox.select(); });
searchBox.addEventListener('input', () => { errBox.textContent = ''; suggest(); });
searchBox.addEventListener('keydown', e => {
  if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
    if (!hits.length) return;
    e.preventDefault();
    hitOn = e.key === 'ArrowDown' ? Math.min(hits.length - 1, hitOn + 1) : Math.max(-1, hitOn - 1);
    drawHits(searchBox.value.trim());
  } else if (e.key === 'Escape') {
    closeHits(); searchBox.value = shownName; searchBox.blur();
  } else if (e.key === 'Enter') {
    e.preventDefault();
    enter(e.shiftKey);
  }
});
searchBox.addEventListener('blur', () => { closeHits(); if (!errBox.textContent) searchBox.value = shownName; });
document.getElementById('search-form').addEventListener('submit', e => e.preventDefault());

// a second compound beside what is on screen, so a reagent can be brought in
function addBeside(smiles) {
  if (!current || !current.spec || movie) { errBox.textContent = 'there is nothing to add it to; it is loaded on its own instead'; return tryLoad({ kind: 'smiles', name: smiles, smiles }); }
  try {
    const other = nodeFromSmiles(smiles);
    const c = combine(current.spec, other);
    // the newcomer starts to the right of everything as seen on screen, in its own relaxed shape,
    // level with the middle of what is there
    const own = relax(buildNode(other), enumerate(buildNode(other)));
    camera.matrixWorld.extractBasis(camRight, camUp, camBack);
    const along = p => camRight.dot(new THREE.Vector3(...p));
    const oldMid = current.geo.pos.reduce((a, p) => a.add(new THREE.Vector3(...p)), new THREE.Vector3()).divideScalar(current.geo.pos.length);
    const newMid = own.pos.reduce((a, p) => a.add(new THREE.Vector3(...p)), new THREE.Vector3()).divideScalar(own.pos.length);
    const gap = Math.max(...current.geo.pos.map(along)) - Math.min(...own.pos.map(along)) + 2.8;
    const move3 = oldMid.clone().sub(newMid);
    move3.addScaledVector(camRight, gap - camRight.dot(move3));
    const start = [...current.geo.pos.map(p => [...p]), ...own.pos.map(p => [p[0] + move3.x, p[1] + move3.y, p[2] + move3.z])];
    loadSpec(c.spec, start);
    current.reaction = null;
    showHint();
    errBox.textContent = '';
    moveCamera(frameFor(current.node, current.geo, null, { dir: camera.position.clone().sub(controls.target).normalize(), up: camera.up }), true);
    return true;
  } catch (err) { errBox.textContent = err.message; return false; }
}

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
  if (arrows) arrows.dropIfStale();
  lastSwitch = performance.now();
  showStructure(true);
}
function togglePause() {
  if (movie) {
    // space holds the film where it is, and lets it go on
    const now = performance.now();
    if (movie.pausedAt === null) movie.pausedAt = now;
    else { movie.start += now - movie.pausedAt; movie.pausedAt = null; }
    updatePanel();
    return;
  }
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
  if (movie) {
    const made = movie.out.made.length, broken = movie.out.broken.length;
    // the products are not named until the film has shown them
    const lines = [movie.pausedAt === null ? '<b>the move</b> is playing (space holds it)' : '<b>the move</b> is held (space goes on)', `${made} bond${made === 1 ? '' : 's'} made, ${broken} broken`];
    document.getElementById('status').innerHTML = lines.map(l => `<div>${l}</div>`).join('');
    document.getElementById('steps').style.visibility = 'hidden';
    return;
  }
  const { r, families } = current;
  const n = r.structures.length;
  const stepping = showResonance && families.length > 1;
  document.getElementById('steps').style.visibility = stepping ? 'visible' : 'hidden';
  // the engine shares out the pool; a loose anion keeps its own electrons outside it, and they count
  // too (a chloride's eight: the SN2 products showed 14 where the reagents showed 22)
  const own = current.node.fixed.reduce((sum, f) => sum + (f ? f.lone : 0), 0);
  const lines = [`<b>${r.pool + own}</b> valence electrons` + (current.spec && current.spec.charge ? ` · charge <b>${(current.spec.charge > 0 ? '+' : MINUS) + Math.abs(current.spec.charge)}</b>` : ''), `<b>${n}</b> contributing structure${n === 1 ? '' : 's'}`];
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
  if (movie) stepMovie(t);
  if (current) {
    if (fade) {
      const f = Math.min(1, (t - fade.start) / FADE);
      const e = f * f * (3 - 2 * f);
      fade.going.forEach(c => setOpacity(c, 1 - e));
      fade.coming.forEach(c => setOpacity(c, e));
      if (f >= 1) finishFade();
    }
    const family = current.families[familyIdx];
    const held = (arrows && arrows.count()) || movie;   // arrows drawn: the picture stays put
    if (showResonance && !paused && family.length > 1 && !fade && !held && t - lastSwitch > DWELL) {
      memberIdx = (memberIdx + 1) % family.length;
      lastSwitch = t;
      showStructure(true);
    }
    placeSprites(dt);
    updateTip();
  }
  if (arrows) arrows.update();
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
    // the box names what is on screen now, if the list has it
    nameOnScreen(current ? 'edited molecule' : 'nothing on screen');
    subBox.textContent = '';
    errBox.textContent = '';
    // say so when the octet leaves no room for a bond the edit asked for (a fifth neighbour on C)
    cutNote = '';
    if (current && current.spec) {
      const st = current.r.structures[current.main];
      const cut = current.r.edges.filter(([i, j], k) => st.orders[k] === 0 && !current.node.fixed[i] && !current.node.fixed[j]);
      const nm = i => `${current.node.atoms[i].el}${i}`;
      if (cut.length) cutNote = `No room in the octet for ${cut.length === 1 ? 'the bond' : 'the bonds'} ${cut.map(([i, j]) => `${nm(i)}-${nm(j)}`).join(', ')}: drawn hashed, as no bond.`;
    }
  },
  setNote: text => { document.getElementById('note').textContent = cutNote ? `${cutNote} ${text}` : text; },
});
document.getElementById('fill').addEventListener('click', e => { editor.fill(); e.target.blur(); });

// ---------- moves: curved arrows, and the short film from before to after ----------
arrows = attachArrows({
  scene, camera, domElement: renderer.domElement,
  controls: () => controls,
  getCurrent: () => current,
  getStructure: () => shownStructure(),
  pairSprites: () => tracked.filter(t => t.kind === 'pair').map(t => ({ sprite: t.sprite, index: t.index })),
  bondMeshes: () => { const out = []; structureGroup.traverse(o => { if (o.isMesh && o.userData.edge !== undefined && o.visible) out.push(o); }); return out; },
  atomMeshes: () => atomGroup.children.filter(o => o.userData.atom !== undefined),
  isAllowed: () => !!current && !!current.spec && !movie && !move,
  setNote: text => { document.getElementById('note').textContent = text; },
  onCommit: makeMove,
  onChange: count => { document.getElementById('go').disabled = !count; document.getElementById('undo').disabled = !count; },
});
document.getElementById('go').addEventListener('click', e => { arrows.commit(); e.target.blur(); });
document.getElementById('undo').addEventListener('click', e => { arrows.clear(); e.target.blur(); });
document.getElementById('show').addEventListener('click', e => {
  e.target.blur();
  if (!current || !current.reaction || movie) return;
  // the textbook arrows are drawn on the structure where every tail exists; switch to it if needed
  const list = resolveArrows(current.spec, current.reaction.arrows);
  const fits = st => list.every(a => ('pair' in a.from ? st.lone[a.from.pair] >= 2 : (() => { const k = current.r.edges.findIndex(e => (e[0] === a.from.bond[0] && e[1] === a.from.bond[1]) || (e[0] === a.from.bond[1] && e[1] === a.from.bond[0])); return k >= 0 && st.orders[k] >= 1; })()));
  if (!fits(shownStructure())) {
    const q = current.r.structures.findIndex(fits);
    if (q >= 0) {
      const fi = current.families.findIndex(f => f.includes(q));
      showResonance = true; resonanceBox.checked = true;
      familyIdx = fi; memberIdx = current.families[fi].indexOf(q); lastSwitch = performance.now();
      showStructure(true);
    }
  }
  arrows.set(list);
});
function showHint() {
  document.getElementById('show').style.display = current && current.reaction ? '' : 'none';
}
arrows.note();

function makeMove(list) {
  const { node, r } = current;
  const out = commit(node, shownStructure(), list, r);
  if (out.kind === 'refused') { arrows.note(`not allowed: ${out.reason}`); arrows.flash(out.atoms); return; }
  if (out.kind === 'resonance') {
    const s = out.applied.structure;
    const q = r.structures.findIndex(t => t.orders.every((o, k) => o === s.orders[k]) && t.lone.every((v, i) => v === s.lone[i]) && t.charge.every((c, i) => c === s.charge[i]));
    arrows.clear();
    if (q >= 0) {
      const fi = current.families.findIndex(f => f.includes(q));
      showResonance = true; resonanceBox.checked = true;
      familyIdx = fi; memberIdx = current.families[fi].indexOf(q); lastSwitch = performance.now();
      showStructure(true);
    }
    arrows.note('no bond made or broken: that is resonance, another picture of the same molecule, shown now');
    return;
  }
  playMove(out);
}

// The move as a short film (item 15): the weights blend from all-before to all-after, so every bond
// order and charge slides from its old value to its new one; the shape follows, relaxed step by
// step; at the end the products are loaded as the molecule on screen.
function playMove(out) {
  const film = filmSteps(current.node, out, current.geo, 21);
  const { node, steps, roles, endPos } = film;
  arrows.clear();
  finishFade();
  const products = pieceFormulas(out.spec, out.after);
  // the film begins on the picture the arrows were drawn on and ends on the picture the products will
  // show, so nothing jumps when it starts or when the products are loaded (a ring of five structures
  // turned from one picture to their average the moment the film began; moves review, 7 Nov 2026)
  const lift = (st, edges) => {
    const at = new Map(edges.map((e, k) => [`${e[0]}-${e[1]}`, k]));
    return { ...st, orders: out.hybrid.edges.map(e => { const k = at.get(`${e[0]}-${e[1]}`); return k === undefined ? 0 : st.orders[k]; }) };
  };
  const { main, families } = mainAndFamilies(out.after, out.rAfter);
  const first = lift(shownStructure(), current.r.edges);
  const last = lift(out.rAfter.structures[showResonance ? families[0][0] : main], out.rAfter.edges);
  movie = { node, steps, roles, first, last, start: performance.now(), out, products, shown: null, endPos, pausedAt: null };
  // the camera draws back at the start, if it must, so that nothing leaves the picture during the film
  const every = { atoms: steps.flatMap(() => node.atoms) };
  const everyGeo = { pos: steps.flatMap(st => st.geo.pos) };
  if (outOfView(everyGeo, every)) moveCamera(frameFor(every, everyGeo, null, { dir: camera.position.clone().sub(controls.target).normalize(), up: camera.up }), true);
  current = { ...current, node, r: steps[0].r, geo: steps[0].geo };
  buildAtoms();
  stepMovie(performance.now());
  updatePanel();
  arrows.note('watch the bonds, pairs and charges change; hover an atom for its numbers');
}
function movieMix() { return movie && movie.shown ? movie.mix : 0; }
function stepMovie(t) {
  if (movie.pausedAt !== null) return;
  const f = Math.min(1, (t - movie.start) / MOVIE);
  const last = movie.steps.length - 1;
  // a fifth of the time for the approach, then the blend, then the products settling, each eased
  const ease = v => v * v * (3 - 2 * v);
  const x = f < 0.2 ? ease(f / 0.2) : f < 0.85 ? 1 + ease((f - 0.2) / 0.65) * (last - 2) : (last - 1) + ease((f - 0.85) / 0.15);
  const k = Math.min(last - 1, Math.floor(x));
  const u = x - k;
  const a = movie.steps[k], b = movie.steps[k + 1];
  const pos = a.geo.pos.map((p, i) => p.map((v, q) => v + (b.geo.pos[i][q] - v) * u));
  const mix = a.mix + (b.mix - a.mix) * u;
  // the drawing: from the picture on screen to the average of the reagents' structures during the
  // approach, then the weights blend from before to after, then from that average to the picture the
  // products will show during the settle
  const A = movie.steps[0].r.avg, B = movie.steps[last].r.avg;
  const [P, Q, w] = f < 0.2 ? [movie.first, A, ease(f / 0.2)] : f < 0.85 ? [A, B, mix] : [B, movie.last, ease((f - 0.85) / 0.15)];
  const blend = key => P[key].map((v, i) => v + (Q[key][i] - v) * w);
  const orders = blend('orders'), lone = blend('lone');
  // the electrons an atom counts, from the pairs and bonds drawn (a leaving bromide counts its eight
  // all the way, not the nothing a loose ion counts in the engine's books)
  const count = lone.map(v => v);
  movie.node.edges.forEach(([i, j], k) => { count[i] += 2 * orders[k]; count[j] += 2 * orders[k]; });
  movie.shown = { orders, lone, charge: blend('charge'), count, cuts: 0, weight: 1 };
  movie.mix = mix;
  current.r = a.r;
  current.geo = { ...a.geo, pos };
  atomGroup.children.forEach(o => { if (o.userData.atom !== undefined) o.position.set(...pos[o.userData.atom]); });
  swapStructure(buildStructure(), false);
  if (f >= 1) {
    const out = movie.out, endPos = movie.endPos;
    movie = null;
    loadSpec(out.spec, endPos, true);
    // the products may have drawn together or come apart: fit the picture to them, gently
    moveCamera(frameFor(current.node, current.geo, null, { dir: camera.position.clone().sub(controls.target).normalize(), up: camera.up }), true);
    nameOnScreen('after the move');
    subBox.textContent = `after the move: ${pieceFormulas(out.spec, out.after).join(', ')}`;
    errBox.textContent = '';
    current.reaction = null; showHint();
    const made = out.made.length, broken = out.broken.length;
    arrows.note(`done: ${made} bond${made === 1 ? '' : 's'} made, ${broken} broken. Drag from a pair of dots or a bond to draw the next arrow`);
  }
}
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

// the page opens on SF6
tryLoad({ kind: 'smiles', name: 'sulfur hexafluoride', smiles: 'FS(F)(F)(F)(F)F' }, false);
onScreen('sulfur hexafluoride');
requestAnimationFrame(animate);

// a handle for scripted tests, so a script can find an atom on the screen and click it
window.electronFlow = {
  THREE, camera, canvas: renderer.domElement,
  positions: () => (current ? current.geo.pos : []),
  pairs: () => tracked.filter(t => t.kind === 'pair').map(t => ({ index: t.index, pos: t.sprite.position.toArray() })),
  arrows: () => (arrows ? arrows.arrows() : []),
  playing: () => !!movie,
  smiles: () => (current && current.spec ? toSmiles(current.node, current.r.structures[current.main]) : ''),
  // hold the film at a fraction of its length (0 to 1), for frame-by-frame checks; null lets it go on
  holdFilm: f => {
    if (!movie) return false;
    const now = performance.now();
    if (f === null) { if (movie.pausedAt !== null) { movie.start += now - movie.pausedAt; movie.pausedAt = null; } return true; }
    movie.pausedAt = null; movie.start = now - f * MOVIE; stepMovie(now);
    if (movie) movie.pausedAt = now;
    return true;
  },
  namesReady: () => !!names,
};
