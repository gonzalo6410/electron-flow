// Electron Flow viewer. Renders what the core computes; nothing chemical is decided here.
import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';
import { CASES, runCase } from './core/cases.js';
import { relax, lonePairDirections } from './core/geometry.js';
import { ELEMENTS } from './core/elements.js';

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
scene.add(new THREE.HemisphereLight(0xffffff, 0x888888, 1.1));
const dir = new THREE.DirectionalLight(0xffffff, 1.2);
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

const sphereGeo = new THREE.SphereGeometry(1, 32, 24);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 16, 1);
const matCache = new Map();
function mat(color, opacity = 1) {
  const key = `${color}-${opacity}`;
  if (!matCache.has(key)) matCache.set(key, new THREE.MeshLambertMaterial({ color, transparent: opacity < 1, opacity }));
  return matCache.get(key);
}
const atomRadius = el => 0.15 + 0.15 * ELEMENTS[el].radius;

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
  const seg = 0.16, gap = 0.12;
  const n = Math.max(1, Math.floor((L + gap) / (seg + gap)));
  const total = n * seg + (n - 1) * gap;
  const start = (L - total) / 2;
  const u = d.clone().normalize();
  for (let k = 0; k < n; k++) {
    const s = start + k * (seg + gap);
    const p = va.clone().add(u.clone().multiplyScalar(s));
    const q = va.clone().add(u.clone().multiplyScalar(s + seg));
    group.add(cylinder(p.toArray(), q.toArray(), radius, material));
  }
}

function label(text, pos, size = 0.5, color = '#000') {
  const c = document.createElement('canvas');
  c.width = 256; c.height = 128;
  const g = c.getContext('2d');
  g.font = '56px Helvetica, Arial, sans-serif';
  g.textAlign = 'center'; g.textBaseline = 'middle';
  g.lineWidth = 10; g.strokeStyle = 'rgba(255,255,255,0.85)'; g.lineJoin = 'round';
  g.strokeText(text, 128, 64);
  g.fillStyle = color;
  g.fillText(text, 128, 64);
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false }));
  sp.position.set(...pos);
  sp.scale.set(size * 2, size, 1);
  return sp;
}

// ---------- state ----------
let current = null;   // {c, node, r, geo, lpDirsByStructure}
let mode = 'cycle';
let structIdx = 0;
let lastSwitch = 0;
const molecule = new THREE.Group();
scene.add(molecule);

function load(idx) {
  const c = CASES[idx];
  const { node, r } = runCase(c);
  const geo = relax(node, r);
  current = { c, node, r, geo };
  structIdx = 0;
  lastSwitch = performance.now();
  // frame the molecule: camera a little above the mean plane, far enough to see everything
  const extent = Math.max(1.5, ...geo.pos.map(p => Math.hypot(...p)));
  const dist = extent / Math.tan((camera.fov / 2) * Math.PI / 180) * 1.5 + 1;
  camera.position.set(dist * 0.38, dist * 0.42, dist * 0.82); // off every axis, so nothing lines up behind anything
  controls.target.set(0, 0, 0);
  controls.update();
  document.getElementById('rules').textContent =
    'Rules: octet never exceeded; contributing structures are those with the fewest octet holes; sigma bonds break only where the octet forces it; equal weights; geometry from electron groups and springs.';
  draw();
}

// what to draw: a single structure (cycle) or the average
function picture() {
  const { r } = current;
  if (mode === 'cycle') {
    const st = r.structures[structIdx];
    return { orders: st.orders, charge: st.charge, lone: st.lone, count: st.count, single: true };
  }
  return { orders: r.avg.orders, charge: r.avg.charge, lone: r.avg.lone, count: r.avg.count, single: false };
}

function draw() {
  molecule.clear();
  const { node, r, geo } = current;
  const pos = geo.pos;
  const pic = picture();
  const sym = node.atoms.map(a => a.el);

  // atoms
  sym.forEach((el, i) => {
    const m = new THREE.Mesh(sphereGeo, mat(ELEMENTS[el].color));
    const rad = atomRadius(el);
    m.position.set(...pos[i]);
    m.scale.setScalar(rad);
    molecule.add(m);
    const q = pic.charge[i];
    const text = el + (Math.abs(q) > 1e-9 ? ' ' + frac(q) : '');
    molecule.add(label(text, [pos[i][0], pos[i][1] + rad + 0.22, pos[i][2]], 0.36));
  });

  // bonds
  r.edges.forEach(([i, j], k) => {
    const o = pic.orders[k];
    const loose = node.els[i].ionic || node.els[j].ionic;
    if (o <= 1e-9) {
      if (loose) dashed(pos[i], pos[j], 0.025, mat(0xbbbbbb), molecule);
      return;
    }
    const full = Math.floor(o + 1e-9), fr = o - full;
    const lines = full + (fr > 1e-9 ? 1 : 0);
    // the extra lines go on the side where the neighbours of the two atoms sit (inside a ring),
    // and in the plane of the bond and its first neighbour when the neighbours balance (ethylene)
    const a = new THREE.Vector3(...pos[i]), b = new THREE.Vector3(...pos[j]);
    const bd = b.clone().sub(a).normalize();
    const mid = a.clone().add(b).multiplyScalar(0.5);
    const perpOf = v => v.clone().sub(bd.clone().multiplyScalar(v.dot(bd)));
    let sum = new THREE.Vector3(), first = null;
    for (const [p, q] of r.edges) {
      const other = p === i ? q : q === i ? p : p === j ? q : q === j ? p : -1;
      if (other < 0 || other === i || other === j) continue;
      const v = perpOf(new THREE.Vector3(...pos[other]).sub(mid));
      sum.add(v);
      if (!first && v.length() > 1e-6) first = v;
    }
    let ref = sum.length() > 0.3 ? sum : first;
    if (!ref) ref = perpOf(Math.abs(bd.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0));
    const off = ref.clone().normalize().multiplyScalar(0.17);
    // shorten to the sphere surfaces
    const ra = atomRadius(sym[i]) * 0.9, rb = atomRadius(sym[j]) * 0.9;
    const a2 = a.clone().add(bd.clone().multiplyScalar(ra)), b2 = b.clone().sub(bd.clone().multiplyScalar(rb));
    for (let q = 0; q < lines; q++) {
      const shift = off.clone().multiplyScalar(q - (lines - 1) / 2);
      const pa = a2.clone().add(shift).toArray(), pb = b2.clone().add(shift).toArray();
      const partial = q === lines - 1 && fr > 1e-9;
      if (partial) dashed(pa, pb, 0.055, mat(0x333333), molecule);
      else molecule.add(cylinder(pa, pb, 0.07, mat(0x333333)));
    }
  });

  // lone pairs: whole pairs as two dots, a fractional pair faded
  const lpCount = sym.map((el, i) => (node.els[i].ionic ? 0 : Math.ceil(pic.lone[i] / 2 - 1e-9)));
  const dirs = lonePairDirections(pos, geo.groups, lpCount);
  sym.forEach((el, i) => {
    const rad = atomRadius(el);
    const whole = Math.floor(pic.lone[i] / 2 + 1e-9);
    dirs[i].forEach((d, q) => {
      const fractional = q >= whole;
      const opacity = fractional ? Math.max(0.15, (pic.lone[i] / 2) - whole) : 1;
      const base = new THREE.Vector3(...pos[i]).add(new THREE.Vector3(...d).multiplyScalar(rad + 0.14));
      const perp = new THREE.Vector3(...d).cross(Math.abs(d[0]) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0)).normalize().multiplyScalar(0.07);
      for (const s of [-1, 1]) {
        const m = new THREE.Mesh(sphereGeo, mat(0x000000, opacity));
        m.position.copy(base).add(perp.clone().multiplyScalar(s));
        m.scale.setScalar(0.05);
        molecule.add(m);
      }
    });
  });

  // readout
  const st = document.getElementById('status');
  const n = r.structures.length;
  st.innerHTML = `${r.pool} valence electrons, ${n} contributing structure${n === 1 ? '' : 's'}, ${r.holes / 2} octet hole${r.holes === 2 ? '' : 's'}` +
    (r.transitionEdges ? '<br>transition: before and after each weighted 1/2' : '') +
    (mode === 'cycle' && n > 1 ? `<br>showing structure ${structIdx + 1} of ${n}` : '') +
    (mode === 'average' ? '<br>weighted average of all contributing structures' : '');

  const atomsT = document.getElementById('atoms');
  atomsT.innerHTML = '<tr class="small"><td></td><td class="num">charge</td><td class="num">count</td><td class="num">lone e</td></tr>' +
    sym.map((el, i) => `<tr><td>${el}${i}</td><td class="num">${frac(pic.charge[i])}</td><td class="num">${node.els[i].ionic ? '' : frac(pic.count[i], false)}</td><td class="num">${node.els[i].ionic ? '' : frac(pic.lone[i], false)}</td></tr>`).join('');
  const bondsT = document.getElementById('bonds');
  bondsT.innerHTML = r.edges.map(([i, j], k) => {
    const o = pic.orders[k];
    const loose = node.els[i].ionic || node.els[j].ionic;
    return `<tr><td>${sym[i]}${i} ${sym[j]}${j}</td><td class="num">${loose && o <= 1e-9 ? 'ion pair' : frac(o, false)}</td></tr>`;
  }).join('');
}

// ---------- controls ----------
const sel = document.getElementById('case');
CASES.forEach((c, i) => { const o = document.createElement('option'); o.value = i; o.textContent = c.name; sel.appendChild(o); });
sel.addEventListener('change', () => load(Number(sel.value)));
const bCycle = document.getElementById('btn-cycle'), bAvg = document.getElementById('btn-avg');
bCycle.addEventListener('click', () => { mode = 'cycle'; bCycle.classList.add('on'); bAvg.classList.remove('on'); draw(); });
bAvg.addEventListener('click', () => { mode = 'average'; bAvg.classList.add('on'); bCycle.classList.remove('on'); draw(); });

const DWELL = 1000; // ms per structure; equal weights, so equal time
function animate(t) {
  requestAnimationFrame(animate);
  if (current && mode === 'cycle' && current.r.structures.length > 1 && t - lastSwitch > DWELL) {
    structIdx = (structIdx + 1) % current.r.structures.length;
    lastSwitch = t;
    draw();
  }
  controls.update();
  renderer.render(scene, camera);
}
load(0);
requestAnimationFrame(animate);
