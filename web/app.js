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
let structIdx = 0;
let lastSwitch = 0;
let fade = null;      // {from, to, start}
const atomGroup = new THREE.Group();   // spheres and symbols, fixed per case
let structureGroup = new THREE.Group(); // bonds, charges, lone pairs, per structure
scene.add(atomGroup);
scene.add(structureGroup);
const faceCamera = []; // {sprite, index, kind} sprites repositioned every frame

function load(idx) {
  const c = CASES[idx];
  const { node, r } = runCase(c);
  const geo = relax(node, r);
  const maxPairs = node.atoms.map((a, i) => node.els[i].ionic ? 0 : Math.max(...r.structures.map(st => Math.floor(st.lone[i] / 2)), Math.ceil(r.avg.lone[i] / 2 - 1e-9)));
  const lpDirs = lonePairDirections(geo.pos, geo.groups, maxPairs);
  current = { c, node, r, geo, lpDirs };
  structIdx = 0;
  lastSwitch = performance.now();
  fade = null;
  document.getElementById('rules').textContent =
    'Rules: the octet is never exceeded; contributing structures are the ones with the fewest octet holes; sigma bonds break only where the octet forces it; equal weights; geometry from electron groups and springs.';
  // frame the molecule
  const extent = Math.max(1.5, ...geo.pos.map(p => Math.hypot(...p)));
  const dist = extent / Math.tan((camera.fov / 2) * Math.PI / 180) * 1.5 + 1;
  camera.position.set(dist * 0.38, dist * 0.42, dist * 0.82);
  controls.target.set(0, 0, 0);
  controls.update();
  buildAtoms();
  swapStructure(buildStructure(picture()), false);
  readout();
}

function picture() {
  const { r } = current;
  if (mode === 'cycle') {
    const st = r.structures[structIdx];
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
  const g = new THREE.Group();
  const { node, r, geo, lpDirs } = current;
  const pos = geo.pos;
  const sym = node.atoms.map(a => a.el);
  const bondMat = solid(0x3a3a3a);
  const ghostMat = solid(0xd0d0d0);
  const ionMat = solid(0xc8c8c8);

  // bonds
  r.edges.forEach(([i, j], k) => {
    const o = pic.orders[k];
    const loose = node.els[i].ionic || node.els[j].ionic;
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
    const full = Math.floor(o + 1e-9), fr = o - full;
    const lines = full + (fr > 1e-9 ? 1 : 0);
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
    for (let q = 0; q < lines; q++) {
      const shift = off.clone().multiplyScalar(q - (lines - 1) / 2);
      const pa = a2.clone().add(shift).toArray(), pb = b2.clone().add(shift).toArray();
      const partial = q === lines - 1 && fr > 1e-9;
      if (partial) dashed(pa, pb, 0.055, bondMat, g);
      else g.add(cylinder(pa, pb, 0.07, bondMat));
    }
  });

  // lone pairs: whole pairs full, a fractional pair faded
  sym.forEach((el, i) => {
    if (node.els[i].ionic) return;
    const rad = atomRadius(el);
    const pairs = pic.lone[i] / 2;
    const whole = Math.floor(pairs + 1e-9);
    const shown = Math.ceil(pairs - 1e-9);
    lpDirs[i].slice(0, shown).forEach((d, q) => {
      const opacity = q < whole ? 1 : Math.max(0.2, pairs - whole);
      const sp = pairSprite(opacity);
      sp.userData.baseOpacity = opacity;
      sp.position.set(pos[i][0] + d[0] * (rad + 0.17), pos[i][1] + d[1] * (rad + 0.17), pos[i][2] + d[2] * (rad + 0.17));
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

function readout() {
  const { node, r } = current;
  const pic = picture();
  const sym = node.atoms.map(a => a.el);
  const n = r.structures.length;
  const st = document.getElementById('status');
  st.innerHTML = `<b>${r.pool}</b> valence electrons · <b>${n}</b> contributing structure${n === 1 ? '' : 's'} · <b>${r.holes / 2}</b> octet hole${r.holes === 2 ? '' : 's'}` +
    (r.transitionEdges ? '<br>transition: before and after each weighted 1/2' : '') +
    (mode === 'cycle' && n > 1 ? `<br>structure <b>${structIdx + 1}</b> of ${n}` : '') +
    (mode === 'average' ? '<br>weighted average of all contributing structures' : '');
  const label = i => `${sym[i]}<sub>${i}</sub>`;
  document.getElementById('atoms').innerHTML =
    '<tr><th></th><th>charge</th><th>count</th><th>lone e</th></tr>' +
    sym.map((el, i) => node.els[i].ionic
      ? `<tr><td class="sym">${label(i)}</td><td>${frac(pic.charge[i])}</td><td class="dim">ion</td><td class="dim"></td></tr>`
      : `<tr><td class="sym">${label(i)}</td><td>${frac(pic.charge[i])}</td><td>${frac(pic.count[i], false)}</td><td>${frac(pic.lone[i], false)}</td></tr>`).join('');
  document.getElementById('bonds').innerHTML =
    '<tr><th></th><th>order</th></tr>' +
    r.edges.map(([i, j], k) => {
      const o = pic.orders[k];
      const loose = node.els[i].ionic || node.els[j].ionic;
      return `<tr><td class="sym">${label(i)} – ${label(j)}</td><td${loose ? ' class="dim"' : ''}>${loose && o <= 1e-9 ? 'ion pair' : frac(o, false)}</td></tr>`;
    }).join('');
}

// ---------- controls ----------
const sel = document.getElementById('case');
CASES.forEach((c, i) => { const o = document.createElement('option'); o.value = i; o.textContent = c.name; sel.appendChild(o); });
sel.addEventListener('change', () => load(Number(sel.value)));
const bCycle = document.getElementById('btn-cycle'), bAvg = document.getElementById('btn-avg');
function setMode(m) {
  mode = m;
  bCycle.classList.toggle('on', m === 'cycle');
  bAvg.classList.toggle('on', m === 'average');
  lastSwitch = performance.now();
  swapStructure(buildStructure(picture()), true);
  readout();
}
bCycle.addEventListener('click', () => setMode('cycle'));
bAvg.addEventListener('click', () => setMode('average'));

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
    if (mode === 'cycle' && current.r.structures.length > 1 && !fade && t - lastSwitch > DWELL) {
      structIdx = (structIdx + 1) % current.r.structures.length;
      lastSwitch = t;
      swapStructure(buildStructure(picture()), true);
      readout();
    }
    placeSprites();
  }
  controls.update();
  renderer.render(scene, camera);
}
load(0);
requestAnimationFrame(animate);
