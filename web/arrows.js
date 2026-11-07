// Curved arrows on the page: drag from a pair of dots or from a bond to an atom or to a bond, and a
// curved arrow is drawn. Arrows wait until Enter (or the button) makes the move; Escape clears them;
// a click on an arrow removes it. What the arrows do is decided in core/moves.js; this file only
// draws them and turns the pointer into tails and heads. It talks to the viewer through a small api
// object (attachArrows), like the editor does.
import * as THREE from 'three';
import { arrowShape, tailShape, describeArrow } from './core/moves.js';

const COLOR = 0x1f4fa3;        // arrows in a dark blue, apart from the grey of bonds
const TUBE = 0.042;            // radius of the arrow's line
const HEAD_LEN = 0.26, HEAD_R = 0.1;
const key = (i, j) => `${Math.min(i, j)}-${Math.max(i, j)}`;

export function attachArrows(api) {
  // api: { scene, camera, domElement, controls(), getCurrent(), getStructure(), pairSprites(),
  //        bondMeshes(), atomMeshes(), isAllowed(), setNote(text), onCommit(arrows), onChange(count) }
  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  const group = new THREE.Group();
  api.scene.add(group);
  const material = new THREE.MeshLambertMaterial({ color: COLOR, transparent: true, opacity: 0.92 });
  const ghostMaterial = new THREE.MeshLambertMaterial({ color: COLOR, transparent: true, opacity: 0.45 });
  const pending = [];          // { arrow, dot }  dot: which of the atom's dot pairs the arrow leaves from
  let drag = null;             // { tail: { from, dot }, point: Vector3, head: arrow head | null }
  const flashes = [];          // red halos on the atoms a refusal names

  // ---------- picking ----------
  function cast(event) {
    const rect = api.domElement.getBoundingClientRect();
    pointer.set(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
    raycaster.setFromCamera(pointer, api.camera);
  }
  // what is under the pointer: a pair of dots, a bond, or an atom
  function pick({ pairs = true, bonds = true, atoms = true } = {}) {
    if (pairs) {
      const sprites = api.pairSprites();
      const hit = raycaster.intersectObjects(sprites.map(p => p.sprite), false)[0];
      if (hit) {
        const p = sprites.find(q => q.sprite === hit.object);
        const dot = sprites.filter(q => q.index === p.index).indexOf(p);
        return { kind: 'pair', index: p.index, dot };
      }
    }
    const candidates = [];
    if (bonds) for (const m of api.bondMeshes()) candidates.push(m);
    if (atoms) for (const m of api.atomMeshes()) candidates.push(m);
    const hit = raycaster.intersectObjects(candidates, false)[0];
    if (!hit) return null;
    if (hit.object.userData.atom !== undefined) return { kind: 'atom', index: hit.object.userData.atom };
    return { kind: 'bond', edge: hit.object.userData.edge };
  }

  // ---------- points in space ----------
  const positions = () => api.getCurrent().geo.pos.map(p => new THREE.Vector3(...p));
  const radiusOf = i => { const m = api.atomMeshes()[i]; return m ? m.scale.x : 0.3; };
  const camBack = () => new THREE.Vector3().subVectors(api.camera.position, api.controls().target).normalize();
  const edgeIndex = (i, j) => api.getCurrent().r.edges.findIndex(e => key(e[0], e[1]) === key(i, j));
  const orderNow = (i, j) => { const k = edgeIndex(i, j); return k < 0 ? 0 : api.getStructure().orders[k]; };

  function midpoint(i, j) {
    const pos = positions();
    return pos[i].clone().add(pos[j]).multiplyScalar(0.5);
  }
  // where an arrow leaves from: the dot pair it was dragged from (or the atom's first pair), or the
  // middle of the bond
  function tailPoint(from, dot) {
    if ('pair' in from) {
      const mine = api.pairSprites().filter(p => p.index === from.pair);
      const p = mine[Math.min(dot || 0, mine.length - 1)];
      if (p) return p.sprite.position.clone();
      const pos = positions();
      return pos[from.pair].clone().addScaledVector(camBack(), radiusOf(from.pair) * 1.2);
    }
    return midpoint(from.bond[0], from.bond[1]);
  }
  // The ends and the bow of an arrow. An arrow bows across its chord in the plane of the screen,
  // away from the middle of the molecule, and at least half an angstrom so a short arrow still
  // reads as a curve. It points at an atom's rim, or at the middle of a bond that is there; a bond
  // yet to be made is pointed at by its far atom, as on paper. An arrow from a bond leaves from the
  // bond's side rather than through it.
  function shape(from, dot, to, centre, free) {
    const pos = positions();
    const back = camBack();
    const p0 = tailPoint(from, dot);
    let aim, atom = -1;
    if (free) aim = free.clone();
    else if ('atom' in to) { atom = to.atom; aim = pos[atom].clone(); }
    else {
      const [k, l] = to.bond;
      if (orderNow(k, l) > 1e-9) aim = midpoint(k, l);
      else {
        const own = 'pair' in from ? from.pair : (from.bond.includes(k) ? k : l);
        atom = own === k ? l : k;
        aim = pos[atom].clone();
      }
    }
    const chord = aim.clone().sub(p0);
    const L = chord.length() || 0.1;
    let side = new THREE.Vector3().crossVectors(back, chord).normalize();
    if (side.lengthSq() < 1e-6) side = new THREE.Vector3(0, 1, 0);
    const mid = p0.clone().add(aim).multiplyScalar(0.5);
    if (side.dot(mid.clone().sub(centre)) < 0) side.negate();
    if ('bond' in from) p0.addScaledVector(side, 0.1).addScaledVector(back, 0.1);
    const ctrl = mid.clone().addScaledVector(side, Math.max(0.32 * L, 0.5)).addScaledVector(back, 0.12 * L + 0.1);
    let p1;
    if (free) p1 = aim;
    else if (atom >= 0) p1 = aim.clone().add(ctrl.clone().sub(aim).normalize().multiplyScalar(radiusOf(atom) * 1.12));
    else p1 = aim.clone().addScaledVector(side, 0.1).addScaledVector(back, 0.1);
    return { p0, ctrl, p1 };
  }

  // ---------- tails and heads into arrows ----------
  function tailOf(hit) {
    if (hit.kind === 'pair') return { pair: hit.index };
    const [i, j] = api.getCurrent().r.edges[hit.edge];
    return { bond: [i, j] };
  }
  // the head an arrow from `from` means when released on `hit`
  function headOf(from, hit) {
    const edges = api.getCurrent().r.edges;
    const pos = positions();
    if ('pair' in from) {
      const i = from.pair;
      if (hit.kind === 'atom') return hit.index === i ? null : { bond: [i, hit.index] };
      if (hit.kind === 'bond') { const [k, l] = edges[hit.edge]; return k === i || l === i ? { bond: [k, l] } : null; }
      return null;
    }
    const [i, j] = from.bond;
    if (hit.kind === 'atom') {
      const k = hit.index;
      if (k === i || k === j) return { atom: k };
      // a third atom: the new bond goes from whichever end of the old bond is nearer to it
      const di = pos[i].distanceTo(pos[k]), dj = pos[j].distanceTo(pos[k]);
      return { bond: [di <= dj ? i : j, k] };
    }
    if (hit.kind === 'bond') { const [k, l] = edges[hit.edge]; return { bond: [k, l] }; }
    return null;
  }

  // ---------- drawing ----------
  function arrowMesh({ p0, ctrl, p1 }, mat) {
    const c = new THREE.QuadraticBezierCurve3(p0, ctrl, p1);
    const g = new THREE.Group();
    const total = c.getLength();
    const tEnd = Math.max(0.05, 1 - HEAD_LEN / Math.max(total, 1e-6));
    const pts = [];
    for (let q = 0; q <= 24; q++) pts.push(c.getPoint((q / 24) * tEnd));
    g.add(new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), 24, TUBE, 10, false), mat));
    const tip = c.getPoint(1), base = c.getPoint(tEnd);
    const cone = new THREE.Mesh(new THREE.ConeGeometry(HEAD_R, HEAD_LEN, 16), mat);
    cone.position.copy(base).add(tip).multiplyScalar(0.5);
    cone.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), tip.clone().sub(base).normalize());
    g.add(cone);
    return g;
  }
  function centreOf() {
    const pos = positions();
    const c = new THREE.Vector3();
    pos.forEach(p => c.add(p));
    return c.multiplyScalar(1 / Math.max(1, pos.length));
  }

  // rebuilt every frame: the dots move with the camera, so the arrows follow them
  function update() {
    group.children.forEach(ch => ch.traverse(o => { if (o.geometry) o.geometry.dispose(); }));
    group.clear();
    if (!api.getCurrent()) { drag = null; pending.length = 0; return; }
    const centre = centreOf();
    for (const p of pending) {
      const m = arrowMesh(shape(p.arrow.from, p.dot, p.arrow.to, centre, null), material);
      m.userData.pending = p;
      group.add(m);
    }
    if (drag) {
      const s = shape(drag.tail.from, drag.tail.dot, drag.head, centre, drag.head ? null : drag.point);
      group.add(arrowMesh(s, drag.head ? material : ghostMaterial));
    }
    const now = performance.now();
    for (let k = flashes.length - 1; k >= 0; k--) {
      const f = flashes[k];
      const t = (now - f.start) / 1600;
      const mesh = api.atomMeshes()[f.index];
      if (t >= 1 || !mesh) { api.scene.remove(f.mesh); f.mesh.material.dispose(); flashes.splice(k, 1); continue; }
      f.mesh.material.opacity = 0.35 * (1 - t);
      f.mesh.position.copy(mesh.position);
      f.mesh.scale.setScalar(mesh.scale.x * (1.5 + 0.3 * t));
    }
  }

  // a red halo, fading, on the atoms a refusal names
  function flash(atoms) {
    for (const i of atoms || []) {
      const mesh = api.atomMeshes()[i];
      if (!mesh) continue;
      const halo = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), new THREE.MeshBasicMaterial({ color: 0xc0392b, transparent: true, opacity: 0.35, depthWrite: false }));
      halo.position.copy(mesh.position);
      halo.scale.setScalar(mesh.scale.x * 1.5);
      api.scene.add(halo);
      flashes.push({ mesh: halo, index: i, start: performance.now() });
    }
  }

  // ---------- notes ----------
  function note(extra) {
    if (extra) { api.setNote(extra); return; }
    if (!pending.length) { api.setNote('to move electrons, drag from a pair of dots or from a bond to an atom or a bond'); return; }
    const node = api.getCurrent().node;
    const words = pending.map(p => describeArrow(node, p.arrow)).join('; ');
    api.setNote(`${pending.length} arrow${pending.length === 1 ? '' : 's'}: ${words}. Enter makes the move; Escape clears; click an arrow to remove it`);
  }

  // ---------- pointer ----------
  const dragPlane = new THREE.Plane();
  api.domElement.addEventListener('pointerdown', e => {
    if (!api.isAllowed() || e.button !== 0) return;
    cast(e);
    // a click on an arrow removes it
    const arrowHit = raycaster.intersectObjects(group.children.filter(c => c.userData.pending).flatMap(c => c.children), false)[0];
    if (arrowHit) {
      const p = arrowHit.object.parent.userData.pending;
      pending.splice(pending.indexOf(p), 1);
      note();
      api.onChange(pending.length);
      e.stopImmediatePropagation();
      return;
    }
    const hit = pick({ atoms: false });
    if (!hit) return;
    const from = tailOf(hit);
    const bad = tailShape(api.getCurrent().node, api.getStructure(), from);
    if (bad) { note(bad.reason); flash(bad.atoms); e.stopImmediatePropagation(); return; }
    const p0 = tailPoint(from, hit.dot);
    dragPlane.setFromNormalAndCoplanarPoint(camBack(), p0);
    drag = { tail: { from, dot: hit.dot || 0 }, point: p0.clone(), head: null };
    api.controls().enabled = false;
    try { api.domElement.setPointerCapture(e.pointerId); } catch (_) { /* fine without */ }
    e.stopImmediatePropagation();
  }, true);

  api.domElement.addEventListener('pointermove', e => {
    if (!drag) return;
    cast(e);
    const hit = pick({ pairs: false });
    drag.head = hit ? headOf(drag.tail.from, hit) : null;
    const p = new THREE.Vector3();
    if (raycaster.ray.intersectPlane(dragPlane, p)) drag.point.copy(p);
  });

  function endDrag(e) {
    if (!drag) return;
    const d = drag; drag = null;
    api.controls().enabled = true;
    try { api.domElement.releasePointerCapture(e.pointerId); } catch (_) { /* already released */ }
    if (!d.head) { note(); return; }
    const arrow = { from: d.tail.from, to: d.head };
    const bad = arrowShape(api.getCurrent().node, api.getStructure(), [...pending.map(p => p.arrow), arrow]);
    if (bad) { note(bad.reason); flash(bad.atoms); return; }
    pending.push({ arrow, dot: d.tail.dot });
    note();
    api.onChange(pending.length);
  }
  api.domElement.addEventListener('pointerup', endDrag);
  api.domElement.addEventListener('pointercancel', endDrag);

  window.addEventListener('keydown', e => {
    if (['INPUT', 'SELECT', 'TEXTAREA', 'BUTTON'].includes(e.target.tagName) || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === 'Enter' && pending.length) { commit(); e.preventDefault(); }
    else if (e.key === 'Escape' && (pending.length || drag)) { clear(); }
  });

  function commit() {
    if (!pending.length || !api.isAllowed()) return;
    api.onCommit(pending.map(p => p.arrow));
  }
  function clear() {
    pending.length = 0;
    drag = null;
    api.controls().enabled = true;
    note();
    api.onChange(0);
  }
  // arrows given ready-made (the textbook arrows of a reaction scene)
  function set(list) {
    pending.length = 0;
    const node = api.getCurrent().node;
    const bad = arrowShape(node, api.getStructure(), list);
    if (bad) { note(bad.reason); flash(bad.atoms); return; }
    list.forEach(arrow => pending.push({ arrow, dot: 0 }));
    note();
    api.onChange(pending.length);
  }

  return {
    update, commit, clear, flash, set,
    count: () => pending.length,
    arrows: () => pending.map(p => p.arrow),
    note: text => note(text),
    // the structure on screen changed under the arrows: they no longer point at anything real
    dropIfStale: () => { if (pending.length || drag) clear(); },
  };
}
