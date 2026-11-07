// Editing on the page: click an atom to select it, press an element key to add a bonded atom,
// Delete to remove it, shift-click another atom to bond or unbond, plus and minus for the charge.
// The engine recomputes the structure after every change; this file never sets a bond order.
// It talks to the viewer through a small api object (see attachEditor) so the viewer stays simple.
import * as THREE from 'three';
import { addAtom, removeAtom, toggleBond, setCharge, fillHydrogens, placeNear, neighbours } from './core/edit.js';

// element keys: one letter where one letter is enough, otherwise a letter the element is known by
const KEYS = { c: 'C', n: 'N', o: 'O', f: 'F', s: 'S', p: 'P', b: 'B', h: 'H', l: 'Cl', r: 'Br', i: 'I', k: 'K', a: 'Na', m: 'Mg', x: 'Xe', e: 'Se', g: 'Al', t: 'Si' };

export function attachEditor(api) {
  // api: { scene, camera, domElement, getSpec(), getPositions(), atomMeshes(), loadSpec(spec, start, name), setNote(text), isEditable() }
  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  let selected = null;
  const halo = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 24), new THREE.MeshBasicMaterial({ color: 0x111111, transparent: true, opacity: 0.18, depthWrite: false }));
  halo.visible = false;
  api.scene.add(halo);

  function atomAt(event) {
    const rect = api.domElement.getBoundingClientRect();
    pointer.set(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
    raycaster.setFromCamera(pointer, api.camera);
    const hit = raycaster.intersectObjects(api.atomMeshes(), false)[0];
    return hit ? hit.object.userData.atom : null;
  }

  function select(i) {
    selected = i;
    refreshHalo();
    note();
  }

  function refreshHalo() {
    const pos = api.getPositions();
    if (selected === null || selected >= pos.length) { halo.visible = false; selected = null; return; }
    const mesh = api.atomMeshes()[selected];
    halo.visible = true;
    halo.position.copy(mesh.position);
    halo.scale.setScalar(mesh.scale.x * 1.45);
  }

  function note() {
    if (!api.isEditable()) { api.setNote('transitions cannot be edited yet'); return; }
    if (selected === null) api.setNote('click an atom to edit, or press a letter to start a new atom');
    else {
      const spec = api.getSpec();
      api.setNote(`${spec.atoms[selected].el}${selected} selected: a letter adds an atom to it, Delete removes it, shift-click another atom bonds or unbonds`);
    }
  }

  // a press or drag of less than a few pixels counts as a click; drags orbit the camera
  let down = null;
  api.domElement.addEventListener('pointerdown', e => { down = [e.clientX, e.clientY, e.shiftKey]; });
  api.domElement.addEventListener('pointerup', e => {
    if (!down) return;
    const moved = Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4;
    const shift = down[2];
    down = null;
    if (moved || !api.isEditable()) return;
    const i = atomAt(e);
    if (shift && selected !== null && i !== null && i !== selected) {
      apply(toggleBond(api.getSpec(), selected, i), api.getPositions());
      return;
    }
    select(i);
  });

  // show an edited molecule, relaxed from the given positions; a change for which the engine finds
  // no Lewis structure is refused, and the molecule stays as it was
  function apply(spec, start) {
    try {
      api.loadSpec(spec, start);
      refreshHalo();
      note();
    } catch (err) {
      api.setNote(`that change is not allowed: ${err.message}`);
    }
  }

  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.tagName === 'TEXTAREA') return;
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    if (!api.isEditable()) return;
    const spec = api.getSpec();
    const pos = api.getPositions();
    const key = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (KEYS[key] && !e.shiftKey) {
      addElement(KEYS[key], selected);
      e.preventDefault();
    } else if ((e.key === 'Delete' || e.key === 'Backspace') && selected !== null) {
      const victim = selected;
      const next = removeAtom(spec, victim);
      const gone = new Set(next.removed);
      const heavyNb = neighbours(spec, victim).find(j => !gone.has(j));
      selected = heavyNb === undefined ? null : heavyNb - next.removed.filter(g => g < heavyNb).length;
      apply(next, pos.filter((_, k) => !gone.has(k)));
      e.preventDefault();
    } else if (e.key === '+' || e.key === '=') {
      apply(setCharge(spec, spec.charge + 1), pos);
    } else if (e.key === '-' || e.key === '_') {
      apply(setCharge(spec, spec.charge - 1), pos);
    } else if (e.key === 'Escape') {
      select(null);
    } else if (e.key === 'H' && e.shiftKey) {
      fill();
      e.preventDefault();
    }
  });

  // where a free atom goes: where it was dropped, at the depth of the middle of the molecule; from the
  // keyboard or the table, a little to the right of everything as seen on the screen
  function freePosition(event) {
    const pos = api.getPositions().map(p => new THREE.Vector3(...p));
    const centre = pos.reduce((s, p) => s.add(p), new THREE.Vector3()).divideScalar(Math.max(1, pos.length));
    if (event) {
      const rect = api.domElement.getBoundingClientRect();
      pointer.set(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
      raycaster.setFromCamera(pointer, api.camera);
      const toViewer = api.camera.getWorldDirection(new THREE.Vector3()).negate();
      const hit = raycaster.ray.intersectPlane(new THREE.Plane().setFromNormalAndCoplanarPoint(toViewer, centre), new THREE.Vector3());
      if (hit) return hit.toArray();
    }
    const right = new THREE.Vector3().setFromMatrixColumn(api.camera.matrixWorld, 0);
    const reach = pos.reduce((m, p) => Math.max(m, p.clone().sub(centre).dot(right)), 0);
    return centre.addScaledVector(right, reach + 3).toArray();
  }

  function addElement(el, parent, event = null) {
    const spec = api.getSpec(), pos = api.getPositions();
    const r = addAtom(spec, el, parent);
    const start = [...pos, parent === null ? freePosition(event) : placeNear(pos, parent, el, r.spec)];
    selected = r.index;
    apply(r.spec, start);
  }

  // drops from the periodic table: onto an atom bonds to it, onto empty space makes a free atom there
  api.domElement.addEventListener('dragover', e => { e.preventDefault(); });
  api.domElement.addEventListener('drop', e => {
    e.preventDefault();
    const el = e.dataTransfer.getData('text/element');
    if (!el || !api.isEditable()) return;
    addElement(el, atomAt(e), e);
  });

  // give every atom its usual number of bonds with hydrogens, each placed beside its parent
  function fill() {
    const spec = api.getSpec(), pos = api.getPositions();
    const r = fillHydrogens(spec);
    const start = [...pos];
    for (const a of r.added) start.push(placeNear(start, a.parent, 'H', r.spec));
    apply(r.spec, start);
  }

  return {
    addElement,
    fill,
    selectedAtom() { return selected; },
    onLoaded() { select(null); },
    refresh() { refreshHalo(); },
    clear() { selected = null; apply({ atoms: [], adj: [], charge: 0 }, []); },
  };
}
