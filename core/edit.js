// Editing a molecule: pure operations on a spec { atoms: [{el, ion?}], adj: [[i, j]], charge }.
// Every operation returns a new spec (and, where atoms are added, a suggested position for each new
// atom so the geometry can start near its parent). The engine then recomputes everything else:
// bond orders, lone pairs and charges are never set by hand (ledger items 7 and 12). Unsaturation
// comes from the electron count: build a skeleton, fill hydrogens, take two away, and the engine
// finds the double bond.

import { element } from './elements.js';

// how many bonds an atom usually carries, used only by fillHydrogens
const USUAL_BONDS = { H: 1, B: 3, C: 4, N: 3, O: 2, F: 1, Si: 4, P: 3, S: 2, Cl: 1, Br: 1, I: 1, Se: 2, Al: 3 };

export function emptySpec() { return { atoms: [], adj: [], charge: 0 }; }

// an edited molecule is one pool of electrons with one overall charge: the charges a SMILES wrote on
// its atoms are dropped, so plus and minus on the page mean what they say
const plain = atoms => atoms.map(a => (a.ion ? { el: a.el, ion: a.ion } : { el: a.el }));

// a spectator anion (a loose Cl-) that the student bonds to something stops being a spectator: it
// becomes an ordinary atom and the engine works out its bonds and charge with the rest
const joined = (atoms, i) => (atoms[i].ion ? { el: atoms[i].el } : atoms[i]);

export function neighbours(spec, i) {
  return spec.adj.filter(([a, b]) => a === i || b === i).map(([a, b]) => (a === i ? b : a));
}

// add an atom of element `el` bonded to atom `parent` (or free, when parent is null)
export function addAtom(spec, el, parent = null) {
  element(el); // throws on an unknown element
  const atoms = [...plain(spec.atoms), { el }];
  const adj = [...spec.adj];
  const index = atoms.length - 1;
  if (parent !== null) { adj.push([parent, index]); atoms[parent] = joined(atoms, parent); }
  return { spec: { atoms, adj, charge: spec.charge }, index };
}

// remove an atom; by default the hydrogens that hang only on it go with it
export function removeAtom(spec, i, withHydrogens = true) {
  const gone = new Set([i]);
  if (withHydrogens && spec.atoms[i].el !== 'H') {
    for (const j of neighbours(spec, i)) if (spec.atoms[j].el === 'H' && neighbours(spec, j).length === 1) gone.add(j);
  }
  const keep = spec.atoms.map((_, k) => !gone.has(k));
  const newIndex = []; let q = 0;
  keep.forEach(k => newIndex.push(k ? q++ : -1));
  const atoms = plain(spec.atoms.filter((_, k) => keep[k]));
  const adj = spec.adj.filter(([a, b]) => keep[a] && keep[b]).map(([a, b]) => [newIndex[a], newIndex[b]]);
  return { atoms, adj, charge: spec.charge, removed: [...gone].sort((a, b) => a - b) };
}

// make or unmake a bond between two existing atoms (ring closure, or cutting a ring); a pair that is
// only the tie holding a spectator anion beside its partner becomes a real bond
export function toggleBond(spec, i, j) {
  if (i === j) return spec;
  const key = ([a, b]) => `${Math.min(a, b)}-${Math.max(a, b)}`;
  const want = key([i, j]);
  const has = spec.adj.some(e => key(e) === want);
  const atoms = plain(spec.atoms);
  const tie = has && !!(atoms[i].ion || atoms[j].ion);
  let adj = spec.adj;
  if (!has) adj = [...spec.adj, [Math.min(i, j), Math.max(i, j)]];
  else if (!tie) adj = spec.adj.filter(e => key(e) !== want);
  if (!has || tie) { atoms[i] = joined(atoms, i); atoms[j] = joined(atoms, j); }
  return { atoms, adj, charge: spec.charge };
}

export function setCharge(spec, charge) { return { atoms: plain(spec.atoms), adj: spec.adj, charge }; }

// give every heavy atom its usual number of bonds by adding hydrogens; returns the new spec and the
// indices of the atoms added, with the parent of each
export function fillHydrogens(spec) {
  let out = { atoms: [...spec.atoms], adj: [...spec.adj], charge: spec.charge };
  const added = [];
  spec.atoms.forEach((a, i) => {
    if (a.ion || element(a.el).ionic) return;
    const usual = USUAL_BONDS[a.el];
    if (usual === undefined) return;
    const have = neighbours(spec, i).length;
    for (let q = have; q < usual; q++) {
      const r = addAtom(out, 'H', i);
      out = r.spec;
      added.push({ index: r.index, parent: i });
    }
  });
  return { spec: out, added };
}

// a starting position for a new atom next to its parent: away from the parent's other neighbours,
// at a plausible bond length (a free atom is placed by the page, which knows where the student is
// looking from)
export function placeNear(positions, parentIndex, elNew, spec) {
  const r = element(elNew).radius;
  const p = positions[parentIndex];
  const d = r + element(spec.atoms[parentIndex].el).radius;
  const nbrs = neighbours(spec, parentIndex).filter(j => j < positions.length);
  let away = [0, 0, 0];
  for (const j of nbrs) {
    const v = [positions[j][0] - p[0], positions[j][1] - p[1], positions[j][2] - p[2]];
    const l = Math.hypot(...v) || 1;
    away = [away[0] - v[0] / l, away[1] - v[1] / l, away[2] - v[2] / l];
  }
  let l = Math.hypot(...away);
  if (l < 1e-6) {
    // no neighbours, or neighbours that cancel out (a linear pair): go sideways
    const v = nbrs.length ? [positions[nbrs[0]][0] - p[0], positions[nbrs[0]][1] - p[1], positions[nbrs[0]][2] - p[2]] : [1, 0, 0];
    away = Math.abs(v[0]) < Math.abs(v[1]) ? [0, -v[2], v[1]] : [-v[1], v[0], 0];
    l = Math.hypot(...away) || 1;
    if (!nbrs.length) { away = [1, 0, 0]; l = 1; }
  }
  // a touch of tilt so three new atoms on one parent do not all start on the same line
  const tilt = [0.17 * Math.sin(positions.length), 0.17 * Math.cos(positions.length * 1.3), 0.11];
  const dir = [away[0] / l + tilt[0], away[1] / l + tilt[1], away[2] / l + tilt[2]];
  const dl = Math.hypot(...dir);
  return [p[0] + d * dir[0] / dl, p[1] + d * dir[1] / dl, p[2] + d * dir[2] / dl];
}
