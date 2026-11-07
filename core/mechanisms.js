// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Mechanisms: one entry per mechanism, its steps in order (ledger items 72 and 85). Each step is a
// reaction scene of core/reactions.js, with its own reagents and textbook arrows, so every step is
// still checked on its own (test/reactions.js). What links two steps is the carry-over: the products
// of one step that are reagents of the next stay where they are; what the next step does not use
// leaves; what it adds comes in. Every step plays forward and backward.
//
// What leaves or comes in goes as a whole group: atoms joined by bonds or by the tie of a loose ion
// (a sodium chloride, a sodium amide), so nothing on screen is ever a bare ion (item 16).

import { REACTIONS } from './reactions.js';
import { element } from './elements.js';
import { buildNode, enumerate } from './lewis.js';
import { relax } from './geometry.js';
import { add, sub, scale, dot } from './vectors.js';

const scene = name => {
  const rx = REACTIONS.find(x => x.name === name);
  if (!rx) throw new Error(`no reaction scene called ${name}`);
  return rx;
};
const step = (name, label) => ({ ...scene(name), label });

export const MECHANISMS = [
  { name: 'SN2: hydroxide and chloromethane', steps: [step('SN2: hydroxide and chloromethane', 'hydroxide comes in from behind, chloride leaves')] },
  {
    name: 'SN1: tert-butyl bromide and water',
    steps: [
      step('SN1 step 1: tert-butyl bromide lets go of bromide', 'tert-butyl bromide lets go of bromide'),
      step('SN1 step 2: water adds to the cation', 'water adds to the cation'),
      step('SN1 step 3: water takes the proton', 'water takes the proton'),
    ],
  },
  { name: 'E2: hydroxide and ethyl bromide', steps: [step('E2: hydroxide and ethyl bromide', 'hydroxide takes a proton, the C=C forms, bromide leaves')] },
  { name: 'acid and base: acetic acid and hydroxide', steps: [step('acid and base: acetic acid and hydroxide', 'hydroxide takes the acid proton')] },
  {
    name: 'aldol addition: acetaldehyde and hydroxide',
    steps: [
      step('enolate: hydroxide takes the alpha proton of acetaldehyde', 'hydroxide takes the alpha proton'),
      step('aldol addition: the enolate adds to acetaldehyde', 'the enolate adds to acetaldehyde'),
    ],
  },
  {
    name: 'Michael addition: acetaldehyde enolate and methyl vinyl ketone',
    steps: [
      step('enolate: hydroxide takes the alpha proton of acetaldehyde', 'hydroxide takes the alpha proton'),
      step('Michael addition: the enolate adds to methyl vinyl ketone', 'the enolate adds to methyl vinyl ketone'),
    ],
  },
  {
    name: 'benzyne: chlorobenzene and sodium amide',
    steps: [
      step('benzyne step 1: amide takes the proton next to Cl, chloride leaves', 'amide takes the proton next to Cl, chloride leaves'),
      step('benzyne step 2: amide adds to benzyne', 'amide adds to benzyne'),
      step('benzyne step 3: the ring carbon takes a proton from ammonia (aniline)', 'the ring carbon takes a proton from ammonia'),
    ],
  },
  { name: 'Diels Alder: butadiene and ethylene', steps: [step('Diels Alder: butadiene and ethylene', 'three pairs move round the ring at once')] },
];

// ---------- the carry-over from one step to the next ----------

const loose = (spec, i) => spec.atoms[i].ion === -1 || !!element(spec.atoms[i].el).ionic;

// Groups (atoms joined by bonds or ties) and, inside each, molecules (atoms joined by bonds; a loose
// ion is a molecule of its own).
function parts(spec) {
  const n = spec.atoms.length;
  const unions = bondsOnly => {
    const parent = spec.atoms.map((_, i) => i);
    const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
    spec.adj.forEach(([i, j]) => { if (!bondsOnly || (!loose(spec, i) && !loose(spec, j))) parent[find(i)] = find(j); });
    const by = new Map();
    for (let i = 0; i < n; i++) { const r = find(i); if (!by.has(r)) by.set(r, []); by.get(r).push(i); }
    return [...by.values()];
  };
  const mols = unions(true);
  const molOf = new Array(n);
  mols.forEach((m, q) => m.forEach(i => { molOf[i] = q; }));
  const nbrs = spec.atoms.map(() => []);
  spec.adj.forEach(([i, j]) => { if (!loose(spec, i) && !loose(spec, j)) { nbrs[i].push(j); nbrs[j].push(i); } });
  const charge = m => m.reduce((s, i) => s + (spec.atoms[i].written || 0), 0);
  // a key no renaming of atoms can change: each atom's element, its bond count and its neighbours'
  // elements, sorted, and the molecule's charge
  const molKey = m => m.map(i => `${spec.atoms[i].el}${nbrs[i].length}:${nbrs[i].map(j => spec.atoms[j].el).sort().join('')}`).sort().join(' ') + `|${charge(m)}`;
  const groups = unions(false).map(members => {
    const ms = [...new Set(members.map(i => molOf[i]))];
    return { members, mols: ms, key: ms.map(q => molKey(mols[q])).sort().join(' + '), charge: charge(members) };
  });
  return { mols, molKey: mols.map(molKey), groups, nbrs };
}

// A renaming of the atoms of molecule `b` (in spec B) onto molecule `a` (in spec A) that keeps every
// element and every bond, or null.
function iso(A, pa, a, B, pb, b) {
  if (a.length !== b.length) return null;
  // visit b's atoms so that each one after the first has a neighbour already placed
  const order = [], seen = new Set();
  for (const s of b) {
    if (seen.has(s)) continue;
    const queue = [s]; seen.add(s);
    while (queue.length) { const i = queue.shift(); order.push(i); for (const j of pb.nbrs[i]) if (!seen.has(j)) { seen.add(j); queue.push(j); } }
  }
  const map = new Map(), used = new Set();
  const fits = (i, c) => A.atoms[c].el === B.atoms[i].el && pa.nbrs[c].length === pb.nbrs[i].length &&
    pb.nbrs[i].every(j => !map.has(j) || pa.nbrs[c].includes(map.get(j)));
  const go = q => {
    if (q === order.length) return true;
    const i = order[q];
    for (const c of a) {
      if (used.has(c) || !fits(i, c)) continue;
      map.set(i, c); used.add(c);
      if (go(q + 1)) return true;
      map.delete(i); used.delete(c);
    }
    return false;
  };
  return go(0) ? map : null;
}

// From the spec on screen (`from`, the products of a step) to the spec of the next step's reagents
// (`to`): which atom of `from` each atom of `to` is (-1 for an atom that comes in), and the groups
// that leave (atoms of `from`) and come in (atoms of `to`). Groups carry over whole: a group of `to`
// is matched to a group of `from` with the same molecules, and each molecule renamed onto its match.
export function carryOver(from, to) {
  const pf = parts(from), pt = parts(to);
  const map = new Array(to.atoms.length).fill(-1);
  const usedGroup = new Set();
  const arriving = [];
  // biggest groups first, so a substrate is matched before a stray small group could take its place
  const order = pt.groups.map((g, q) => q).sort((x, y) => pt.groups[y].members.length - pt.groups[x].members.length);
  for (const gq of order) {
    const g = pt.groups[gq];
    let done = false;
    for (let fq = 0; fq < pf.groups.length && !done; fq++) {
      const h = pf.groups[fq];
      if (usedGroup.has(fq) || h.key !== g.key) continue;
      // rename each molecule of g onto a molecule of h
      const local = new Map(), usedMol = new Set();
      const ok = g.mols.every(tm => {
        for (const fm of h.mols) {
          if (usedMol.has(fm) || pf.molKey[fm] !== pt.molKey[tm]) continue;
          const m = iso(from, pf, pf.mols[fm], to, pt, pt.mols[tm]);
          if (m) { usedMol.add(fm); m.forEach((c, i) => local.set(i, c)); return true; }
        }
        return false;
      });
      if (!ok) continue;
      local.forEach((c, i) => { map[i] = c; });
      usedGroup.add(fq);
      done = true;
    }
    if (!done) arriving.push(g.members);
  }
  const leaving = pf.groups.filter((_, fq) => !usedGroup.has(fq)).map(h => h.members);
  const chargeOf = (spec, list) => list.reduce((s, m) => s + m.reduce((t, i) => t + (spec.atoms[i].written || 0), 0), 0);
  return { map, leaving, arriving, leavingCharge: chargeOf(from, leaving), arrivingCharge: chargeOf(to, arriving) };
}

// The shape the next step starts from: the atoms that carry over where they are on screen, each group
// that comes in in its own shape set down beside them (to the right as seen, `right` being the
// screen's right in the scene's coordinates), then all of it relaxed together. Returns the carry-over
// and the positions of `to`'s atoms.
export function nextShape(from, fromPos, to, right = [1, 0, 0]) {
  const carry = carryOver(from, to);
  const start = carry.map.map(c => (c >= 0 ? [...fromPos[c]] : null));
  const kept = carry.map.map((c, i) => i).filter(i => carry.map[i] >= 0);
  const along = p => dot(p, right);
  let edge = kept.length ? Math.max(...kept.map(i => along(start[i]))) : 0;
  const mid = kept.length ? scale(kept.reduce((s, i) => add(s, start[i]), [0, 0, 0]), 1 / kept.length) : [0, 0, 0];
  for (const members of carry.arriving) {
    const at = new Map(members.map((i, q) => [i, q]));
    const sub1 = {
      atoms: members.map(i => ({ ...to.atoms[i] })),
      adj: to.adj.filter(([i, j]) => at.has(i) && at.has(j)).map(([i, j]) => [at.get(i), at.get(j)]),
      charge: members.reduce((s, i) => s + (to.atoms[i].written || 0), 0),
    };
    const node = buildNode(sub1);
    const own = relax(node, enumerate(node)).pos;
    const ownMid = scale(own.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / own.length);
    const gap = edge - Math.min(...own.map(along)) + 2.8;
    let shift = sub(mid, ownMid);
    shift = add(shift, scale(right, gap - dot(shift, right)));
    members.forEach((i, q) => { start[i] = add(own[q], shift); });
    edge = Math.max(edge, ...members.map(i => along(start[i])));
  }
  const node = buildNode(to);
  const pos = relax(node, enumerate(node), { start }).pos;
  return { carry, start, pos };
}
