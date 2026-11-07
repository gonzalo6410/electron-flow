// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Names. The search box takes a name (common or systematic) or a SMILES; the names come from PubChem
// through the index that tools/build-names.mjs writes (web/names.json.gz). This file holds what does
// not depend on the page: comparing names, searching the index, the counter-ions a bare ion gets, and
// a key for a molecule's skeleton so the page can name what is on screen after an edit or a move.

import { Molecule } from 'openchemlib';
import { pieces, buildNode } from './lewis.js';
import { element } from './elements.js';

// How names are compared: lowercase, letters and digits only ("tert-Butanol" and "tertbutanol" match).
export function normName(s) { return s.toLowerCase().replace(/[^a-z0-9]/g, ''); }

// Charge balance (ledger item 16): a SMILES whose pieces add up to a charge gets counter-ions, a
// sodium for each negative charge and a chloride for each positive one, written after it so each ties
// to its partner (core/smiles.js). Returns the SMILES to read and what was added, for the status line.
export function balanceSmiles(smiles) {
  const text = smiles.trim();
  let m;
  try { m = Molecule.fromSmiles(text); } catch (e) { return { smiles: text, added: null }; }
  let q = 0;
  for (let i = 0; i < m.getAllAtoms(); i++) q += m.getAtomCharge(i);
  if (!q) return { smiles: text, added: null };
  const ion = q < 0 ? 'Na+' : 'Cl-';
  return { smiles: text + `.[${ion}]`.repeat(Math.abs(q)), added: { ion, count: Math.abs(q) } };
}

// The skeleton of a molecule as the engine sees it: which atoms, which are joined (ties to loose ions
// left out, so a loose ion is a piece of its own), and the charge each piece carries. Bond orders,
// lone pairs and charges inside a piece are left to the engine, so every structure of a node, and
// every way of writing its SMILES, gives the same key. Canonical through openchemlib's ID code, with
// every hydrogen an atom of its own.
export function skeletonKey(node) {
  const m = new Molecule(node.atoms.length, node.edges.length);
  const ids = node.atoms.map(a => m.addAtom(Molecule.getAtomicNoFromLabel(a.el)));
  const deg = node.atoms.map(() => 0);
  for (const [i, j] of node.edges) {
    if (node.fixed[i] || node.fixed[j]) continue;
    m.addBond(ids[i], ids[j], Molecule.cBondTypeSingle);
    deg[i]++; deg[j]++;
  }
  ids.forEach((id, i) => m.setAtomAbnormalValence(id, deg[i]));
  const charges = pieces(node).map(p => {
    const valence = p.members.reduce((s, i) => s + (node.fixed[i] ? 0 : element(node.atoms[i].el).valence), 0);
    const fixedCharge = p.members.reduce((s, i) => s + (node.fixed[i] ? node.fixed[i].charge : 0), 0);
    return valence - p.pool + fixedCharge;
  }).sort((a, b) => a - b);
  return m.getIDCode() + '|' + charges.join(',');
}

// 64-bit FNV-1a of the key, in base 36: what the index stores in place of the key itself.
export function keyHash(key) {
  let h = 0xcbf29ce484222325n;
  const p = 0x100000001b3n, mask = (1n << 64n) - 1n;
  for (let i = 0; i < key.length; i++) { h ^= BigInt(key.charCodeAt(i)); h = (h * p) & mask; }
  return h.toString(36);
}

// The index, ready to search. data.entries: [SMILES or null, names (first is shown), skeleton hash or
// the reason it has no SMILES, 0 for a common chemical, PubChem CID].
export function makeIndex(data) {
  const entries = data.entries.map(([smiles, names, key, tier, cid]) => ({ smiles, names, key, tier, cid, name: names[0] }));
  const flat = []; // { n: normalised name, e: entry, name }
  entries.forEach(e => e.names.forEach(name => { const n = normName(name); if (n) flat.push({ n, e, name }); }));
  flat.sort((a, b) => (a.n < b.n ? -1 : a.n > b.n ? 1 : 0));
  const byKey = new Map();
  for (const e of entries) if (e.smiles && !byKey.has(e.key)) byKey.set(e.key, e);
  return { entries, flat, byKey };
}

// Up to `limit` entries for what was typed: an entry whose own name (the one shown) equals it, then
// one with another name equal to it ("acetate" is the ion before it is a synonym of acetic acid),
// then own names starting with it, other names starting with it, and names with a word starting with
// it (from three letters on); within each, common chemicals first, then the shorter name, then the lower PubChem
// number. Each hit says which name matched.
export function search(index, text, limit = 8) {
  const q = normName(text);
  if (!q) return [];
  const { flat } = index;
  let lo = 0, hi = flat.length;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (flat[mid].n < q) lo = mid + 1; else hi = mid; }
  const best = new Map(); // entry -> { tier, name, n }
  const offer = (x, tier) => {
    const had = best.get(x.e);
    if (!had || tier < had.tier || (tier === had.tier && x.n.length < had.n.length)) best.set(x.e, { tier, name: x.name, n: x.n });
  };
  const own = x => x.name === x.e.name;
  for (let k = lo; k < flat.length && flat[k].n.startsWith(q); k++) {
    const x = flat[k];
    offer(x, x.n === q ? (own(x) ? 0 : 1) : (own(x) ? 2 : 3));
  }
  // inside a name, only from the start of one of its words ("butanol" finds tert-butanol; "thf" does
  // not find moth flakes)
  const atWord = x => {
    const words = x.name.toLowerCase().split(/[^a-z0-9]+/).filter(Boolean);
    for (let w = 1; w < words.length; w++) if (words.slice(w).join('').startsWith(q)) return true;
    return false;
  };
  if (q.length >= 3 && best.size < limit * 4) {
    for (const x of flat) if (!x.n.startsWith(q) && x.n.includes(q) && atWord(x)) offer(x, 4);
  }
  return [...best.entries()]
    .sort(([ea, a], [eb, b]) => a.tier - b.tier || ea.tier - eb.tier || a.n.length - b.n.length || ea.cid - eb.cid)
    .slice(0, limit)
    .map(([e, h]) => ({ entry: e, matched: h.name, exact: h.tier <= 1 }));
}

// The entry whose skeleton is the node's, or null.
export function nameFor(index, node) {
  if (!index) return null;
  try { return index.byKey.get(keyHash(skeletonKey(node))) || null; } catch (e) { return null; }
}

// Names for what is on screen, given its spec: the name of the whole if the list has it; else, when
// every atom carries its own charge (a scene, the products of a move), one name per molecule, a
// molecule being the atoms joined by bonds or by the tie of a loose ion (so the sodium and chloride
// left by an SN2 read as sodium chloride), with the formula for a molecule the list does not have.
// Null when nothing could be named.
export function namesFor(index, spec) {
  if (!index || !spec || !spec.atoms.length) return null;
  try {
    const whole = nameFor(index, buildNode(spec));
    if (whole) return [whole.name];
    if (!spec.atoms.every(a => Number.isInteger(a.written))) return null;
    const n = spec.atoms.length;
    const parent = spec.atoms.map((_, i) => i);
    const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
    spec.adj.forEach(([i, j]) => { parent[find(i)] = find(j); });
    const groups = new Map();
    for (let i = 0; i < n; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
    if (groups.size < 2) return null;
    let named = 0;
    const out = [...groups.values()].map(members => {
      const at = new Map(members.map((i, q) => [i, q]));
      const sub = {
        atoms: members.map(i => ({ ...spec.atoms[i] })),
        adj: spec.adj.filter(([i]) => at.has(i)).map(([i, j]) => [at.get(i), at.get(j)]),
        charge: members.reduce((t, i) => t + spec.atoms[i].written, 0),
      };
      const hit = nameFor(index, buildNode(sub));
      if (hit) { named++; return hit.name; }
      return formula(sub);
    });
    return named ? out : null;
  } catch (e) { return null; }
}

// a formula with its charge, carbon and hydrogen first: C4H9+
function formula(spec) {
  const counts = {};
  spec.atoms.forEach(a => { counts[a.el] = (counts[a.el] || 0) + 1; });
  const els = Object.keys(counts).sort((a, b) => {
    const rank = e => (e === 'C' ? 0 : e === 'H' ? 1 : 2);
    return rank(a) - rank(b) || (a < b ? -1 : 1);
  });
  const q = spec.charge;
  const mark = q === 0 ? '' : Math.abs(q) === 1 ? (q > 0 ? '+' : '-') : `${Math.abs(q)}${q > 0 ? '+' : '-'}`;
  return els.map(e => e + (counts[e] > 1 ? counts[e] : '')).join('') + mark;
}

// Why an entry has no SMILES, in plain words.
export function reasonText(key) {
  if (key.startsWith('el:')) return `it contains ${key.slice(3)}, which is not in the engine's table yet`;
  if (key === 'slow') return 'it has too many structures for the engine to finish quickly yet';
  if (key === 'ring') return 'no Lewis structure fits its ring';
  return 'the engine could not draw it';
}
