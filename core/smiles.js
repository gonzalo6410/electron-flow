// Molecules from SMILES, using openchemlib (BSD-3) only to read the text: which atoms, which atoms are
// joined, the overall charge. Bond orders, lone pairs and formal charges are recomputed by our own
// engine (core/lewis.js), so what the SMILES says about those is ignored.
//
// Counter-ions written as separate pieces ([Na+].[O-]...) are tied loosely to a partner (item 15):
// a cation to the atom carrying a negative charge in the SMILES (the first one, if several),
// failing that to the most electronegative atom; an anion the other way round.

import { Molecule } from 'openchemlib';
import { ELEMENTS, element } from './elements.js';

export function nodeFromSmiles(smiles) {
  const text = smiles.trim();
  if (!text) throw new Error('empty SMILES');
  let m;
  try { m = Molecule.fromSmiles(text); } catch (e) { throw new Error('could not read that SMILES'); }
  m.addImplicitHydrogens();
  const n = m.getAllAtoms();
  if (!n) throw new Error('no atoms');
  const atoms = [];
  for (let i = 0; i < n; i++) {
    const el = m.getAtomLabel(i);
    if (!ELEMENTS[el]) throw new Error(`element ${el} is not in the table yet`);
    atoms.push({ el, written: m.getAtomCharge(i) });
  }
  const adj = [];
  for (let b = 0; b < m.getAllBonds(); b++) adj.push([m.getBondAtom(0, b), m.getBondAtom(1, b)]);
  const charge = atoms.reduce((s, a) => s + a.written, 0);

  // tie loose ions to partners: each cation to its own negatively written atom (the first free one),
  // each loose anion to a positively written atom, failing that to the most or least electronegative
  const bonded = new Set(adj.flat());
  const taken = new Set();
  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    const loose = element(a.el).ionic || (!bonded.has(i) && n > 1 && a.written !== 0);
    if (!loose || bonded.has(i)) continue;
    const wantSign = a.written > 0 ? -1 : 1;
    const candidates = atoms.map((b, j) => j).filter(j => j !== i && Math.sign(atoms[j].written) === wantSign && !element(atoms[j].el).ionic && !atoms[j].loose);
    let partner = candidates.find(j => !taken.has(j));
    if (partner === undefined && candidates.length) partner = candidates[0];
    if (partner === undefined) {
      let best = -1, bestEn = wantSign < 0 ? -1 : 99;
      atoms.forEach((b, j) => {
        if (j === i || element(b.el).ionic || b.el === 'H') return;
        const en = element(b.el).en;
        if (wantSign < 0 ? en > bestEn : en < bestEn) { best = j; bestEn = en; }
      });
      partner = best;
    }
    if (partner >= 0) { adj.push([i, partner]); taken.add(partner); }
  }
  // a single-atom piece written with a negative charge is a spectator anion
  return { atoms: atoms.map((a, i) => ({ el: a.el, ion: (!bonded.has(i) && a.written < 0 && !element(a.el).ionic) ? -1 : 0 })), adj, charge };
}

// A SMILES for a molecule as the engine sees it in one structure: the atoms, the bonds with that
// structure's whole orders (a pair at order 0 is left out), and that structure's formal charges.
// Hydrogens are written out only where they carry a charge or are alone.
export function toSmiles(node, structure) {
  const m = new Molecule(node.atoms.length + 8, node.edges.length + 8);
  const ids = node.atoms.map(a => m.addAtom(Molecule.getAtomicNoFromLabel(a.el)));
  node.edges.forEach(([i, j], k) => {
    const o = structure.orders[k];
    if (o <= 0) return;
    const b = m.addBond(ids[i], ids[j]);
    m.setBondType(b, o >= 3 ? Molecule.cBondTypeTriple : o >= 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeSingle);
  });
  structure.charge.forEach((q, i) => { if (q) m.setAtomCharge(ids[i], Math.round(q)); });
  return m.toSmiles();
}
