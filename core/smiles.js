// Molecules from SMILES, using openchemlib (BSD-3) only to read the text: which atoms, which atoms are
// joined, the overall charge. Bond orders, lone pairs and formal charges are recomputed by our own
// engine (core/lewis.js), so what the SMILES says about those is ignored. Isotope marks and stereo
// marks are read and dropped. Written charges are kept per atom only so that each separate piece
// keeps the electrons it was written with (core/lewis.js, pieces).
//
// Counter-ions written as separate pieces ([Na+].[O-]...) are tied loosely to a partner (item 15):
// a cation to an atom carrying a negative charge in the SMILES (the first one not already taken by
// another ion, so the two Na+ of carbonate sit by two different O), failing that to the most
// electronegative atom; an anion the other way round.

import { Molecule } from 'openchemlib';
import { ELEMENTS, element } from './elements.js';

export function nodeFromSmiles(smiles) {
  const text = smiles.trim();
  if (!text) throw new Error('empty SMILES');
  // openchemlib is forgiving about an unmatched bracket or a bond with no atom after it; a student
  // who typed one should hear about it rather than see a different molecule
  if (count(text, '(') !== count(text, ')')) throw new Error('could not read that SMILES: unmatched parenthesis');
  if (count(text, '[') !== count(text, ']')) throw new Error('could not read that SMILES: unmatched bracket');
  if (/[=#\-]$/.test(text) || /^[=#]/.test(text)) throw new Error('could not read that SMILES: a bond with nothing on one side');
  let m;
  try { m = Molecule.fromSmiles(text); } catch (e) { throw new Error('could not read that SMILES'); }
  if (m.isFragment()) throw new Error('could not read that SMILES: a wildcard atom (*) is not a molecule');
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

  // a single-atom piece written with the charge that fills its shell (Cl-, O2-, H-) is a spectator
  // anion; a lone [C-] or [O-] is not, it is an atom short of its octet and the engine treats it so
  const bonded = new Set(adj.flat());
  const spectator = atoms.map((a, i) => {
    const e = element(a.el);
    return !bonded.has(i) && !e.ionic && a.written < 0 && a.written === e.valence - e.ceiling;
  });

  // tie loose ions to partners
  const taken = new Set();
  for (let i = 0; i < n; i++) {
    const a = atoms[i];
    const loose = (element(a.el).ionic && !bonded.has(i)) || spectator[i];
    if (!loose || n < 2) continue;
    const wantSign = element(a.el).ionic ? -1 : 1;
    const charged = (b, j) => j !== i && Math.sign(b.written) === wantSign && !element(b.el).ionic;
    let partner = atoms.findIndex((b, j) => charged(b, j) && !taken.has(j));
    if (partner < 0) partner = atoms.findIndex(charged);
    if (partner < 0) {
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
  return {
    atoms: atoms.map((a, i) => ({ el: a.el, written: a.written, ion: spectator[i] ? -1 : 0 })),
    adj,
    charge,
  };
}

function count(text, ch) {
  let c = 0;
  for (const x of text) if (x === ch) c++;
  return c;
}

// A SMILES for a molecule as the engine sees it in one structure, written so that reading it back
// gives the same molecule: the atoms, every adjacent pair as a bond with that structure's order (a
// pair the octet ceiling has cut, like the two spare F of SF6, is written as a single bond, since the
// reader only uses a SMILES for which atoms are joined and how much charge each piece carries), that
// structure's formal charges, and the ties to loose ions as separate pieces. Hydrogens are written
// out only where they carry a charge or are alone. openchemlib would add hydrogens of its own to any
// atom with fewer bonds than it expects (the methyl radical would come out as methane, the nitrogen
// dioxide radical as nitrous acid, a bare carbon as methane), so every such atom is told that its
// bonds are all it has, and is then written in brackets: [CH3], [O]N=O, [C].
export function toSmiles(node, structure) {
  const m = new Molecule(node.atoms.length + 8, node.edges.length + 8);
  const ids = node.atoms.map(a => m.addAtom(Molecule.getAtomicNoFromLabel(a.el)));
  const bonds = node.atoms.map(() => 0);
  node.edges.forEach(([i, j], k) => {
    if (node.fixed[i] || node.fixed[j]) return; // a tie to a loose ion, not a bond
    const o = Math.max(1, Math.round(structure.orders[k]));
    const b = m.addBond(ids[i], ids[j]);
    m.setBondType(b, o >= 3 ? Molecule.cBondTypeTriple : o >= 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeSingle);
    bonds[i] += o; bonds[j] += o;
  });
  structure.charge.forEach((q, i) => { if (q) m.setAtomCharge(ids[i], Math.round(q)); });
  ids.forEach((id, i) => {
    if (m.getImplicitHydrogens(id) <= 0) return;
    m.setAtomAbnormalValence(id, bonds[i]);
    const still = m.getImplicitHydrogens(id); // a charged atom is allowed one bond more or less
    if (still > 0) m.setAtomAbnormalValence(id, Math.max(0, bonds[i] - still));
  });
  return m.toSmiles();
}
