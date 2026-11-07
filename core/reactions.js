// Reaction scenes: the reagents of a mechanism step placed together, with the arrows that make the
// step, for the picker and for the tests. The arrows name atoms by their place in the SMILES:
// 'a2' is the third heavy atom written, 'h2' a hydrogen on it. Ledger item 5.
//
// The arrows are the textbook ones; the engine checks them like any others (core/moves.js), so a
// scene is also a test that the engine agrees with the textbook on these steps.

export const REACTIONS = [
  {
    name: 'SN2: hydroxide and chloromethane',
    smiles: '[OH-].[Na+].CCl',
    arrows: [{ from: { pair: 'a0' }, to: { bond: ['a0', 'a2'] } }, { from: { bond: ['a2', 'a3'] }, to: { atom: 'a3' } }],
    products: ['CH4O', 'Cl-', 'Na+'],
  },
  {
    name: 'SN1 step 1: tert-butyl bromide lets go of bromide',
    smiles: 'CC(C)(C)Br',
    arrows: [{ from: { bond: ['a1', 'a4'] }, to: { atom: 'a4' } }],
    products: ['C4H9+', 'Br-'],
  },
  {
    name: 'SN1 step 2: water adds to the cation',
    smiles: 'C[C+](C)C.[Br-].O',
    arrows: [{ from: { pair: 'a5' }, to: { bond: ['a5', 'a1'] } }],
    products: ['C4H11O+', 'Br-'],
  },
  {
    name: 'SN1 step 3: water takes the proton',
    smiles: 'CC(C)(C)[OH2+].[Br-].O',
    arrows: [{ from: { pair: 'a6' }, to: { bond: ['a6', 'h4'] } }, { from: { bond: ['a4', 'h4'] }, to: { atom: 'a4' } }],
    products: ['C4H10O', 'H3O+', 'Br-'],
  },
  {
    name: 'E2: hydroxide and ethyl bromide',
    smiles: '[OH-].[Na+].CCBr',
    arrows: [
      { from: { pair: 'a0' }, to: { bond: ['a0', 'h2'] } },
      { from: { bond: ['a2', 'h2'] }, to: { bond: ['a2', 'a3'] } },
      { from: { bond: ['a3', 'a4'] }, to: { atom: 'a4' } },
    ],
    products: ['H2O', 'C2H4', 'Br-', 'Na+'],
  },
  {
    name: 'acid and base: acetic acid and hydroxide',
    smiles: 'CC(=O)O.[OH-].[Na+]',
    arrows: [{ from: { pair: 'a4' }, to: { bond: ['a4', 'h3'] } }, { from: { bond: ['a3', 'h3'] }, to: { atom: 'a3' } }],
    products: ['C2H3O2-', 'H2O', 'Na+'],
  },
  {
    name: 'enolate: hydroxide takes the alpha proton of acetaldehyde',
    smiles: 'CC=O.[OH-].[Na+]',
    arrows: [
      { from: { pair: 'a3' }, to: { bond: ['a3', 'h0'] } },
      { from: { bond: ['a0', 'h0'] }, to: { bond: ['a0', 'a1'] } },
      { from: { bond: ['a1', 'a2'] }, to: { atom: 'a2' } },
    ],
    products: ['C2H3O-', 'H2O', 'Na+'],
  },
  {
    name: 'aldol addition: the enolate adds to acetaldehyde',
    smiles: 'C=C[O-].[Na+].CC=O',
    arrows: [{ from: { bond: ['a0', 'a1'] }, to: { bond: ['a0', 'a5'] } }, { from: { bond: ['a5', 'a6'] }, to: { atom: 'a6' } }],
    products: ['C4H7O2-', 'Na+'],
  },
  {
    name: 'Michael addition: the enolate adds to methyl vinyl ketone',
    smiles: 'C=C[O-].[Na+].C=CC(C)=O',
    arrows: [
      { from: { bond: ['a0', 'a1'] }, to: { bond: ['a0', 'a4'] } },
      { from: { bond: ['a4', 'a5'] }, to: { bond: ['a5', 'a6'] } },
      { from: { bond: ['a6', 'a8'] }, to: { atom: 'a8' } },
    ],
    products: ['C6H9O2-', 'Na+'],
  },
  {
    name: 'Diels Alder: butadiene and ethylene',
    smiles: 'C=CC=C.C=C',
    arrows: [
      { from: { bond: ['a0', 'a1'] }, to: { bond: ['a0', 'a4'] } },
      { from: { bond: ['a4', 'a5'] }, to: { bond: ['a5', 'a3'] } },
      { from: { bond: ['a2', 'a3'] }, to: { bond: ['a1', 'a2'] } },
    ],
    products: ['C6H10'],
  },
];

// Turn the names in a scene's arrows into atom indices of its spec (a spec read from the scene's
// SMILES: heavy atoms first in the order written, hydrogens after, each next to its parent). The
// same name always means the same atom; 'h2.2' is a second hydrogen on heavy atom 2.
export function resolveArrows(spec, arrows) {
  const heavy = spec.atoms.map((a, i) => (a.el === 'H' ? -1 : i)).filter(i => i >= 0);
  const atom = name => {
    if (typeof name === 'number') return name;
    const [q, which = '1'] = name.slice(1).split('.').map(Number);
    const parent = heavy[q];
    if (parent === undefined) throw new Error(`no heavy atom ${q} in this scene`);
    if (name[0] === 'a') return parent;
    const hs = spec.adj.filter(([a, b]) => a === parent || b === parent).map(([a, b]) => (a === parent ? b : a)).filter(j => spec.atoms[j].el === 'H');
    const h = hs[which - 1];
    if (h === undefined) throw new Error(`no hydrogen ${which} on heavy atom ${q}`);
    return h;
  };
  return arrows.map(a => ({
    from: 'pair' in a.from ? { pair: atom(a.from.pair) } : { bond: a.from.bond.map(atom) },
    to: 'atom' in a.to ? { atom: atom(a.to.atom) } : { bond: a.to.bond.map(atom) },
  }));
}

// The pieces of a spec by its bonds (not the loose ties of ions), each as a formula with its charge:
// 'CH4O', 'Cl-', 'Na+'. For tests and for the status line after a move.
export function pieceFormulas(spec, node) {
  const n = spec.atoms.length;
  const parent = spec.atoms.map((_, i) => i);
  const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  const loose = i => (node ? !!node.fixed[i] : spec.atoms[i].ion === -1 || ['Li', 'Na', 'K', 'Mg', 'Ca', 'Be'].includes(spec.atoms[i].el));
  spec.adj.forEach(([i, j]) => { if (!loose(i) && !loose(j)) parent[find(i)] = find(j); });
  const groups = new Map();
  for (let i = 0; i < n; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
  const order = ['C', 'H'];
  return [...groups.values()].map(members => {
    const counts = {};
    members.forEach(i => { counts[spec.atoms[i].el] = (counts[spec.atoms[i].el] || 0) + 1; });
    const els = Object.keys(counts).sort((a, b) => {
      const ia = order.indexOf(a), ib = order.indexOf(b);
      if (ia >= 0 || ib >= 0) return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
      return a < b ? -1 : 1;
    });
    const q = members.reduce((s, i) => s + (spec.atoms[i].written || 0), 0);
    const mark = q === 0 ? '' : q === 1 ? '+' : q === -1 ? '-' : q > 0 ? `${q}+` : `${-q}-`;
    return els.map(e => e + (counts[e] > 1 ? counts[e] : '')).join('') + mark;
  }).sort();
}
