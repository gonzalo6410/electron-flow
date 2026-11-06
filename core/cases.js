import { buildNode, enumerate, hybrid } from './lewis.js';

// Known-answer cases from the ledger (item 7 and item 28). Expected values are his statements
// where tagged [STATED]; the rest follow from the rules and are checked by hand in the test output.
// Positions are rough starting guesses only; the geometry module relaxes them.

export const CASES = [
  {
    name: 'SF6',
    atoms: [{ el: 'S' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }],
    adj: [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6]],
    charge: 0,
    expect: { structures: 15, holes: 0, charge: { 0: 2, 1: -1 / 3 }, orders: { 0: 2 / 3 } }, // [STATED]
  },
  {
    name: 'POCl3',
    atoms: [{ el: 'P' }, { el: 'O' }, { el: 'Cl' }, { el: 'Cl' }, { el: 'Cl' }],
    adj: [[0, 1], [0, 2], [0, 3], [0, 4]],
    charge: 0,
    // [STATED 5 Nov 2026] one structure, P+ with four single bonds and O-. A free Cl- is a sigma bond
    // broken, which is another node, not resonance within this one (ledger item 23).
    expect: { structures: 1, holes: 0, charge: { 0: 1, 1: -1 } },
  },
  {
    name: 'acetate',
    atoms: [{ el: 'C' }, { el: 'O' }, { el: 'O' }, { el: 'C' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2], [0, 3], [3, 4], [3, 5], [3, 6]],
    charge: -1,
    expect: { structures: 2, holes: 0, charge: { 1: -0.5, 2: -0.5 }, orders: { 0: 1.5, 1: 1.5 } },
  },
  {
    name: 'allyl cation',
    atoms: [{ el: 'C' }, { el: 'C' }, { el: 'C' }, { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [1, 2], [0, 3], [0, 4], [1, 5], [2, 6], [2, 7]],
    charge: 1,
    expect: { structures: 2, holes: 2, charge: { 0: 0.5, 1: 0, 2: 0.5 }, orders: { 0: 1.5, 1: 1.5 } },
  },
  {
    name: 'benzene',
    atoms: [
      { el: 'C' }, { el: 'C' }, { el: 'C' }, { el: 'C' }, { el: 'C' }, { el: 'C' },
      { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'H' },
    ],
    adj: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11]],
    charge: 0,
    expect: { structures: 2, holes: 0, orders: { 0: 1.5, 1: 1.5 } },
  },
  {
    name: 'water',
    atoms: [{ el: 'O' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2]],
    charge: 0,
    expect: { structures: 1, holes: 0, charge: { 0: 0 } },
  },
  {
    name: 'NaOH',
    atoms: [{ el: 'Na' }, { el: 'O' }, { el: 'H' }],
    adj: [[0, 1], [1, 2]],
    charge: 0,
    expect: { structures: 1, holes: 0, charge: { 0: 1, 1: -1 } },
  },
  {
    name: 'chloromethane',
    atoms: [{ el: 'C' }, { el: 'Cl' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2], [0, 3], [0, 4]],
    charge: 0,
    expect: { structures: 1, holes: 0 },
  },
  {
    name: 'SN2 transition: HO- + CH3Cl',
    // before: HO- (Na+ beside it) and CH3Cl. after: CH3OH and Cl- (Na+ beside it).
    // The hybrid is the union of edges; the forming O-C and breaking C-Cl come out at order 1/2.
    atoms: [{ el: 'O' }, { el: 'H' }, { el: 'C' }, { el: 'Cl' }, { el: 'H' }, { el: 'H' }, { el: 'H' }, { el: 'Na' }],
    before: [[0, 1], [2, 3], [2, 4], [2, 5], [2, 6], [0, 7]],
    after: [[0, 1], [0, 2], [2, 4], [2, 5], [2, 6], [3, 7]],
    charge: 0,
    expect: { structures: 2, holes: 0, charge: { 0: -0.5, 3: -0.5, 2: 0 }, orders: { 'O0-C2': 0.5, 'C2-Cl3': 0.5 }, count: { 2: 8 } }, // [STATED]
  },
  {
    name: 'methyl cation',
    atoms: [{ el: 'C' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2], [0, 3]],
    charge: 1,
    expect: { structures: 1, holes: 2, charge: { 0: 1 } },
  },
  {
    name: 'borane',
    atoms: [{ el: 'B' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2], [0, 3]],
    charge: 0,
    expect: { structures: 1, holes: 2, charge: { 0: 0 } },
  },
  {
    name: 'nitrate with Na+',
    atoms: [{ el: 'N' }, { el: 'O' }, { el: 'O' }, { el: 'O' }, { el: 'Na' }],
    adj: [[0, 1], [0, 2], [0, 3], [1, 4]],
    charge: 0,
    expect: { structures: 3, holes: 0, charge: { 0: 1, 1: -2 / 3 }, orders: { 0: 4 / 3 } },
  },
  {
    name: 'formamide',
    // O=C-N (no charges) and O(-)-C=N(+): both full octets, no cuts, so both contribute equally
    atoms: [{ el: 'C' }, { el: 'O' }, { el: 'N' }, { el: 'H' }, { el: 'H' }, { el: 'H' }],
    adj: [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5]],
    charge: 0,
    expect: { structures: 2, holes: 0, charge: { 1: -0.5, 2: 0.5 }, orders: { 0: 1.5, 1: 1.5 } },
  },
  {
    name: 'carbon dioxide',
    // O=C=O plus the two O(-)-C#O(+) forms: three equal structures under pure formalism
    atoms: [{ el: 'C' }, { el: 'O' }, { el: 'O' }],
    adj: [[0, 1], [0, 2]],
    charge: 0,
    expect: { structures: 3, holes: 0, charge: { 0: 0, 1: 0 }, orders: { 0: 2, 1: 2 } },
  },
  {
    name: 'nitrous oxide',
    // N#N(+)-O(-) and N(-)=N(+)=O; the N(2-)-N(+)#O(+) form is dropped by the magnitude rule
    atoms: [{ el: 'N' }, { el: 'N' }, { el: 'O' }],
    adj: [[0, 1], [1, 2]],
    charge: 0,
    expect: { structures: 2, holes: 0, charge: { 0: -0.5, 1: 1, 2: -0.5 }, orders: { 0: 2.5, 1: 1.5 } },
  },
];

// Build and enumerate a case; a case with before/after is a transition hybrid.
export function runCase(c) {
  if (c.before) {
    const a = buildNode({ atoms: c.atoms, adj: c.before, charge: c.charge });
    const b = buildNode({ atoms: c.atoms, adj: c.after, charge: c.charge });
    const r = hybrid(a, b);
    return { node: { ...a, edges: r.edges }, r };
  }
  const node = buildNode(c);
  return { node, r: enumerate(node) };
}
