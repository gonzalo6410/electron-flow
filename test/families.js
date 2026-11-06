// node test/families.js   families of equivalent structures for every case, and for molecules
// chosen to trip the grouping: two identical substituents on a ring (o-, m- and p-xylene differ),
// chains with methyl groups (the hydrogen folding), spectator ions beside a symmetric anion, and
// skeletons with far more than 5000 symmetries (S2F10 has 28800), where an earlier version that
// listed only the first 5000 split one family of 100 into two of 50.
import { CASES, runCase } from '../core/cases.js';
import { structureFamilies, skeletonSymmetries } from '../core/symmetry.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { nodeFromSmiles } from '../core/smiles.js';

const extra = {
  'phenol': { atoms: 'C C C C C C O H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [6, 7], [1, 8], [2, 9], [3, 10], [4, 11], [5, 12]], charge: 0 },
  'cyclohexane': { atoms: 'C C C C C C H H H H H H H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [0, 7], [1, 8], [1, 9], [2, 10], [2, 11], [3, 12], [3, 13], [4, 14], [4, 15], [5, 16], [5, 17]], charge: 0 },
  'pentadienyl cation': { atoms: 'C C C C C H H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [0, 5], [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [4, 11]], charge: 1 },
};
const smiles = {
  'o-xylene': 'Cc1ccccc1C',                 // the two Kekule pictures differ: one has the double bond between the substituted carbons
  'm-xylene': 'Cc1cccc(C)c1',               // a mirror carries one Kekule picture onto the other
  'p-xylene': 'Cc1ccc(C)cc1',
  'hydroquinone': 'Oc1ccc(O)cc1',           // 2 Kekule, 4 ortho and 2 para charge-separated, 2 doubly separated
  'terephthalate (2 Na+)': '[O-]C(=O)c1ccc(C(=O)[O-])cc1.[Na+].[Na+]',
  '2,2-dimethylbutane': 'CC(C)(C)CC',       // the methyl hydrogens fold into C/3 labels
  'dodecane': 'CCCCCCCCCCCC',
  'neopentane': 'CC(C)(C)C',
  'oxalate (2 Na+)': '[O-]C(=O)C(=O)[O-].[Na+].[Na+]',
  'two acetates with Mg2+': 'CC(=O)[O-].CC(=O)[O-].[Mg+2]',
  'triphenylmethyl cation (Cl-)': '[C+](c1ccccc1)(c2ccccc2)c3ccccc3.[Cl-]',
  'IF7': 'FI(F)(F)(F)(F)(F)F',              // 5040 symmetries, 35 structures (three F- cut), one family
  'S2F10': 'FS(F)(F)(F)(F)S(F)(F)(F)(F)F',  // 28800 symmetries; the structures with the S-S pair cut are one family of 100
};
const expected = {
  'SF6': [15], 'POCl3': [1], 'acetate': [2], 'allyl cation': [2], 'benzene': [2], 'nitrate with Na+': [3], 'formamide': [1, 1],
  'carbon dioxide': [1, 2], 'nitrous oxide': [1, 1], 'SN2 transition: HO- + CH3Cl': [1, 1], 'phenol': [2, 2, 1], 'pentadienyl cation': [2, 1],
  'o-xylene': [1, 1], 'm-xylene': [2], 'p-xylene': [2], 'hydroquinone': [4, 2, 2, 2], 'terephthalate (2 Na+)': [8],
  '2,2-dimethylbutane': [1], 'dodecane': [1], 'neopentane': [1], 'oxalate (2 Na+)': [4], 'two acetates with Mg2+': [4],
  'triphenylmethyl cation (Cl-)': [24, 12, 8], 'IF7': [35], 'S2F10': [100, 100],
};
const expectedSymmetries = { 'SF6': 720, 'benzene': 12, 'o-xylene': 2, 'm-xylene': 2, 'p-xylene': 4, 'neopentane': 24, 'dodecane': 2, 'cyclohexane': 12 };
let bad = 0;
const all = [
  ...CASES.map(c => [c.name, () => runCase(c)]),
  ...Object.entries(extra).map(([name, c]) => [name, () => { const node = buildNode(c); return { node, r: enumerate(node) }; }]),
  ...Object.entries(smiles).map(([name, s]) => [name, () => { const node = buildNode(nodeFromSmiles(s)); return { node, r: enumerate(node) }; }]),
];
for (const [name, run] of all) {
  const { node, r } = run();
  const t0 = performance.now();
  const f = structureFamilies(node, r);
  const ms = (performance.now() - t0).toFixed(1);
  const sym = skeletonSymmetries(node, r);
  const sizes = f.families.map(x => x.length).sort((a, b) => b - a);
  const want = expected[name];
  let ok = !want || JSON.stringify(sizes) === JSON.stringify([...want].sort((a, b) => b - a));
  if (expectedSymmetries[name] !== undefined && sym.symmetries !== expectedSymmetries[name]) ok = false;
  if (!ok) bad++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}: ${f.families.length} famil${f.families.length === 1 ? 'y' : 'ies'} of sizes [${sizes.join(' ')}], ${sym.symmetries}${sym.truncated ? '+' : ''} symmetries, ${ms} ms`);
}
console.log(bad ? `${bad} problem(s)` : 'all families as expected');
process.exit(bad ? 1 : 0);
