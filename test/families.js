// node test/families.js   families of equivalent structures for every case
import { CASES, runCase } from '../core/cases.js';
import { structureFamilies } from '../core/symmetry.js';
import { buildNode, enumerate } from '../core/lewis.js';

const extra = {
  'phenol': { atoms: 'C C C C C C O H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [6, 7], [1, 8], [2, 9], [3, 10], [4, 11], [5, 12]], charge: 0 },
  'cyclohexane': { atoms: 'C C C C C C H H H H H H H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [0, 7], [1, 8], [1, 9], [2, 10], [2, 11], [3, 12], [3, 13], [4, 14], [4, 15], [5, 16], [5, 17]], charge: 0 },
  'pentadienyl cation': { atoms: 'C C C C C H H H H H H H'.split(' ').map(el => ({ el })), adj: [[0, 1], [1, 2], [2, 3], [3, 4], [0, 5], [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [4, 11]], charge: 1 },
};
const expected = { 'SF6': [15], 'POCl3': [1], 'acetate': [2], 'allyl cation': [2], 'benzene': [2], 'nitrate with Na+': [3], 'formamide': [1, 1], 'carbon dioxide': [1, 2], 'nitrous oxide': [1, 1], 'SN2 transition: HO- + CH3Cl': [1, 1], 'phenol': [2, 2, 1], 'pentadienyl cation': [2, 1] };
let bad = 0;
const all = [...CASES.map(c => [c.name, () => runCase(c)]), ...Object.entries(extra).map(([name, c]) => [name, () => { const node = buildNode(c); return { node, r: enumerate(node) }; }])];
for (const [name, run] of all) {
  const { node, r } = run();
  const t0 = performance.now();
  const f = structureFamilies(node, r);
  const sizes = f.families.map(x => x.length).sort((a, b) => b - a);
  const want = expected[name];
  const ok = !want || JSON.stringify(sizes) === JSON.stringify([...want].sort((a, b) => b - a));
  if (!ok) bad++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}: ${f.families.length} famil${f.families.length === 1 ? 'y' : 'ies'} of sizes [${sizes.join(' ')}], ${f.symmetries} symmetries${f.truncated ? ' (truncated)' : ''}, ${(performance.now() - t0).toFixed(1)} ms`);
}
console.log(bad ? `${bad} problem(s)` : 'all families as expected');
process.exit(bad ? 1 : 0);
