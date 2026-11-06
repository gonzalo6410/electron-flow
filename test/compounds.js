// node test/compounds.js   every named compound parses, enumerates and relaxes; prints a one-line summary
// with the time of the whole run and of the geometry relaxation alone, then the slowest relaxation.
import { COMPOUNDS } from '../core/compounds.js';
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';
import { structureFamilies } from '../core/symmetry.js';
let bad = 0;
let relaxTotal = 0, slowest = { name: '', ms: 0 };
for (const c of COMPOUNDS) {
  const t0 = performance.now();
  try {
    const spec = nodeFromSmiles(c.smiles);
    const node = buildNode(spec);
    const r = enumerate(node);
    const t1 = performance.now();
    const geo = relax(node, r);
    const relaxMs = performance.now() - t1;
    relaxTotal += relaxMs;
    if (relaxMs > slowest.ms) slowest = { name: c.name, ms: relaxMs };
    const fam = structureFamilies(node, r);
    const sym = node.atoms.map(a => a.el);
    const charges = r.avg.charge.map((q, i) => Math.abs(q) > 1e-9 ? `${sym[i]}${i}:${q.toFixed(2)}` : null).filter(Boolean).join(' ');
    console.log(`ok   ${c.name}: ${node.atoms.length} atoms, ${r.structures.length} structures in ${fam.families.length} families, holes ${r.holes / 2}${r.radicals ? ', radical' : ''}, geometry energy ${geo.energy.toFixed(2)} in ${geo.iterations} steps, ${(performance.now() - t0).toFixed(0)} ms (relax ${relaxMs.toFixed(1)} ms)${charges ? ' | ' + charges : ''}`);
  } catch (e) {
    bad++;
    console.log(`FAIL ${c.name}: ${e.message}`);
  }
}
console.log(`relaxation: ${relaxTotal.toFixed(0)} ms over ${COMPOUNDS.length} compounds, slowest ${slowest.name} at ${slowest.ms.toFixed(1)} ms`);
console.log(bad ? `${bad} problem(s)` : 'all compounds run');
