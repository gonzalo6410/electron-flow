// node test/geometry.js   prints rule-based geometry for every case
import { relax, geometryReport } from '../core/geometry.js';
import { CASES, runCase } from '../core/cases.js';

for (const c of CASES) {
  const { node, r } = runCase(c);
  const t0 = performance.now();
  const geo = relax(node, r);
  console.log(`\n== ${c.name} (${(performance.now() - t0).toFixed(0)} ms)`);
  console.log(geometryReport(node, r, geo));
}
