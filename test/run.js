// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/run.js [--all]   (run from the electron-flow folder)
import { describe } from '../core/lewis.js';
import { CASES, runCase } from '../core/cases.js';

const EPS = 1e-9;
let failures = 0;

for (const c of CASES) {
  const t0 = performance.now();
  const { node, r } = runCase(c);
  const ms = (performance.now() - t0).toFixed(1);
  const sym = node.atoms.map(a => a.el);
  const edgeName = ([i, j]) => `${sym[i]}${i}-${sym[j]}${j}`;
  const problems = [];
  const e = c.expect || {};
  if (e.structures !== undefined && r.structures.length !== e.structures) problems.push(`structures ${r.structures.length}, expected ${e.structures}`);
  if (e.holes !== undefined && r.holes !== e.holes) problems.push(`holes ${r.holes}, expected ${e.holes}`);
  for (const [i, v] of Object.entries(e.charge || {})) if (Math.abs(r.avg.charge[i] - v) > EPS) problems.push(`charge on atom ${i} ${r.avg.charge[i]}, expected ${v}`);
  for (const [k, v] of Object.entries(e.orders || {})) {
    const idx = /^\d+$/.test(k) ? Number(k) : r.edges.findIndex(ed => edgeName(ed) === k);
    if (idx < 0) { problems.push(`no edge ${k}`); continue; }
    if (Math.abs(r.avg.orders[idx] - v) > EPS) problems.push(`order on edge ${k} ${r.avg.orders[idx]}, expected ${v}`);
  }
  for (const [i, v] of Object.entries(e.count || {})) if (Math.abs(r.avg.count[i] - v) > EPS) problems.push(`count on atom ${i} ${r.avg.count[i]}, expected ${v}`);
  console.log(`\n== ${c.name} (${ms} ms) ${problems.length ? 'FAIL' : 'ok'}`);
  console.log(describe(node, r, { all: process.argv.includes('--all') }));
  if (c.hisStatement) console.log('  NOTE: engine default differs from his statement (' + c.hisStatement + '); open call');
  problems.forEach(p => console.log('  PROBLEM: ' + p));
  failures += problems.length;
}
console.log(`\n${failures ? failures + ' problem(s)' : 'all cases pass'}`);
process.exit(failures ? 1 : 0);
