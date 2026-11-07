// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/gradient.js   the exact derivative of the spring energy must agree with finite differences
// on several molecules, at the rough starting positions (where every term is pulling) and at a
// stretched copy of them. Also checks that the relaxed energy is the same whether the energy is
// evaluated through the model or recomputed from the relaxed positions.
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax, energyModel } from '../core/geometry.js';
import { CASES, runCase } from '../core/cases.js';

const smiles = {
  ethane: 'CC', butadiene: 'C=CC=C', allene: 'C=C=C', ketene: 'C=C=O', benzene: 'c1ccccc1',
  cyclohexane: 'C1CCCCC1', 'acetate (Na+)': 'CC(=O)[O-].[Na+]', tropylium: '[CH+]1C=CC=CC=C1.[Cl-]',
  SF6: 'FS(F)(F)(F)(F)F', XeF4: 'F[Xe](F)(F)F', water: 'O', pyridine: 'c1ccncc1', 'methyl cation (Cl-)': '[CH3+].[Cl-]',
};
let bad = 0;
function check(name, node, r) {
  const model = energyModel(node, r);
  // a deterministic rough start: the relaxed positions pulled about by up to 0.3 A, then also squashed
  const geo = relax(node, r);
  let seed = 777;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };
  const starts = [
    new Float64Array(geo.pos.flat().map(v => v + 0.6 * rnd())),
    new Float64Array(geo.pos.flat().map((v, q) => v * (q % 3 === 2 ? 0.7 : 1.15) + 0.2 * rnd())),
  ];
  for (const x of starts) {
    const g = new Float64Array(x.length);
    const E = model.gradient(x, g);
    const h = 1e-6;
    let worst = 0, scale = 0;
    for (let q = 0; q < x.length; q++) {
      const keep = x[q];
      x[q] = keep + h; const ep = model.energy(x);
      x[q] = keep - h; const em = model.energy(x);
      x[q] = keep;
      const num = (ep - em) / (2 * h);
      worst = Math.max(worst, Math.abs(num - g[q]));
      scale = Math.max(scale, Math.abs(num));
    }
    const ok = worst <= 1e-6 * Math.max(1, scale);
    if (!ok) bad++;
    console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}: energy ${E.toFixed(3)}, largest slope ${scale.toFixed(3)}, worst disagreement ${worst.toExponential(2)}`);
  }
  // the energy reported by relax is the energy of the positions it returns
  const again = model.energy(new Float64Array(geo.pos.flat()));
  if (Math.abs(again - geo.energy) > 1e-9) { bad++; console.log(`FAIL ${name}: relaxed energy ${geo.energy} but recomputed ${again}`); }
}
for (const [name, smi] of Object.entries(smiles)) {
  const node = buildNode(nodeFromSmiles(smi));
  check(name, node, enumerate(node));
}
for (const name of ['SN2 transition: HO- + CH3Cl', 'nitrous oxide']) {
  const { node, r } = runCase(CASES.find(c => c.name === name));
  check(name, node, r);
}
console.log(bad ? `${bad} problem(s)` : 'exact and finite-difference slopes agree');
process.exit(bad ? 1 : 0);
