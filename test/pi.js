// node test/pi.js   checks pi normals and lone pair placement on a few molecules
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';
import { piFrame } from '../core/pi.js';
import { CASES, runCase } from '../core/cases.js';
const deg = x => (Math.acos(Math.max(-1, Math.min(1, x))) * 180 / Math.PI).toFixed(0);
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const M = (els, adj, charge) => ({ atoms: els.split(' ').map(el => ({ el })), adj, charge });
const extra = {
  'allene': M('C C C H H H H', [[0, 1], [1, 2], [0, 3], [0, 4], [2, 5], [2, 6]], 0),
  'ketene': M('C C O H H', [[0, 1], [1, 2], [0, 3], [0, 4]], 0),
};
let bad = 0;
function check(name, node, r, mode) {
  const geo = relax(node, r);
  const pic = mode === 'avg' ? { orders: r.avg.orders, lone: r.avg.lone } : { orders: r.structures[0].orders, lone: r.structures[0].lone };
  const f = piFrame(node, r, geo.pos, geo.groups, pic);
  const sym = node.atoms.map(a => a.el);
  const lines = [];
  r.edges.forEach(([i, j], k) => { if (f.normals[k].length) lines.push(`${sym[i]}${i}-${sym[j]}${j} pi x${f.normals[k].length}`); });
  f.lonePairs.forEach((lp, i) => { if (lp.length) lines.push(`${sym[i]}${i} lp ${lp.map(o => o.kind).join(',')}`); });
  console.log(`== ${name} (${mode}): ` + lines.join('; '));
  return { geo, f, sym };
}
// allene: the two pi normals perpendicular, H planes perpendicular
{
  const node = buildNode(extra.allene); const r = enumerate(node);
  const { f, geo } = check('allene', node, r, 'struct');
  const n1 = f.normals[0][0], n2 = f.normals[1][0];
  const ang = deg(dot(n1, n2));
  const p = geo.pos; const sub = (a, b) => a.map((x, q) => x - b[q]); const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; const ln = a => Math.hypot(...a);
  const pa = cr(sub(p[3], p[0]), sub(p[4], p[0])), pb = cr(sub(p[5], p[2]), sub(p[6], p[2]));
  const twist = deg(dot(pa, pb) / (ln(pa) * ln(pb)));
  console.log(`   pi normals at ${ang} degrees; CH2 planes at ${twist} degrees`);
  if (ang !== '90' || (twist !== '90')) bad++;
}
{
  const node = buildNode(extra.ketene); const r = enumerate(node);
  const { f } = check('ketene', node, r, 'struct');
  if (deg(dot(f.normals[0][0], f.normals[1][0])) !== '90') bad++;
}
for (const name of ['nitrous oxide', 'carbon dioxide', 'formamide', 'acetate', 'SF6', 'SN2 transition: HO- + CH3Cl', 'water', 'nitrate with Na+']) {
  const c = CASES.find(x => x.name === name);
  const { node, r } = runCase(c);
  check(name, node, r, 'struct');
  const { f, geo, sym } = check(name, node, r, 'avg');
  if (name === 'formamide') {
    // O lone pairs in the plane: at 90 degrees from the pi normal; N pair along the p orbital
    const nCO = f.normals[0][0];
    const okO = f.lonePairs[1].slice(0, 2).every(o => deg(dot(o.dir, nCO)) === '90');
    console.log(`   O slot pairs at 90 from pi normal: ${okO}; O third pair kind ${f.lonePairs[1][2] && f.lonePairs[1][2].kind}`);
    if (!okO) bad++;
  }
  if (name === 'carbon dioxide') {
    const ang = deg(dot(f.normals[0][0], f.normals[1][0]));
    console.log(`   CO2 pi normals at ${ang} degrees`);
    if (ang !== '90') bad++;
  }
}
// a triple bond with no flat atom at either end lines up with the pi system across the single bond:
// vinylacetylene's C#C has one pi part parallel to the C=C pi part, benzonitrile's C#N one parallel
// to the ring's; before this rule the three lines of such a triple bond were drawn at a random angle
{
  const { nodeFromSmiles } = await import('../core/smiles.js');
  for (const [name, s] of [['vinylacetylene', 'C=CC#C'], ['benzonitrile', 'N#Cc1ccccc1'], ['diacetylene', 'C#CC#C']]) {
    const node = buildNode(nodeFromSmiles(s)); const r = enumerate(node);
    const st = r.structures.find(x => x.charge.every(q => q === 0)) || r.structures[0];
    const geo = relax(node, r);
    const f = piFrame(node, r, geo.pos, geo.groups, { orders: st.orders, lone: st.lone });
    const multi = r.edges.map((e, k) => k).filter(k => st.orders[k] > 1);
    let widest = 0; // largest angle between the first normals of any two multiple bonds
    multi.forEach((a, q) => multi.slice(q + 1).forEach(b => { widest = Math.max(widest, Number(deg(Math.abs(dot(f.normals[a][0], f.normals[b][0]))))); }));
    const ok = widest < 1;
    console.log(`== ${name}: ${multi.length} multiple bonds, first normals at most ${widest} degrees apart`);
    if (!ok) bad++;
  }
}
console.log(bad ? `${bad} problem(s)` : 'pi frame checks pass');
process.exit(bad ? 1 : 0);
