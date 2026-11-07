// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/rings.js   the ring rule (core/lewis.js, 7 Nov 2026): two neighbouring ring atoms may share
// one bond in the plane of a small ring (benzyne's third bond); an atom with two pi parts and no such
// partner is dropped (pyridine's C=N(+)=C; a phenyl anion whose ring would take its lone pair's place)
// unless no structure has a partner (the ring allenes, which are drawn); electrons in the ring's plane
// stay put in every structure (pyridazine keeps its two nitrogen lone pairs); and the shapes that
// follow (benzyne flat, 1,2-cyclohexadiene drawn).
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';

let bad = 0;
const say = (ok, text) => { if (!ok) bad++; console.log(`${ok ? 'ok  ' : 'FAIL'} ${text}`); };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = a => Math.hypot(a[0], a[1], a[2]);

function run(smi) {
  const node = buildNode(nodeFromSmiles(smi));
  try { return { node, r: enumerate(node) }; } catch (e) { return { node, error: e.message }; }
}

// structure counts
const counts = [
  ['benzyne', 'C1=CC#CC=C1', 2, 'the triple bond and its Kekule partner, where the in-plane bond is the only extra bond between the two carbons'],
  ['cyclohexyne', 'C1CCC#CC1', 1, 'one triple bond in the ring'],
  ['2,3-pyridyne', 'c1cc#cnc1', 2, 'like benzyne'],
  ['1,2-naphthalyne', 'C1=CC#CC2=CC=CC=C12', 3, 'like naphthalene'],
  ['pyridine', 'c1ccncc1', 2, 'C=N(+)=C still refused: the N keeps its lone pair in the plane'],
  ['pyridazine', 'c1ccnnc1', 2, 'no N(+)#N(+) form: the nitrogens keep their in-plane lone pairs'],
  ['1,2,3-triazine', 'c1cnnnc1', 2, 'likewise'],
  ['benzene', 'c1ccccc1', 2, 'unchanged'],
  ['naphthalene', 'c1ccc2ccccc2c1', 3, 'unchanged'],
  ['furan', 'c1ccoc1', 5, 'unchanged: the O keeps one lone pair in the plane in all five'],
  ['1,2-cyclohexadiene', 'C1CC=C=CC1', 1, 'its own only structure: no partner anywhere, so the middle carbon keeps its second bond'],
  ['1,2-cyclopentadiene', 'C1CC=C=C1', 1, 'likewise'],
  ['1,2-cycloheptadiene', 'C1CCC=C=CC1', 1, 'likewise'],
  ['1,2,3-cyclohexatriene', 'C1CC=C=C=C1', 1, 'likewise'],
];
for (const [name, smi, want, why] of counts) {
  const { r, error } = run(smi);
  const got = r ? r.structures.length : error;
  say(got === want, `${name}: ${got} structures (want ${want}: ${why})`);
}
// the phenyl anion keeps its lone pair on the bare carbon, in the plane (dropping the pairing part
// outright moved it out of the plane onto a ring CH and broke the third benzyne scene)
{
  const { node, r } = run('[c-]1ccccc1.[Na+]');
  const bare = node.atoms.findIndex((a, i) => a.el === 'C' && !node.edges.some(e => e.includes(i) && node.atoms[e[0] === i ? e[1] : e[0]].el === 'H'));
  say(r.structures.length === 2 && r.structures.every(st => st.charge[bare] === -1 && st.lone[bare] === 2),
    `phenyl anion: ${r.structures.length} structures, the minus and the lone pair on the bare carbon in each`);
}

// 1,2-cyclohexadiene drawn: the ring closes, its bonds near their lengths (precision is the geometric
// view's job, item 67)
{
  const { node, r } = run('C1CC=C=CC1');
  const geo = relax(node, r);
  const p = geo.pos;
  const ds = node.edges.map(([i, j], k) => ({ d: len(sub(p[i], p[j])), o: r.avg.orders[k], h: node.atoms[i].el === 'H' || node.atoms[j].el === 'H' }))
    .filter(e => !e.h);
  const ok = ds.every(e => (e.o === 2 ? e.d > 1.25 && e.d < 1.4 : e.d > 1.45 && e.d < 1.6));
  say(ok, `1,2-cyclohexadiene drawn: ring bonds ${ds.map(e => e.d.toFixed(2)).join(' ')} A`);
}

// benzyne's shape: flat, the in-plane bond between a double and a triple bond in length
{
  const { node, r } = run('C1=CC#CC=C1');
  const geo = relax(node, r);
  const p = geo.pos;
  const heavy = node.atoms.map((a, i) => i).filter(i => node.atoms[i].el !== 'H');
  const c = [0, 1, 2].map(k => heavy.reduce((s, i) => s + p[i][k], 0) / heavy.length);
  const nrm = cross(sub(p[heavy[0]], c), sub(p[heavy[2]], c));
  const u = nrm.map(v => v / len(nrm));
  const off = Math.max(...node.atoms.map((a, i) => Math.abs(dot(sub(p[i], c), u))));
  say(off < 0.02, `benzyne flat: every atom within ${off.toFixed(3)} A of the ring plane`);
  const k = node.edges.findIndex(([i, j], q) => Math.abs(r.avg.orders[q] - 2.5) < 1e-9);
  const d = len(sub(p[node.edges[k][0]], p[node.edges[k][1]]));
  say(d > 1.2 && d < 1.3, `benzyne's in-plane bond (order 2.5) ${d.toFixed(3)} A, between triple and double`);
  const angs = [];
  for (const i of heavy) {
    const nb = node.edges.filter(e => e.includes(i)).map(e => (e[0] === i ? e[1] : e[0])).filter(j => node.atoms[j].el !== 'H');
    const a = sub(p[nb[0]], p[i]), b = sub(p[nb[1]], p[i]);
    angs.push(Math.acos(dot(a, b) / (len(a) * len(b))) * 180 / Math.PI);
  }
  say(angs.every(x => x > 115 && x < 125), `benzyne ring angles ${angs.map(x => x.toFixed(1)).join(' ')} (the ring's 120, give or take)`);
}

if (bad) { console.log(`${bad} failed`); process.exit(1); }
console.log('ring rule checks pass');
