// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/shapes.js   known-answer checks on the relaxed shapes: the hypervalent and odd cases,
// the angle rules (lone pairs, ring sizes) and the short polar single bonds.
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';
import { CASES, runCase } from '../core/cases.js';

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = a => Math.hypot(a[0], a[1], a[2]);
const deg = x => Math.acos(Math.max(-1, Math.min(1, x))) * 180 / Math.PI;
const angle = (p, a, i, b) => { const u = sub(p[a], p[i]), v = sub(p[b], p[i]); return deg(dot(u, v) / (len(u) * len(v))); };
const dist = (p, i, j) => len(sub(p[i], p[j]));
// signed torsion a-b-c-d in degrees
function torsion(p, a, b, c, d) {
  const b1 = sub(p[b], p[a]), b2 = sub(p[c], p[b]), b3 = sub(p[d], p[c]);
  const n1 = cross(b1, b2), n2 = cross(b2, b3);
  const m = cross(n1, b2.map(x => x / len(b2)));
  return Math.atan2(dot(m, n2), dot(n1, n2)) * 180 / Math.PI;
}
// largest distance of the listed atoms from their best plane
function outOfPlane(p, atoms) {
  const c = [0, 1, 2].map(k => atoms.reduce((s, i) => s + p[i][k], 0) / atoms.length);
  let worst = 0;
  // try the plane through the centroid normal to each pair of atom directions, keep the best fit
  let best = Infinity, bestN = null;
  for (let x = 0; x < atoms.length; x++) for (let y = x + 1; y < atoms.length; y++) {
    const nrm = cross(sub(p[atoms[x]], c), sub(p[atoms[y]], c));
    if (len(nrm) < 1e-6) continue;
    const u = nrm.map(v => v / len(nrm));
    const off = Math.max(...atoms.map(i => Math.abs(dot(sub(p[i], c), u))));
    if (off < best) { best = off; bestN = u; }
  }
  worst = best;
  return worst;
}

let bad = 0;
function expect(label, value, want, tol) {
  const ok = Math.abs(value - want) <= tol;
  if (!ok) bad++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: ${value.toFixed(2)} (want ${want} within ${tol})`);
}
function fromSmiles(smi) { const node = buildNode(nodeFromSmiles(smi)); const r = enumerate(node); return { node, r, geo: relax(node, r) }; }
function fromCase(name) { const { node, r } = runCase(CASES.find(c => c.name === name)); return { node, r, geo: relax(node, r) }; }
const sortedAngles = (p, i, nbrs) => { const out = []; for (let a = 0; a < nbrs.length; a++) for (let b = a + 1; b < nbrs.length; b++) out.push(angle(p, nbrs[a], i, nbrs[b])); return out.sort((x, y) => x - y); };

// hypervalent shapes
{
  const { geo } = fromSmiles('FS(F)(F)(F)(F)F'); // S is atom 1
  const angs = sortedAngles(geo.pos, 1, [0, 2, 3, 4, 5, 6]);
  expect('SF6 octahedral: twelve angles of 90', angs[11], 90, 0.5);
  expect('SF6 octahedral: three angles of 180', angs[12], 180, 0.5);
}
{
  const { geo } = fromSmiles('ClP(Cl)(Cl)(Cl)Cl'); // P is atom 1
  const angs = sortedAngles(geo.pos, 1, [0, 2, 3, 4, 5]);
  expect('PCl5 trigonal bipyramid: six angles of 90', angs[5], 90, 0.5);
  expect('PCl5 trigonal bipyramid: three angles of 120', angs[6], 120, 0.5);
  expect('PCl5 trigonal bipyramid: three angles of 120', angs[8], 120, 0.5);
  expect('PCl5 trigonal bipyramid: one angle of 180', angs[9], 180, 0.5);
}
{
  const { geo } = fromSmiles('F[Xe](F)(F)F'); // Xe is atom 1
  const angs = sortedAngles(geo.pos, 1, [0, 2, 3, 4]);
  expect('XeF4 square planar: four angles of 90', angs[3], 90, 0.5);
  expect('XeF4 square planar: two angles of 180', angs[4], 180, 0.5);
  expect('XeF4 square planar: all five atoms in one plane', outOfPlane(geo.pos, [0, 1, 2, 3, 4]), 0, 0.02);
}
{
  const { geo } = fromCase('SN2 transition: HO- + CH3Cl'); // O0 H1 C2 Cl3 H4 H5 H6 Na7
  expect('SN2 state: half bonds on the axis, O-C-Cl', angle(geo.pos, 0, 2, 3), 180, 3);
  expect('SN2 state: flat CH3, H-C-H', angle(geo.pos, 4, 2, 5), 120, 1);
  expect('SN2 state: O-C-H at right angles', angle(geo.pos, 0, 2, 4), 90, 3);
  expect('SN2 state: half bond O-C drawn long', dist(geo.pos, 0, 2), 1.775, 0.02);
}
// allene, benzene, cyclohexane, butadiene
{
  const { geo } = fromSmiles('C=C=C'); // C0 C1 C2 H3 H4 H5 H6
  const p = geo.pos;
  const n1 = cross(sub(p[3], p[0]), sub(p[4], p[0])), n2 = cross(sub(p[5], p[2]), sub(p[6], p[2]));
  expect('allene: CH2 planes at right angles', deg(dot(n1, n2) / (len(n1) * len(n2))), 90, 1);
  expect('allene: straight C=C=C', angle(p, 0, 1, 2), 180, 1);
}
{
  const { geo } = fromSmiles('c1ccccc1');
  expect('benzene: flat', outOfPlane(geo.pos, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]), 0, 0.02);
  expect('benzene: C-C-C 120', angle(geo.pos, 0, 1, 2), 120, 0.5);
  expect('benzene: C-C 1.42 (order 1.5 by the length rule)', dist(geo.pos, 0, 1), 1.4212, 0.01);
}
{
  const { geo } = fromSmiles('C1CCCCC1');
  const t = [0, 1, 2, 3, 4, 5].map(k => torsion(geo.pos, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6));
  expect('cyclohexane chair: ring torsions near 55', Math.min(...t.map(Math.abs)), 55, 6);
  expect('cyclohexane chair: ring torsions near 55', Math.max(...t.map(Math.abs)), 55, 6);
  const alternating = t.every((v, k) => Math.sign(v) !== Math.sign(t[(k + 1) % 6]));
  expect('cyclohexane chair: torsion signs alternate', alternating ? 1 : 0, 1, 0);
  expect('cyclohexane: C-C-C 109.5', angle(geo.pos, 0, 1, 2), 109.47, 0.5);
}
{
  const { geo } = fromSmiles('C=CC=C');
  expect('butadiene: flat', outOfPlane(geo.pos, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), 0, 0.02);
  expect('butadiene: s-trans, C-C-C-C torsion 180', Math.abs(torsion(geo.pos, 0, 1, 2, 3)), 180, 1);
}
// lone pair rule
{
  expect('water H-O-H (two pairs)', angle(fromSmiles('O').geo.pos, 1, 0, 2), 104.5, 0.1);
  expect('ammonia H-N-H (one pair)', angle(fromSmiles('N').geo.pos, 1, 0, 2), 107, 0.1);
  expect('methanol C-O-H (two pairs)', angle(fromSmiles('CO').geo.pos, 0, 1, 5), 104.5, 0.1);
  expect('sulfur dioxide O-S-O (one pair, flat)', angle(fromSmiles('O=S=O').geo.pos, 0, 1, 2), 117.5, 0.1);
  const amide = fromSmiles('C(=O)N'); // C0 O1 N2 H3 H4 H5: the N pair is in the p orbital, no narrowing, flat
  expect('formamide H-N-H stays 120', angle(amide.geo.pos, 4, 2, 5), 120, 0.5);
  expect('formamide: flat', outOfPlane(amide.geo.pos, [0, 1, 2, 3, 4, 5]), 0, 0.02);
  expect('ammonium H-N-H unchanged (no pairs)', angle(fromSmiles('[NH4+].[Cl-]').geo.pos, 2, 0, 3), 109.47, 1.5);
}
// ring rule
{
  const cp = fromSmiles('[CH-]1C=CC=C1.[Na+]'); // C0..C4 ring, Na5, H6..
  const ring = [0, 1, 2, 3, 4];
  expect('cyclopentadienyl anion: flat ring', outOfPlane(cp.geo.pos, ring), 0, 0.03);
  expect('cyclopentadienyl anion: ring angles 108', angle(cp.geo.pos, 1, 2, 3), 108, 1);
  expect('cyclopentadienyl anion: C-C-H 126', angle(cp.geo.pos, 1, 2, 8), 126, 1);
  const trop = fromSmiles('[CH+]1C=CC=CC=C1.[Cl-]'); // C0..C6 ring
  expect('tropylium: flat ring', outOfPlane(trop.geo.pos, [0, 1, 2, 3, 4, 5, 6]), 0, 0.03);
  expect('tropylium: ring angles 128.6', angle(trop.geo.pos, 1, 2, 3), 128.57, 1);
  expect('tropylium: low strain', trop.geo.energy, 0, 1);
  const tropone = fromSmiles('O=C1C=CC=CC=C1');
  expect('tropone: flat ring', outOfPlane(tropone.geo.pos, [1, 2, 3, 4, 5, 6, 7]), 0, 0.03);
  const pyrrole = fromSmiles('c1cc[nH]c1'); // C0 C1 C2 N3 C4
  expect('pyrrole: ring angle at N near 108', angle(pyrrole.geo.pos, 2, 3, 4), 108, 2);
  const cyclopentane = fromSmiles('C1CCCC1');
  expect('cyclopentane: ring angle 108 (below tetrahedral)', angle(cyclopentane.geo.pos, 0, 1, 2), 108, 0.5);
  const cyclohexene = fromSmiles('C1CCC=CC1'); // six-ring: no change from the templates
  expect('cyclohexene: C=C-C keeps near 120 (six-ring, no ring rule)', angle(cyclohexene.geo.pos, 3, 4, 5), 120, 3);
}
// short polar single bonds from a charged third-row atom
{
  const sulfate = fromSmiles('[O-]S(=O)(=O)[O-].[Na+].[Na+]'); // O0 S1 O2 O3 O4
  expect('sulfate S-O at the double-bond length', dist(sulfate.geo.pos, 1, 2), 1.4877, 0.02);
  const dmso = fromSmiles('CS(=O)C'); // C0 S1 O2 C3
  expect('DMSO S-O at the double-bond length', dist(dmso.geo.pos, 1, 2), 1.4877, 0.02);
  expect('DMSO C-S stays single', dist(dmso.geo.pos, 0, 1), 1.81, 0.02);
  const perchlorate = fromSmiles('[O-]Cl(=O)(=O)=O.[Na+]'); // O0 Cl1
  expect('perchlorate Cl-O at the double-bond length', dist(perchlorate.geo.pos, 0, 1), 1.4616, 0.02);
  const ab = fromSmiles('[NH3+][BH3-]'); // second row: the rule does not apply
  expect('ammonia borane N-B stays single', dist(ab.geo.pos, 0, 1), 1.55, 0.02);
  const nitromethane = fromSmiles('C[N+](=O)[O-]'); // N1 O2: order 1.5, not a single bond
  expect('nitromethane N-O by the order rule', dist(nitromethane.geo.pos, 1, 2), 1.2808, 0.02);
}
// loose ions: at their ionic radius plus the partner's covalent radius plus 0.7, never as close as a
// bond, and started where they do not sit on a hydrogen: in the plane of a flat cation between two
// of its bonds, through a face of a tetrahedral one
{
  const mecl = fromSmiles('[CH3+].[Cl-]'); // C0 Cl1 H2 H3 H4
  expect('methyl cation: Cl- sits 3.27 from the carbon', dist(mecl.geo.pos, 0, 1), 3.27, 0.05);
  expect('methyl cation: Cl- in the plane of the CH3 (out of plane)', outOfPlane(mecl.geo.pos, [0, 1, 2, 3, 4]), 0, 0.05);
  expect('methyl cation: Cl- at least 2.5 from every hydrogen', Math.min(...[2, 3, 4].map(h => dist(mecl.geo.pos, 1, h))), 3.0, 0.5);
  const nh4cl = fromSmiles('[NH4+].[Cl-]'); // N0 Cl1 H2..H5
  expect('ammonium: Cl- sits 3.2 from the nitrogen', dist(nh4cl.geo.pos, 0, 1), 3.2, 0.1);
  expect('ammonium: Cl- through a face, three hydrogens at 3.0 and one beyond 4', [2, 3, 4, 5].map(h => dist(nh4cl.geo.pos, 1, h)).sort((a, b) => a - b)[3], 4.2, 0.2);
  const nacl = fromSmiles('[Na+].[Cl-]');
  expect('sodium chloride: the two ions 2.83 apart', dist(nacl.geo.pos, 0, 1), 2.83, 0.05);
  const acetate = fromSmiles('CC(=O)[O-].[Na+]'); // C0 C1 O2 O3 Na4
  expect('sodium acetate: Na+ between 2.4 and 3.3 from its oxygen (the pushes decide)', dist(acetate.geo.pos, 3, 4), 2.85, 0.45);
}
console.log(bad ? `${bad} problem(s)` : 'all shapes as expected');
process.exit(bad ? 1 : 0);
