// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/film.js   the shapes of the film of every reaction scene (core/moves.js, filmSteps): the
// substrate does not swing about (the first round of the moves review saw chloromethane turn about 90
// degrees at the start), no atom jumps between two steps, the SN2 film keeps O, C and Cl on a line,
// and in the lineup before the film every other molecule keeps its own shape and passes no atom of
// another molecule closer than it will be in the approach.
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';
import { commit, filmSteps, biggestPiece, filmCore } from '../core/moves.js';
import { bestFit, rotate, add } from '../core/vectors.js';
import { REACTIONS, resolveArrows } from '../core/reactions.js';

let bad = 0;
const say = (ok, text) => { if (!ok) bad++; console.log(`${ok ? 'ok  ' : 'FAIL'} ${text}`); };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = a => Math.hypot(a[0], a[1], a[2]);
const angle = (p, a, i, b) => { const u = sub(p[a], p[i]), v = sub(p[b], p[i]); return Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (len(u) * len(v))))) * 180 / Math.PI; };

// the turn of a set of atoms as a whole between two geometries: the largest angle any of their
// directions from their centre turns through (atoms near the centre left out)
function turnOf(atoms, p, q) {
  const c = P => [0, 1, 2].map(k => atoms.reduce((s, i) => s + P[i][k], 0) / atoms.length);
  const cp = c(p), cq = c(q);
  let worst = 0;
  for (const i of atoms) {
    const u = sub(p[i], cp), v = sub(q[i], cq);
    if (len(u) < 0.5 || len(v) < 0.5) continue;
    worst = Math.max(worst, Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (len(u) * len(v))))) * 180 / Math.PI);
  }
  return worst;
}

for (const rx of REACTIONS) {
  const spec = nodeFromSmiles(rx.smiles), node = buildNode(spec), r = enumerate(node);
  const arrows = resolveArrows(spec, rx.arrows);
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  let out = null;
  for (const st of [...r.structures].sort((a, b) => spread(a) - spread(b))) { out = commit(node, st, arrows, r); if (out.kind === 'move') break; }
  const geo = relax(node, r);
  const { node: fnode, steps, lineup } = filmSteps(node, out, geo, 21);
  const heavy = biggestPiece(node).filter(i => node.atoms[i].el !== 'H' && !node.fixed[i]);
  // printed, not judged: a molecule may have to change its own shape (butadiene turning s-cis,
  // acetaldehyde's C=O coming into line with the forming enolate)
  const approachTurn = heavy.length >= 2 ? turnOf(heavy, steps[0].geo.pos, steps[lineup + 1].geo.pos) : 0;
  // the biggest jump of an atom between two steps of the blend (the settle at the end is left out:
  // a loose ion may drift a long way to its new partner there)
  let jump = 0, jumpWhere = '';
  for (let q = 1; q < steps.length - 2; q++) {
    steps[q].geo.pos.forEach((p, i) => {
      const x = len(sub(p, steps[q + 1].geo.pos[i]));
      if (!fnode.fixed[i] && x > jump) { jump = x; jumpWhere = `${node.atoms[i].el}${i + 1} at step ${q} of ${steps.length - 1}, lineup ends at ${lineup}`; }
    });
  }
  console.log(`     ${rx.name}: the heavy atoms of the biggest molecule turn ${approachTurn.toFixed(0)} degrees by the first blended frame`);
  say(jump < 1.2, `   no atom moves more than ${jump.toFixed(2)} A between two steps of the lineup and the blend (under 1.2; ${jumpWhere})`);
  // the substrate stays where it was on screen: its heavy atoms that stay together drift little from
  // where they started (the enolate scene's acetaldehyde drifted 1.9 A and turned on screen before
  // each step was laid onto the one before; moves review, 7 Nov 2026)
  const core = filmCore(node, out.after).filter(i => node.atoms[i].el !== 'H');
  let drift = 0;
  for (const st of steps) for (const i of core) drift = Math.max(drift, len(sub(st.geo.pos[i], steps[0].geo.pos[i])));
  say(drift < 1.0, `   the substrate's own heavy atoms drift at most ${drift.toFixed(2)} A from where they were (under 1.0; the benzyne ring opens its angles by itself)`);
  // the lineup: no bond is ever shorter or longer than it is on screen and in the approach (by more
  // than 0.1 A), a loose ion's tie by no more than 1.5 A (it may go round wide to clear), and no two atoms of different molecules come closer than the closer of where
  // they start and where the approach puts them, less 0.2 A (moving atoms straight from start to end
  // shrank a turning hydroxide and could pass one molecule through another)
  {
    const group = (() => {
      const parent = node.atoms.map((_, i) => i);
      const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
      node.edges.forEach(([i, j]) => { parent[find(i)] = find(j); });
      return node.atoms.map((_, i) => find(i));
    })();
    const d = (P, i, j) => len(sub(P[i], P[j]));
    const A = steps[0].geo.pos, Z = steps[lineup].geo.pos;
    let stretch = 0, tieStretch = 0, worstBond = '', close = Infinity, closeWhere = '';
    for (let q = 1; q <= lineup; q++) {
      const P = steps[q].geo.pos;
      node.edges.forEach(([i, j]) => {
        const x = d(P, i, j), lo = Math.min(d(A, i, j), d(Z, i, j)), hi = Math.max(d(A, i, j), d(Z, i, j));
        if (node.fixed[i] || node.fixed[j]) { tieStretch = Math.max(tieStretch, lo - x, x - hi); return; }
        if (Math.max(lo - x, x - hi) > stretch) { stretch = Math.max(lo - x, x - hi); worstBond = `${node.atoms[i].el}${i + 1}-${node.atoms[j].el}${j + 1} ${lo.toFixed(2)}..${hi.toFixed(2)} got ${x.toFixed(2)} at ${q}${node.fixed[i] || node.fixed[j] ? ' tie' : ''}`; }
      });
      for (let i = 0; i < node.atoms.length; i++) for (let j = i + 1; j < node.atoms.length; j++) {
        if (group[i] === group[j]) continue;
        // a loose ion only has to keep 2.2 A clear (or as close as it is at either end)
        const ion = node.fixed[i] || node.fixed[j];
        const end = Math.min(d(A, i, j), d(Z, i, j));
        const gap = ion ? d(P, i, j) - Math.min(2.2, end) + 0.2 : d(P, i, j) - end;
        if (gap < close) { close = gap; closeWhere = `${node.atoms[i].el}${i + 1}-${node.atoms[j].el}${j + 1} ${d(P, i, j).toFixed(2)} A`; }
      }
    }
    say(stretch < 0.1, `   lineup: bonds stay between their lengths at the two ends (outside by at most ${Math.max(0, stretch).toFixed(3)} A ${worstBond})`);
    say(tieStretch < 1.5, `   lineup: loose ions stay near their partners (a tie at most ${Math.max(0, tieStretch).toFixed(2)} A beyond its lengths at the two ends; under 1.5)`);
    say(close > -0.2, `   lineup: molecules pass no closer than start or approach, less 0.2 A (closest ${closeWhere}, ${close.toFixed(2)} A beyond)`);
  }
  if (rx.name.startsWith('SN2')) {
    // the fault of the first round: chloromethane swung about 90 degrees to meet the hydroxide
    const ax = st => sub(st.geo.pos[3], st.geo.pos[2]); // C is atom 2 and Cl atom 3 of '[OH-].[Na+].CCl'
    const turns = steps.slice(1, -1).map(st => angle([ax(steps[0]), [0, 0, 0], ax(st)], 0, 1, 2));
    say(Math.max(...turns) < 30, `   the C-Cl axis turns at most ${Math.max(...turns).toFixed(1)} degrees through the film, gradually (under 30; the first round saw about 90 in a few frames)`);
    const [O, C, Cl] = arrows.flatMap(a => ('pair' in a.from ? [a.from.pair] : [])).concat(arrows[1].from.bond);
    const angs = steps.slice(lineup, -1).map(st => angle(st.geo.pos, O, C, Cl));
    say(Math.min(...angs) > 170, `   O-C-Cl stays on a line through the blend: least ${Math.min(...angs).toFixed(1)} degrees`);
  }
}
// bestFit finds a known turn and shift again
{
  const pts = [[0, 0, 0], [1.5, 0, 0], [0, 1.2, 0.3], [-0.4, 0.2, 1.1], [2, 1, -1]];
  const moved = pts.map(p => add(rotate(p, [0.36, 0.48, 0.8], 1.1), [3, -2, 0.5]));
  const back = bestFit(moved, pts);
  const err = Math.max(...moved.map((p, k) => len(sub(back(p), pts[k]))));
  say(err < 1e-9, `bestFit lays a turned and shifted set back on itself (worst miss ${err.toExponential(1)} A)`);
}
if (bad) { console.log(`${bad} failed`); process.exit(1); }
console.log('film checks pass');
