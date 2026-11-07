// node test/film.js   the shapes of the film of every reaction scene (core/moves.js, filmSteps): the
// substrate does not swing about (the first round of the moves review saw chloromethane turn about 90
// degrees at the start), no atom jumps between two steps, and the SN2 film keeps O, C and Cl on a line.
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
  const { node: fnode, steps } = filmSteps(node, out, geo, 21);
  const heavy = biggestPiece(node).filter(i => node.atoms[i].el !== 'H' && !node.fixed[i]);
  // printed, not judged: a molecule may have to change its own shape (butadiene turning s-cis,
  // acetaldehyde's C=O coming into line with the forming enolate)
  const approachTurn = heavy.length >= 2 ? turnOf(heavy, steps[0].geo.pos, steps[2].geo.pos) : 0;
  // the biggest jump of an atom between two steps of the blend (the settle at the end is left out:
  // a loose ion may drift a long way to its new partner there)
  let jump = 0;
  for (let q = 1; q < steps.length - 2; q++) {
    steps[q].geo.pos.forEach((p, i) => { if (!fnode.fixed[i]) jump = Math.max(jump, len(sub(p, steps[q + 1].geo.pos[i]))); });
  }
  console.log(`     ${rx.name}: the heavy atoms of the biggest molecule turn ${approachTurn.toFixed(0)} degrees by the first blended frame`);
  say(jump < 1.2, `   no atom moves more than ${jump.toFixed(2)} A between two of the 21 steps (under 1.2)`);
  // the substrate stays where it was on screen: its heavy atoms that stay together drift little from
  // where they started (the enolate scene's acetaldehyde drifted 1.9 A and turned on screen before
  // each step was laid onto the one before; moves review, 7 Nov 2026)
  const core = filmCore(node, out.after).filter(i => node.atoms[i].el !== 'H');
  let drift = 0;
  for (const st of steps) for (const i of core) drift = Math.max(drift, len(sub(st.geo.pos[i], steps[0].geo.pos[i])));
  say(drift < 1.0, `   the substrate's own heavy atoms drift at most ${drift.toFixed(2)} A from where they were (under 1.0; the benzyne ring opens its angles by itself)`);
  if (rx.name.startsWith('SN2')) {
    // the fault of the first round: chloromethane swung about 90 degrees to meet the hydroxide
    const ax = st => sub(st.geo.pos[3], st.geo.pos[2]); // C is atom 2 and Cl atom 3 of '[OH-].[Na+].CCl'
    const turns = steps.slice(1, -1).map(st => angle([ax(steps[0]), [0, 0, 0], ax(st)], 0, 1, 2));
    say(Math.max(...turns) < 30, `   the C-Cl axis turns at most ${Math.max(...turns).toFixed(1)} degrees through the film, gradually (under 30; the first round saw about 90 in a few frames)`);
    const [O, C, Cl] = arrows.flatMap(a => ('pair' in a.from ? [a.from.pair] : [])).concat(arrows[1].from.bond);
    const angs = steps.slice(1, -1).map(st => angle(st.geo.pos, O, C, Cl));
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
