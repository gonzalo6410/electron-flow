// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/mechanisms.js   every mechanism, step by step: each step's textbook arrows make a move;
// the products of one step carry over into the next step's reagents (what carries over includes a
// molecule the step made or changed); whatever leaves or comes in is uncharged as a whole (item 16);
// the atoms that carry over barely move when the next step's scene settles; and the film of every
// step can be played backwards, its last frame being the products and its first the reagents.
import { buildNode, enumerate } from '../core/lewis.js';
import { nodeFromSmiles } from '../core/smiles.js';
import { commit, filmSteps } from '../core/moves.js';
import { resolveArrows } from '../core/reactions.js';
import { MECHANISMS, carryOver, nextShape } from '../core/mechanisms.js';
import { relax } from '../core/geometry.js';
import { element } from '../core/elements.js';

let bad = 0;
const check = (name, cond) => { console.log(`${cond ? 'ok  ' : 'FAIL'} ${name}`); if (!cond) bad++; };
const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);

// the textbook move of a scene, from the given spec and positions
function textbook(rx, spec) {
  const node = buildNode(spec);
  const r = enumerate(node);
  const arrows = resolveArrows(spec, rx.arrows);
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  for (const st of [...r.structures].sort((a, b) => spread(a) - spread(b))) {
    const out = commit(node, st, arrows, r);
    if (out.kind === 'move') return { node, r, out };
  }
  return { node, r, out: null };
}

check('eight mechanisms, every scene in at least one', MECHANISMS.length === 8);
for (const m of MECHANISMS) {
  console.log(`  ${m.name}: ${m.steps.length} step${m.steps.length === 1 ? '' : 's'}`);
  let spec = nodeFromSmiles(m.steps[0].smiles);
  let pos = null;
  for (let s = 0; s < m.steps.length; s++) {
    const rx = m.steps[s];
    const { node, r, out } = textbook(rx, spec);
    check(`   step ${s + 1} (${rx.label}): a move`, !!out);
    if (!out) break;
    const geo = relax(node, r, pos ? { start: pos } : {});
    const film = filmSteps(node, out, geo, 11);
    // played backwards the film starts where the forward one ends and ends where it started
    const back = [...film.steps].reverse();
    const startsOnProducts = back[0].geo.pos.every((p, i) => dist(p, film.endPos[i]) < 1e-9);
    const endsOnReagents = back[back.length - 1].geo.pos.every((p, i) => dist(p, geo.pos[i]) < 1e-9);
    check('      backwards: from the products to the reagents', startsOnProducts && endsOnReagents);
    if (s === m.steps.length - 1) break;
    // into the next step
    const next = nodeFromSmiles(m.steps[s + 1].smiles);
    const shape = nextShape(out.spec, film.endPos, next);
    const { carry } = shape;
    const touched = new Set([...out.made, ...out.broken].flat());
    check(`      carries into step ${s + 2}: ${carry.leaving.length} group(s) leave, ${carry.arriving.length} come in`, carry.map.some(c => c >= 0 && touched.has(c)));
    check('      what leaves and what comes in is uncharged', carry.leavingCharge === 0 && carry.arrivingCharge === 0);
    const kept = new Set(carry.map.filter(c => c >= 0));
    const covered = out.spec.atoms.every((_, i) => kept.has(i) !== carry.leaving.some(g => g.includes(i)));
    check('      every atom on screen either carries over or leaves', covered && kept.size === carry.map.filter(c => c >= 0).length);
    let moved = 0, movedIon = 0;
    carry.map.forEach((c, i) => {
      if (c < 0) return;
      const d = dist(shape.pos[i], film.endPos[c]);
      if (next.atoms[i].ion === -1 || element(next.atoms[i].el).ionic) movedIon = Math.max(movedIon, d); else moved = Math.max(moved, d);
    });
    check(`      the atoms that carry over move ${moved.toFixed(2)} A as the next scene settles (loose ions ${movedIon.toFixed(2)} A)`, moved < 0.5 && movedIon < 2.5);
    spec = next; pos = shape.pos;
  }
}

// the carry-over on its own: a scene carries over into itself whole
{
  const a = nodeFromSmiles('CC(C)(C)[OH2+].[Br-].O');
  const c = carryOver(a, nodeFromSmiles('O.[Br-].CC(C)(C)[OH2+]'));
  check('a scene written in another order carries over whole', c.leaving.length === 0 && c.arriving.length === 0 && c.map.every(x => x >= 0));
}

console.log(bad ? `${bad} FAILED` : 'all mechanisms run');
process.exit(bad ? 1 : 0);
