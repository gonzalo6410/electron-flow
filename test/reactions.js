// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/reactions.js   every reaction scene: its textbook arrows make a move, and the products
// are what the textbook says (by formula and charge)
import { buildNode, enumerate } from '../core/lewis.js';
import { nodeFromSmiles } from '../core/smiles.js';
import { commit, transitionFrames } from '../core/moves.js';
import { REACTIONS, resolveArrows, pieceFormulas } from '../core/reactions.js';
import { relax } from '../core/geometry.js';

let bad = 0;
const check = (name, cond) => { console.log(`${cond ? 'ok  ' : 'FAIL'} ${name}`); if (!cond) bad++; };

for (const rx of REACTIONS) {
  const spec = nodeFromSmiles(rx.smiles);
  const node = buildNode(spec);
  const r = enumerate(node);
  const arrows = resolveArrows(spec, rx.arrows);
  // the arrows are drawn on the structure where every tail exists; the main structure first
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  const ordered = [...r.structures].sort((a, b) => spread(a) - spread(b));
  let out = null;
  for (const st of ordered) { out = commit(node, st, arrows, r); if (out.kind === 'move') break; }
  const products = out.kind === 'move' ? pieceFormulas(out.spec, out.after) : [];
  check(`${rx.name}: a move, products ${products.join(' + ')}`, out.kind === 'move' && products.join() === [...rx.products].sort().join());
  if (out.kind !== 'move') { console.log('   ', out.reason); continue; }
  // the movie: every frame has a geometry, each relaxed from the last
  const t0 = Date.now();
  const frames = transitionFrames(node, out.after, 11);
  let pos = relax(node, r).pos;
  const movieNode = { ...node, edges: out.hybrid.edges, fixed: node.fixed.map((f, i) => (node.els[i].ionic ? f : null)) };
  for (const f of frames) pos = relax(movieNode, f.r, { start: pos }).pos;
  const ms = Date.now() - t0;
  check(`   movie: 11 frames relaxed in ${ms} ms`, ms < 1500);
}

console.log(bad ? `${bad} FAILED` : 'all reactions run');
process.exit(bad ? 1 : 0);
