// node test/moves.js   curved arrows on a structure: resonance, moves to a new node, refusals,
// and the transition hybrid, on the basic mechanisms (ledger items 5, 11, 12, 15)
import { buildNode, enumerate } from '../core/lewis.js';
import { nodeFromSmiles, toSmiles } from '../core/smiles.js';
import { applyArrows, commit, specAfter, combine, withCharges, transitionFrames, describeArrow } from '../core/moves.js';
import { relax } from '../core/geometry.js';

let bad = 0;
const check = (name, cond) => { console.log(`${cond ? 'ok  ' : 'FAIL'} ${name}`); if (!cond) bad++; };
const near = (a, b, eps = 1e-6) => Math.abs(a - b) < eps;
const sym = (spec, i) => spec.atoms[i].el;
const find = (spec, el, nth = 0) => spec.atoms.map((a, i) => (a.el === el ? i : -1)).filter(i => i >= 0)[nth];
const nbrs = (spec, i) => spec.adj.filter(([a, b]) => a === i || b === i).map(([a, b]) => (a === i ? b : a));
const hOn = (spec, i) => nbrs(spec, i).find(j => sym(spec, j) === 'H');
const ekey = ([i, j]) => `${Math.min(i, j)}-${Math.max(i, j)}`;
const hasEdge = (spec, i, j) => spec.adj.some(e => ekey(e) === ekey([i, j]));
const main = r => { const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0); let m = r.structures[0]; for (const st of r.structures) if (spread(st) < spread(m)) m = st; return m; };
const orderOf = (node, r, i, j) => { const k = node.edges.findIndex(e => ekey(e) === ekey([i, j])); return k < 0 ? null : r.avg.orders[k]; };
const run = spec => { const node = buildNode(spec); return { node, r: enumerate(node) }; };

// ---------- SN2: hydroxide and chloromethane ----------
{
  const spec = nodeFromSmiles('[OH-].[Na+].CCl');
  const O = find(spec, 'O'), C = find(spec, 'C'), Cl = find(spec, 'Cl'), Na = find(spec, 'Na');
  const { node, r } = run(spec);
  check('SN2 start: one structure', r.structures.length === 1);
  const arrows = [{ from: { pair: O }, to: { bond: [O, C] } }, { from: { bond: [C, Cl] }, to: { atom: Cl } }];
  console.log('  arrows: ' + arrows.map(a => describeArrow(node, a)).join('; '));
  const out = commit(node, r.structures[0], arrows);
  check('SN2 is a move', out.kind === 'move');
  check('SN2 makes O-C and breaks C-Cl', out.made.length === 1 && ekey(out.made[0]) === ekey([O, C]) && out.broken.length === 1 && ekey(out.broken[0]) === ekey([C, Cl]));
  check('after: chloride is a spectator anion with charge -1', out.spec.atoms[Cl].ion === -1 && out.spec.atoms[Cl].written === -1);
  check('after: oxygen neutral, carbon neutral', out.spec.atoms[O].written === 0 && out.spec.atoms[C].written === 0);
  check('after: sodium follows the charge to the chloride', hasEdge(out.spec, Na, Cl) && !hasEdge(out.spec, Na, O));
  check('after: one structure, no holes', out.rAfter.structures.length === 1 && out.rAfter.holes === 0);
  const h = out.hybrid;
  check('hybrid: two structures', h.structures.length === 2);
  check('hybrid: O-C and C-Cl at one half [STATED]', near(orderOf(h, h, O, C), 0.5) && near(orderOf(h, h, C, Cl), 0.5));
  check('hybrid: -1/2 on O and on Cl, carbon counts eight [STATED]', near(h.avg.charge[O], -0.5) && near(h.avg.charge[Cl], -0.5) && near(h.avg.count[C], 8));
  const frames = transitionFrames(node, out.after, 5);
  check('frames: O-C goes 0, 1/4, 1/2, 3/4, 1', frames.every((f, q) => near(orderOf(f.r, f.r, O, C), q / 4)));
  const geo = relax({ ...node, edges: h.edges }, h);
  const d = (i, j) => Math.hypot(...geo.pos[i].map((x, t) => x - geo.pos[j][t]));
  check('hybrid geometry: O and Cl on opposite sides of carbon (angle above 170)', (() => {
    const u = geo.pos[O].map((x, t) => x - geo.pos[C][t]), v = geo.pos[Cl].map((x, t) => x - geo.pos[C][t]);
    const cos = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (d(O, C) * d(C, Cl));
    return Math.acos(cos) * 180 / Math.PI > 170;
  })());

  // refusals
  const one = commit(node, r.structures[0], [arrows[0]]);
  check('attack without the bond leaving is refused: carbon over its ceiling', one.kind === 'refused' && /C \d+ would hold 10 electrons; 8 is the most/.test(one.reason) && one.atoms.includes(C));
  const far = commit(node, r.structures[0], [{ from: { pair: O }, to: { atom: C } }]);
  check('a lone pair cannot land on another atom', far.kind === 'refused' && /can only become a bond to O/.test(far.reason));
  const metal = commit(node, r.structures[0], [{ from: { pair: Na }, to: { bond: [Na, C] } }]);
  check('sodium has no pair to move', metal.kind === 'refused' && /has given its electrons away/.test(metal.reason));
  const toMetal = commit(node, r.structures[0], [{ from: { pair: O }, to: { bond: [O, Na] } }]);
  check('nothing bonds to sodium', toMetal.kind === 'refused' && /takes no bond/.test(toMetal.reason));
  const twice = commit(node, r.structures[0], [{ from: { bond: [C, Cl] }, to: { atom: Cl } }, { from: { bond: [C, Cl] }, to: { atom: Cl } }]);
  check('a single bond has one pair to move', twice.kind === 'refused' && /no pair left/.test(twice.reason));
  const none = commit(node, r.structures[0], []);
  check('no arrows: nothing', none.kind === 'nothing');
}

// ---------- SN1 in three steps: tert-butyl bromide, water ----------
{
  let spec = nodeFromSmiles('CC(C)(C)Br');
  const Br = find(spec, 'Br');
  const Cq = nbrs(spec, Br)[0];
  let { node, r } = run(spec);
  let out = commit(node, r.structures[0], [{ from: { bond: [Cq, Br] }, to: { atom: Br } }]);
  check('SN1 step 1 is a move that breaks C-Br', out.kind === 'move' && out.broken.length === 1 && !out.made.length);
  check('step 1: bromide is a spectator anion tied to the carbocation carbon', out.spec.atoms[Br].ion === -1 && hasEdge(out.spec, Br, Cq));
  check('step 1: the cation has a hole and keeps its charge on that carbon', out.rAfter.holes === 2 && main(out.rAfter).charge[Cq] === 1);
  check('step 1 hybrid: C-Br at one half, charge -1/2 on Br', near(orderOf(out.hybrid, out.hybrid, Cq, Br), 0.5) && near(out.hybrid.avg.charge[Br], -0.5));

  // bring in water
  const water = nodeFromSmiles('O');
  const c = combine(out.spec, water);
  spec = c.spec;
  const Ow = c.offset + find(water, 'O');
  ({ node, r } = run(spec));
  check('cation and bromide and water: the cation stays a cation (pieces keep their electrons)', r.holes === 2 && main(r).charge[Cq] === 1 && main(r).charge[Br] === -1);
  out = commit(node, r.structures[0], [{ from: { pair: Ow }, to: { bond: [Ow, Cq] } }]);
  check('SN1 step 2: water makes the bond; oxygen now +1, no holes', out.kind === 'move' && out.made.length === 1 && out.spec.atoms[Ow].written === 1 && out.rAfter.holes === 0);
  check('step 2: bromide now beside the oxonium oxygen', hasEdge(out.spec, Br, Ow));

  // a second water takes the proton
  const c2 = combine(out.spec, water);
  spec = c2.spec;
  const Ow2 = c2.offset + find(water, 'O');
  const Hplus = hOn(spec, Ow);
  ({ node, r } = run(spec));
  out = commit(node, r.structures[0], [{ from: { pair: Ow2 }, to: { bond: [Ow2, Hplus] } }, { from: { bond: [Ow, Hplus] }, to: { atom: Ow } }]);
  check('SN1 step 3: proton moves to the second water', out.kind === 'move' && out.made.length === 1 && out.broken.length === 1 && hasEdge(out.spec, Ow2, Hplus) && !hasEdge(out.spec, Ow, Hplus));
  check('step 3: alcohol oxygen neutral, hydronium oxygen +1, bromide beside the hydronium', out.spec.atoms[Ow].written === 0 && out.spec.atoms[Ow2].written === 1 && hasEdge(out.spec, Br, Ow2));
}

// ---------- E2: ethyl bromide and hydroxide ----------
{
  const spec = nodeFromSmiles('[OH-].[Na+].CCBr');
  const O = find(spec, 'O'), Br = find(spec, 'Br'), Na = find(spec, 'Na');
  const Ca = nbrs(spec, Br)[0];
  const Cb = nbrs(spec, Ca).find(j => sym(spec, j) === 'C');
  const Hb = hOn(spec, Cb);
  const { node, r } = run(spec);
  const arrows = [
    { from: { pair: O }, to: { bond: [O, Hb] } },
    { from: { bond: [Cb, Hb] }, to: { bond: [Cb, Ca] } },
    { from: { bond: [Ca, Br] }, to: { atom: Br } },
  ];
  const out = commit(node, r.structures[0], arrows);
  check('E2 is a move: one bond made (O-H), two broken (C-H, C-Br)', out.kind === 'move' && out.made.length === 1 && out.broken.length === 2);
  check('E2 after: C=C double, no holes, one structure', out.rAfter.holes === 0 && out.rAfter.structures.length === 1 && orderOf(out.after, out.rAfter, Ca, Cb) === 2);
  check('E2 after: water neutral, bromide a spectator with sodium beside it', out.spec.atoms[O].written === 0 && out.spec.atoms[Br].ion === -1 && hasEdge(out.spec, Na, Br));
  check('E2 hybrid: C-C at 1.5, O-H at 1/2, C-Br at 1/2', near(orderOf(out.hybrid, out.hybrid, Ca, Cb), 1.5) && near(orderOf(out.hybrid, out.hybrid, O, Hb), 0.5) && near(orderOf(out.hybrid, out.hybrid, Ca, Br), 0.5));
  // without the base: the C-H pair cannot just leave as a proton
  const noBase = commit(node, r.structures[0], arrows.slice(1));
  check('E2 without the base is refused: the hydrogen would be left with nothing', noBase.kind === 'refused' && /left with no electrons at all/.test(noBase.reason) && noBase.atoms.includes(Hb));
}

// ---------- acid and base: acetic acid and hydroxide ----------
{
  const spec = nodeFromSmiles('CC(=O)O.[OH-].[Na+]');
  const Oh = find(spec, 'O', 1); // the acid's OH oxygen (written order: C, C, =O, O, then the hydroxide O)
  const Ob = find(spec, 'O', 2);
  const H = hOn(spec, Oh);
  const Na = find(spec, 'Na');
  const { node, r } = run(spec);
  const out = commit(node, r.structures[0], [{ from: { pair: Ob }, to: { bond: [Ob, H] } }, { from: { bond: [Oh, H] }, to: { atom: Oh } }]);
  check('proton transfer is a move', out.kind === 'move' && hasEdge(out.spec, Ob, H) && !hasEdge(out.spec, Oh, H));
  check('acetate after: two structures, each O at -1/2', out.rAfter.structures.length === 2 && near(out.rAfter.avg.charge[Oh], -0.5));
  check('sodium moves to the acetate', nbrs(out.spec, Na).every(j => sym(out.spec, j) === 'O' && j !== Ob));
}

// ---------- Diels Alder: butadiene and ethylene ----------
{
  const diene = nodeFromSmiles('C=CC=C'), ene = nodeFromSmiles('C=C');
  const c = combine(diene, ene);
  const spec = c.spec;
  const [c1, c2, c3, c4] = [0, 1, 2, 3];
  const [e1, e2] = [c.offset, c.offset + 1];
  const { node, r } = run(spec);
  const st = r.structures.find(s => s.orders[node.edges.findIndex(e => ekey(e) === ekey([c1, c2]))] === 2);
  const out = commit(node, st, [
    { from: { bond: [c1, c2] }, to: { bond: [c1, e1] } },
    { from: { bond: [e1, e2] }, to: { bond: [e2, c4] } },
    { from: { bond: [c3, c4] }, to: { bond: [c2, c3] } },
  ]);
  check('Diels Alder is a move: two bonds made, none broken', out.kind === 'move' && out.made.length === 2 && !out.broken.length);
  check('cyclohexene after: one structure, no holes, no charges, C2=C3', out.rAfter.structures.length === 1 && out.rAfter.holes === 0 && out.spec.atoms.every(a => a.written === 0) && orderOf(out.after, out.rAfter, c2, c3) === 2);
  check('hybrid: the two new bonds at one half, C2-C3 at 1.5', near(orderOf(out.hybrid, out.hybrid, c1, e1), 0.5) && near(orderOf(out.hybrid, out.hybrid, c2, c3), 1.5));
}

// ---------- aldol: enolate of acetaldehyde adds to acetaldehyde ----------
{
  const enolate = nodeFromSmiles('C=C[O-].[Na+]'), aldehyde = nodeFromSmiles('CC=O');
  const c = combine(enolate, aldehyde);
  const spec = c.spec;
  const Ca = 0, Cc = 1, Oe = 2, Na = find(spec, 'Na');
  const Cald = c.offset + 1, Oald = c.offset + 2;
  const { node, r } = run(spec);
  const st = r.structures.find(s => s.orders[node.edges.findIndex(e => ekey(e) === ekey([Ca, Cc]))] === 2);
  check('enolate shows both its structures', r.structures.length === 2 && !!st);
  const out = commit(node, st, [
    { from: { bond: [Ca, Cc] }, to: { bond: [Ca, Cald] } },
    { from: { bond: [Cald, Oald] }, to: { atom: Oald } },
  ]);
  check('aldol addition is a move with one bond made', out.kind === 'move' && out.made.length === 1 && ekey(out.made[0]) === ekey([Ca, Cald]));
  check('alkoxide after: the charge sits on the aldehyde oxygen, the enolate oxygen is a ketone', out.spec.atoms[Oald].written === -1 && out.spec.atoms[Oe].written === 0 && orderOf(out.after, out.rAfter, Cc, Oe) === 2);
  check('sodium moves to the alkoxide oxygen', hasEdge(out.spec, Na, Oald));
  // from the other enolate structure (C-C single, O-): the C=C pair is not there to move
  const other = r.structures.find(s => s !== st);
  const bad2 = commit(node, other, [{ from: { bond: [Ca, Cc] }, to: { bond: [Ca, Cald] } }, { from: { bond: [Cald, Oald] }, to: { atom: Oald } }]);
  check('on the C-C single structure the same arrows move the sigma pair: a move that breaks C-C (the lone pair was the right tail)', bad2.kind === 'move' && bad2.broken.length === 1 && ekey(bad2.broken[0]) === ekey([Ca, Cc]));
}

// ---------- Michael: enolate adds to methyl vinyl ketone ----------
{
  const enolate = nodeFromSmiles('C=C[O-].[Na+]'), mvk = nodeFromSmiles('C=CC(C)=O');
  const c = combine(enolate, mvk);
  const spec = c.spec;
  const Ca = 0, Cc = 1;
  const Cbeta = c.offset, Calpha = c.offset + 1, Ccarb = c.offset + 2, Ok = c.offset + 4;
  const { node, r } = run(spec);
  const st = r.structures.find(s => s.orders[node.edges.findIndex(e => ekey(e) === ekey([Ca, Cc]))] === 2
    && s.orders[node.edges.findIndex(e => ekey(e) === ekey([Cbeta, Calpha]))] === 2
    && s.orders[node.edges.findIndex(e => ekey(e) === ekey([Ccarb, Ok]))] === 2);
  check('MVK and the enolate: the all-neutral-ketone structure exists', !!st);
  const out = commit(node, st, [
    { from: { bond: [Ca, Cc] }, to: { bond: [Ca, Cbeta] } },
    { from: { bond: [Cbeta, Calpha] }, to: { bond: [Calpha, Ccarb] } },
    { from: { bond: [Ccarb, Ok] }, to: { atom: Ok } },
  ]);
  check('Michael addition is a move with one bond made', out.kind === 'move' && out.made.length === 1);
  check('after: a new enolate (two structures, -1/2 on its oxygen and alpha carbon)', out.rAfter.structures.length === 2 && near(out.rAfter.avg.charge[Ok], -0.5) && near(out.rAfter.avg.charge[Calpha], -0.5));
}

// ---------- resonance stays inside the node ----------
{
  const spec = nodeFromSmiles('NC=O');
  const N = find(spec, 'N'), C = find(spec, 'C'), O = find(spec, 'O');
  const { node, r } = run(spec);
  const neutral = r.structures.find(s => s.charge.every(q => q === 0));
  const out = commit(node, neutral, [{ from: { pair: N }, to: { bond: [N, C] } }, { from: { bond: [C, O] }, to: { atom: O } }]);
  check('formamide arrows are resonance', out.kind === 'resonance');
  const other = r.structures.find(s => s !== neutral);
  check('and give the other contributing structure exactly', out.applied.structure.orders.join() === other.orders.join() && out.applied.structure.charge.join() === other.charge.join() && out.applied.structure.lone.join() === other.lone.join());
  const half = commit(node, neutral, [{ from: { pair: N }, to: { bond: [N, C] } }]);
  check('one arrow alone puts carbon over its ceiling', half.kind === 'refused' && /C \d+ would hold 10/.test(half.reason));

  const benzene = nodeFromSmiles('c1ccccc1');
  const b = run(benzene);
  const k = b.r.structures[0];
  const dbl = b.node.edges.map((e, q) => (k.orders[q] === 2 ? e : null)).filter(Boolean);
  const sgl = b.node.edges.map((e, q) => (k.orders[q] === 1 && sym(benzene, e[0]) === 'C' && sym(benzene, e[1]) === 'C' ? e : null)).filter(Boolean);
  const next = ([i, j]) => sgl.find(e => e.includes(j) && !e.includes(i)) || sgl.find(e => e.includes(i) && !e.includes(j));
  const arrows = dbl.map(e => ({ from: { bond: e }, to: { bond: next(e) } }));
  const flip = commit(b.node, k, arrows);
  check('the Kekulé flip is resonance', flip.kind === 'resonance');
  const flipped = flip.applied.structure;
  check('and lands on the other Kekulé structure', b.r.structures.some(s => s !== k && s.orders.join() === flipped.orders.join()));
}

// ---------- SF6: a fluoride leaving is a move; the spare pair taking its turn is resonance ----------
{
  const spec = { atoms: [{ el: 'S' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }, { el: 'F' }], adj: [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6]], charge: 0 };
  const { node, r } = run(spec);
  check('SF6: fifteen structures', r.structures.length === 15);
  const text = toSmiles(node, main(r));
  const back = run(nodeFromSmiles(text));
  check(`SF6 written as ${text} reads back to fifteen structures`, back.r.structures.length === 15);
  const st = r.structures[0];
  const spare = node.edges.map((e, k) => (st.orders[k] === 0 ? e[1] : -1)).filter(i => i >= 0);
  const bonded = node.edges.map((e, k) => (st.orders[k] === 1 ? e[1] : -1)).filter(i => i >= 0);
  const turn = commit(node, st, [{ from: { pair: spare[0] }, to: { bond: [0, spare[0]] } }, { from: { bond: [0, bonded[0]] }, to: { atom: bonded[0] } }]);
  check('a spare fluoride taking its turn is resonance', turn.kind === 'resonance');
  const leave = commit(node, st, [{ from: { bond: [0, bonded[0]] }, to: { atom: bonded[0] } }]);
  check('a fluoride leaving for good is a move to SF5+ with F-', leave.kind === 'move' && leave.broken.length === 1 && leave.spec.atoms[bonded[0]].ion === -1);
  check('SF5+ after: five structures, sulfur still 2+, no holes', leave.rAfter.structures.length === 5 && near(leave.rAfter.avg.charge[0], 2) && leave.rAfter.holes === 0);
}

// ---------- a loose ion attacks: bromide and the tert-butyl cation ----------
{
  const spec = nodeFromSmiles('C[C+](C)C.[Br-]');
  const Br = find(spec, 'Br'), Cq = 1;
  const { node, r } = run(spec);
  check('bromide starts as a loose ion tied to the cation', node.fixed[Br] && hasEdge(spec, Br, Cq));
  const out = commit(node, r.structures[0], [{ from: { pair: Br }, to: { bond: [Br, Cq] } }]);
  check('its pair makes the bond: a move with one bond made', out.kind === 'move' && out.made.length === 1 && ekey(out.made[0]) === ekey([Br, Cq]));
  check('after: tert-butyl bromide, neutral, bromine an ordinary atom again', out.spec.atoms[Br].ion === 0 && out.spec.atoms.every(a => a.written === 0) && out.rAfter.holes === 0);
  check('hybrid: C-Br at one half', near(orderOf(out.hybrid, out.hybrid, Cq, Br), 0.5));
}

// ---------- a picture the molecule does not have: ethylene pulled apart ----------
{
  const spec = nodeFromSmiles('C=C');
  const { node, r } = run(spec);
  const out = commit(node, r.structures[0], [{ from: { bond: [0, 1] }, to: { atom: 1 } }]);
  check('C=C to C(+)-C(-) is refused: holes where the molecule has none', out.kind === 'refused' && /2 more electrons missing/.test(out.reason));
}

console.log(bad ? `${bad} FAILED` : 'all passed');
process.exit(bad ? 1 : 0);
