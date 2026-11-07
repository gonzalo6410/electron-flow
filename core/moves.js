// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Moves: curved arrows drawn on one Lewis structure of a node, and what they do. Ledger items 11,
// 12, 15, 16.
//
// An arrow moves one pair of electrons. Its tail is a lone pair on an atom, or a bond between two
// atoms; its head is an atom (the pair becomes a lone pair there) or two atoms (the pair becomes one
// more bond between them). The usual rules of the formalism: a lone pair can only become a bond to
// its own atom; a bond's pair can only become a lone pair on one of its own atoms, or a bond from
// one of its own atoms to a third atom.
//
//   arrow = { from: { pair: i } | { bond: [i, j] }, to: { atom: k } | { bond: [k, l] } }
//
// applyArrows(node, structure, arrows) applies every arrow at once to the structure they were drawn
// on and returns the structure that results, over the node's pairs plus any pair the arrows made.
// It refuses, with a plain reason and the atom concerned, when a pair is moved that is not there,
// when an arrow breaks the rules above, or when an atom would hold more electrons than its ceiling
// (the octet is a ceiling, item 7): attacking a carbon without a bond leaving is the classic refusal.
//
// If no sigma bond is made or broken (no pair of the node drops to order 0 and no new pair appears)
// the arrows describe resonance inside the node: nothing to commit (item 11). Otherwise they describe
// a move to a new node: the same atoms, the adjacency with the broken pairs removed and the made
// pairs added, and each atom carrying the charge the arrows left it with, so every piece keeps its
// own electrons (core/lewis.js, pieces). An atom left with no bond, a full shell and a charge of -1
// is a spectator anion from then on; a spectator anion that made a bond is an ordinary atom again.
// Loose ions are tied again afterwards, each to the atom of opposite charge with the most charge
// (item 16: no bare ions), so a sodium follows the charge from the hydroxide to the chloride.
//
// commit(spec, node, structure, arrows) does all of that and, for a move, returns the new spec, the
// node after, and the transition hybrid between before and after (core/lewis.js, hybrid).
// transitionFrames(before, after, steps) gives the hybrid at a run of mixes for the mini movie
// (item 15): weights blended from all-before to all-after.

import { element } from './elements.js';
import { buildNode, enumerate, hybrid } from './lewis.js';
import { relax } from './geometry.js';
import { bestFit, rotate, sub, add, scale, dot, cross, len, norm, anyPerp } from './vectors.js';

const key = (i, j) => `${Math.min(i, j)}-${Math.max(i, j)}`;

// a plain name for an atom in a reason: the element and its number in the molecule, counting from 1
export function atomName(node, i) { return `${node.atoms[i].el} ${i + 1}`; }

export function applyArrows(node, structure, arrows, opts = {}) {
  const n = node.atoms.length;
  const { els, fixed } = node;
  const edges = node.edges.map(e => [e[0], e[1]]);
  const index = new Map(edges.map((e, k) => [key(e[0], e[1]), k]));
  const orders = [...structure.orders];
  const lone = [...structure.lone];
  const refuse = (reason, atoms = []) => ({ ok: false, reason, atoms });
  const metal = i => els[i].ionic;
  const onIon = i => fixed[i] && !metal(i); // a spectator anion: it does hold pairs
  const madePairs = [];

  // what each tail asks for, counted against the structure as drawn, so the order of the arrows
  // does not matter
  const pairsOut = new Array(n).fill(0);
  const bondsOut = new Map();
  for (const a of arrows) {
    if (!a || !a.from || !a.to) return refuse('an arrow needs a tail and a head');
    if ('pair' in a.from) {
      const i = a.from.pair;
      if (metal(i)) return refuse(`${atomName(node, i)} has given its electrons away; it has no pair to move`, [i]);
      pairsOut[i] += 2;
      if (pairsOut[i] > lone[i]) return refuse(`${atomName(node, i)} has no lone pair left to move`, [i]);
      if (!('bond' in a.to)) return refuse(`a lone pair on ${atomName(node, i)} can only become a bond to ${atomName(node, i)}`, [i]);
      const [k, l] = a.to.bond;
      if (k === l) return refuse('a bond needs two different atoms');
      if (k !== i && l !== i) return refuse(`a lone pair on ${atomName(node, i)} can only become a bond to ${atomName(node, i)}`, [i, k, l]);
      const other = k === i ? l : k;
      if (metal(other)) return refuse(`${atomName(node, other)} takes no bond; it stays beside the charge`, [other]);
    } else if ('bond' in a.from) {
      const [i, j] = a.from.bond;
      const kk = index.get(key(i, j));
      if (kk === undefined || structure.orders[kk] < 1) return refuse(`there is no bond between ${atomName(node, i)} and ${atomName(node, j)} to move`, [i, j]);
      bondsOut.set(kk, (bondsOut.get(kk) || 0) + 1);
      if (bondsOut.get(kk) > structure.orders[kk]) return refuse(`the bond between ${atomName(node, i)} and ${atomName(node, j)} has no pair left to move`, [i, j]);
      if ('atom' in a.to) {
        const k = a.to.atom;
        if (k !== i && k !== j) return refuse(`the pair of the ${atomName(node, i)}-${atomName(node, j)} bond can only become a lone pair on ${atomName(node, i)} or ${atomName(node, j)}`, [i, j, k]);
      } else if ('bond' in a.to) {
        const [k, l] = a.to.bond;
        if (k === l) return refuse('a bond needs two different atoms');
        const shared = [k, l].filter(x => x === i || x === j);
        if (shared.length === 2) return refuse('that arrow goes back where it came from');
        if (shared.length !== 1) return refuse(`the pair of the ${atomName(node, i)}-${atomName(node, j)} bond can only become a bond from ${atomName(node, i)} or ${atomName(node, j)}`, [i, j, k, l]);
        const other = k === shared[0] ? l : k;
        if (metal(other)) return refuse(`${atomName(node, other)} takes no bond; it stays beside the charge`, [other]);
      } else return refuse('an arrow head must point at an atom or between two atoms');
    } else return refuse('an arrow tail must start at a lone pair or a bond');
  }

  if (opts.shapeOnly) return { ok: true, shapeOnly: true };

  // apply
  for (const a of arrows) {
    if ('pair' in a.from) lone[a.from.pair] -= 2;
    else orders[index.get(key(...a.from.bond))] -= 1;
    if ('atom' in a.to) lone[a.to.atom] += 2;
    else {
      const [k, l] = a.to.bond;
      const id = key(k, l);
      if (!index.has(id)) {
        index.set(id, edges.length);
        edges.push([Math.min(k, l), Math.max(k, l)]);
        orders.push(0);
        madePairs.push(edges[edges.length - 1]);
      }
      orders[index.get(id)] += 1;
    }
  }

  // ceilings
  const bondsAt = new Array(n).fill(0);
  edges.forEach(([i, j], k) => { bondsAt[i] += orders[k]; bondsAt[j] += orders[k]; });
  const count = els.map((e, i) => 2 * bondsAt[i] + lone[i]);
  for (let i = 0; i < n; i++) {
    if (metal(i)) continue;
    if (count[i] > els[i].ceiling) {
      return refuse(`${atomName(node, i)} would hold ${count[i]} electrons; ${els[i].ceiling} is the most it can hold`, [i]);
    }
    // no bare ions (item 16): an atom cannot be left with nothing, a proton must be taken by something
    if (count[i] === 0 && structure.count[i] > 0) {
      return refuse(`${atomName(node, i)} would be left with no electrons at all; something has to take it`, [i]);
    }
  }
  const charge = els.map((e, i) => (metal(i) ? fixed[i].charge : e.valence - lone[i] - bondsAt[i]));

  // what changed in the sigma frame: a pair of the node between two ordinary atoms that drops to
  // order 0 is a bond broken; a new pair, or a tether (one end a spectator ion) that now carries a
  // bond, is a bond made. A pair of the node between ordinary atoms that goes from 0 to 1 (a spare
  // F of SF6 taking its turn) is resonance.
  const broken = [], made = [];
  node.edges.forEach((e, k) => {
    const tether = fixed[e[0]] || fixed[e[1]];
    if (tether) { if (orders[k] >= 1) made.push(e); return; }
    if (structure.orders[k] >= 1 && orders[k] === 0) broken.push(e);
  });
  madePairs.forEach((e, q) => { if (orders[node.edges.length + q] >= 1) made.push(e); });
  const touched = new Set(arrows.flatMap(a => [
    ...('pair' in a.from ? [a.from.pair] : a.from.bond),
    ...('atom' in a.to ? [a.to.atom] : a.to.bond),
  ]));

  const result = { orders, lone, charge, count, cuts: 0, weight: 1 };
  return { ok: true, edges, structure: result, broken, made, touched: [...touched], onIon };
}

// The shape rules alone (tails that exist, heads that an arrow may point at), for checking one arrow
// as it is drawn, before its partners are there: null when fine, else the refusal.
export function arrowShape(node, structure, arrows) {
  const out = applyArrows(node, structure, arrows, { shapeOnly: true });
  return out.ok ? null : out;
}

// Can an arrow start here at all? null when it can, else the refusal.
export function tailShape(node, structure, from) {
  if ('pair' in from) {
    const i = from.pair;
    if (node.els[i].ionic) return { ok: false, reason: `${atomName(node, i)} has given its electrons away; it has no pair to move`, atoms: [i] };
    if (structure.lone[i] < 2) return { ok: false, reason: `${atomName(node, i)} has no lone pair to move`, atoms: [i] };
    return null;
  }
  const [i, j] = from.bond;
  const k = node.edges.findIndex(e => key(e[0], e[1]) === key(i, j));
  if (k < 0 || structure.orders[k] < 1) return { ok: false, reason: `there is no bond between ${atomName(node, i)} and ${atomName(node, j)} to move`, atoms: [i, j] };
  return null;
}

// Is this structure (over the node's pairs and nothing beyond) one of the node's contributing ones?
function isContributor(node, r, applied) {
  if (applied.edges.length !== node.edges.length) return false;
  const s = applied.structure;
  return r.structures.some(t => t.orders.every((o, k) => o === s.orders[k]) && t.lone.every((v, i) => v === s.lone[i]) && t.charge.every((q, i) => q === s.charge[i]));
}

// The spec of the node a successful move leads to: the same atoms, the new adjacency, each atom with
// the charge the arrows left it with; loose ions tied again.
export function specAfter(node, applied) {
  const n = node.atoms.length;
  const { els } = node;
  const metal = i => els[i].ionic;
  const brokenKeys = new Set(applied.broken.map(e => key(e[0], e[1])));
  const adj = [];
  applied.edges.forEach(([i, j], k) => {
    if (applied.structure.orders[k] < 1) {
      // a pair of the node at order 0 in this structure stays in the adjacency when it was at order 0
      // as drawn too (a spare F of SF6: bonded in other structures); a bond the arrows broke leaves;
      // a tether is dropped here and tied again below
      const inNode = k < node.edges.length;
      if (inNode && !node.fixed[i] && !node.fixed[j] && !brokenKeys.has(key(i, j))) adj.push([i, j]);
      return;
    }
    if (metal(i) || metal(j)) return;
    adj.push([i, j]);
  });
  const bonded = new Set(adj.flat());
  const drawn = applied.structure.charge;
  let atoms = node.atoms.map((a, i) => {
    const ion = !metal(i) && !bonded.has(i) && drawn[i] < 0 && applied.structure.count[i] === els[i].ceiling ? -1 : 0;
    return { el: a.el, written: metal(i) ? els[i].valence : drawn[i], ion };
  });
  // the charges the arrows left are those of one structure, which may be a poor one of the new node
  // (an aldol leaves C(+)-O(-) where the node has C=O); each piece keeps its electrons either way,
  // so write the charges of the new node's main structure instead, and tie the ions by those
  const loose = atoms.map((a, i) => metal(i) || a.ion === -1);
  const r0 = enumerate(buildNode({ atoms, adj, charge: atoms.reduce((s, a) => s + a.written, 0) }));
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  let best0 = r0.structures[0];
  for (const st of r0.structures) if (spread(st) < spread(best0)) best0 = st;
  atoms = atoms.map((a, i) => (loose[i] ? a : { ...a, written: Math.round(best0.charge[i]) }));
  // tie each loose ion to a partner: the atom of opposite charge with the most charge, not taken by
  // another ion; failing that, the most (or least) electronegative ordinary atom
  const taken = new Set();
  const nonH = i => atoms[i].el !== 'H';
  for (let i = 0; i < n; i++) {
    if (!loose[i]) continue;
    // a loose ion another loose ion has already tied to keeps that one tie (a sodium chloride is one
    // tie, not two on top of each other)
    if (adj.some(([a, b]) => (a === i && loose[b]) || (b === i && loose[a]))) continue;
    const want = atoms[i].written > 0 ? -1 : 1;
    const score = j => (j === i ? -Infinity : Math.sign(atoms[j].written) === want ? Math.abs(atoms[j].written) : -Infinity);
    let partner = -1, best = -Infinity;
    for (let j = 0; j < n; j++) {
      if (metal(j) && metal(i)) continue; // a metal never ties to a metal
      const s = score(j) - (taken.has(j) ? 0.5 : 0);
      if (s > best) { best = s; partner = j; }
    }
    if (best === -Infinity) {
      partner = -1; let bestEn = want < 0 ? -1 : 99;
      for (let j = 0; j < n; j++) {
        if (j === i || loose[j] || !nonH(j)) continue;
        const en = els[j].en;
        if (want < 0 ? en > bestEn : en < bestEn) { partner = j; bestEn = en; }
      }
    }
    if (partner >= 0) { adj.push([Math.min(i, partner), Math.max(i, partner)]); taken.add(partner); }
  }
  const total = atoms.reduce((s, a) => s + a.written, 0);
  return { atoms, adj, charge: total };
}

// Everything a turn needs: refusal, resonance, or a move with the node after and the hybrid.
// `r` is the node's enumeration (made here when not given).
//   resonance: the arrows give one of the node's own contributing structures (a spare F of SF6
//              taking its turn, the Kekulé flip, the amide's two pictures); nothing to commit.
//   move:      a pair of the node dropped to order 0, or a new pair made; the molecule's sigma frame
//              changed, and the picture is not one the node already had.
//   refused:   a rule broken (reasons from applyArrows), or arrows that make no new pair and drop no
//              pair but still miss the node's structures: a picture with holes where the molecule
//              has none (a C=C pulled apart into C(+) and C(-)). Resonance locks nothing in (item 11),
//              so such a picture cannot be committed either.
export function commit(node, structure, arrows, r = null) {
  if (!arrows.length) return { kind: 'nothing', reason: 'no arrows drawn' };
  const applied = applyArrows(node, structure, arrows);
  if (!applied.ok) return { kind: 'refused', reason: applied.reason, atoms: applied.atoms };
  if (!r) r = enumerate(node);
  if (isContributor(node, r, applied)) {
    return { kind: 'resonance', applied, reason: 'no bond made or broken: this is resonance inside the same molecule' };
  }
  if (!applied.made.length && !applied.broken.length) {
    const holes = node.els.reduce((s, e, i) => s + (node.fixed[i] ? 0 : e.ceiling - applied.structure.count[i]), 0);
    const why = holes > r.holes ? `it leaves ${holes - r.holes} more electrons missing from full shells than the molecule's structures do` : 'it breaks one of the molecule\'s rules';
    return { kind: 'refused', applied, reason: `those arrows give a picture that is not one of this molecule's structures: ${why}`, atoms: applied.touched };
  }
  const spec = specAfter(node, applied);
  const after = buildNode(spec);
  let rAfter;
  try { rAfter = enumerate(after); } catch (e) { return { kind: 'refused', reason: `no structure after that move (${e.message})`, atoms: applied.touched }; }
  const h = hybrid(node, after);
  return { kind: 'move', applied, spec, after, rAfter, hybrid: h, broken: applied.broken, made: applied.made };
}

// A spec whose atoms carry charges: the charges of its main structure (least total charge) when it
// has none written. Needed before two molecules share a scene, so each keeps its own electrons: a
// chloride beside a carbocation must stay a chloride beside a carbocation.
export function withCharges(spec) {
  if (spec.atoms.length && spec.atoms.every(a => Number.isInteger(a.written))) return spec;
  if (!spec.atoms.length) return spec;
  const node = buildNode(spec);
  const r = enumerate(node);
  const spread = st => st.charge.reduce((s, q) => s + Math.abs(q), 0);
  let main = r.structures[0];
  for (const st of r.structures) if (spread(st) < spread(main)) main = st;
  return {
    atoms: spec.atoms.map((a, i) => ({ el: a.el, ion: a.ion || 0, written: Math.round(main.charge[i]) })),
    adj: spec.adj,
    charge: spec.charge,
  };
}

// Two molecules in one scene: the second one's atoms after the first's. Returns the spec and the
// offset of the second one's atoms.
export function combine(specA, specB) {
  const a = withCharges(specA), b = withCharges(specB);
  const offset = a.atoms.length;
  return {
    spec: {
      atoms: [...a.atoms, ...b.atoms],
      adj: [...a.adj, ...b.adj.map(([i, j]) => [i + offset, j + offset])],
      charge: a.charge + b.charge,
    },
    offset,
  };
}

// The hybrid at `steps` mixes from all-before to all-after (both ends included), for the mini movie.
export function transitionFrames(before, after, steps = 21) {
  const frames = [];
  for (let s = 0; s < steps; s++) {
    const mix = steps === 1 ? 1 : s / (steps - 1);
    frames.push({ mix, r: hybrid(before, after, mix) });
  }
  return frames;
}

// Arrows as words, for the status line and for tests: "pair on O 1 to O 1-C 3", "C 3-Cl 4 to Cl 4".
export function describeArrow(node, a) {
  const tail = 'pair' in a.from ? `pair on ${atomName(node, a.from.pair)}` : `${atomName(node, a.from.bond[0])}-${atomName(node, a.from.bond[1])}`;
  const head = 'atom' in a.to ? atomName(node, a.to.atom) : `${atomName(node, a.to.bond[0])}-${atomName(node, a.to.bond[1])}`;
  return `${tail} to ${head}`;
}

// What each pair of the film does, for its shapes and its drawing: `changing` when it is a bond on one
// side of the move only (made or broken) between two atoms that take bonds; `tieBefore` and `tieAfter`
// when it ties a loose ion on that side. A bond can become a tie on the same two atoms: the bromide
// of SN1 step 1 leaves its carbon and is then kept beside it (moves review, 7 Nov 2026: that pair
// showed a long hashed line as the film ended). `first` and `last` are the hybrid at the two ends.
export function filmRoles(before, after, first, last) {
  const has = (nd, [i, j]) => nd.edges.some(e => (e[0] === i && e[1] === j) || (e[0] === j && e[1] === i));
  const tie = (nd, e) => has(nd, e) && !!(nd.fixed[e[0]] || nd.fixed[e[1]]);
  return first.edges.map((e, k) => ({
    changing: !before.els[e[0]].ionic && !before.els[e[1]].ionic && (first.avg.orders[k] > 1e-9) !== (last.avg.orders[k] > 1e-9),
    tieBefore: tie(before, e),
    tieAfter: tie(after, e),
  }));
}

// The shapes of the film of a move (item 15), from the geometry on screen: the reagents first drift
// together (the bonds to be made pull a little, the bonds to be broken let go a little) while the
// biggest molecule is held where it is, so the reagent comes round to the substrate rather than the
// substrate swinging about to meet it; then the weights blend from all-before to all-after, the shape
// relaxed at each of `count` mixes from the one before; last the products settle into their own
// shape (a leaving ion drifts off to its partner), so that loading them afterwards moves nothing.
// For the shape only, a bond being made or broken is held between 0.03 and 0.97 inside the blend: at
// zero its two atoms would no longer be held on one line and the molecule would swing (the first and
// last frames turned chloromethane 70 degrees; moves review, 6 Nov 2026). During the film a leaving
// ion is an ordinary atom whose pairs fill up as its bond goes.
// After every step the whole scene is turned and shifted, all of it together, so that the substrate's
// own atoms (filmCore) sit where they were on screen: that changes no distance or angle, and the
// substrate stays still while the rest moves round it. The soft hold alone let acetaldehyde turn 50
// degrees in the approach of the enolate scene and 120 by its end (moves review, 7 Nov 2026).
// Returns the film's node, the steps ({ mix, r, geo }, the first being the geometry on screen), the
// roles of its pairs (filmRoles) and the positions it ends on.
export function filmSteps(before, out, startGeo, count = 21) {
  const frames = transitionFrames(before, out.after, count);
  const node = { ...before, edges: out.hybrid.edges, fixed: before.fixed.map((f, i) => (before.els[i].ionic ? f : null)) };
  const r0 = frames[0].r;
  const roles = filmRoles(before, out.after, r0, frames[frames.length - 1].r);
  const bondEdge = k => !before.els[r0.edges[k][0]].ionic && !before.els[r0.edges[k][1]].ionic;
  const held = k => (r0.transitionEdges[k] && bondEdge(k)) || roles[k].changing;
  // On an atom at one end of a bond being made, a pair that goes is counted, for the shape only, as
  // already a little gone, as much as that bond is held; and on an atom at one end of a bond being
  // broken, a pair that comes as not quite there yet. Else an oxygen about to attack counts its three
  // pairs and the held bond as five groups at the start and four a frame later, and its hydrogen
  // jumps (moves review, 7 Nov 2026). A bond becoming a tie on the same two atoms is a bond being
  // broken here as well.
  const A = r0.avg.lone, B = frames[frames.length - 1].r.avg.lone;
  const ends = made => new Set(r0.edges.flatMap((e, k) => (held(k) && (r0.avg.orders[k] > 1e-9) !== made ? e : [])));
  const making = ends(true), breaking = ends(false);
  const shapeLone = lone => lone.map((v, i) => {
    if (A[i] - B[i] > 1e-9 && making.has(i)) return Math.min(v, A[i] - 0.06);
    if (B[i] - A[i] > 1e-9 && breaking.has(i)) return Math.min(v, B[i] - 0.06);
    return v;
  });
  const shapeOrders = (fr, f) => ({
    ...fr,
    transitionEdges: fr.transitionEdges.map((t, k) => t || roles[k].changing),
    avg: { ...fr.avg, orders: fr.avg.orders.map((o, k) => (held(k) ? f(o) : o)), lone: shapeLone(fr.avg.lone) },
  });
  // each step is laid onto the one before by the substrate's own atoms (not onto the first: a
  // substrate that changes its own shape, ethyl bromide going flat as ethylene, fits its first shape
  // two ways and would flip between them)
  const core = filmCore(before, out.after);
  const weights = core.map(i => (before.atoms[i].el === 'H' ? 0.25 : 1));
  let pos = startGeo.pos;
  const still = geo => {
    if (core.length < 3) return geo;
    const move = bestFit(core.map(i => geo.pos[i]), core.map(i => pos[i]), weights);
    return { ...geo, pos: geo.pos.map(move) };
  };
  const approachShape = shapeOrders(r0, o => (o > 0.5 ? 0.97 : 0.03));
  const approach = still(relax(node, approachShape, { start: pos, anchor: { atoms: biggestPiece(before), k: 1 } }));
  pos = approach.pos;
  // the lineup, before the film proper: from the picture on screen to the approach (lineupFrames)
  const lineup = lineupFrames(before, startGeo.pos, approach.pos, core, LINEUP_FRAMES).map(p => ({ ...approach, pos: p }));
  const steps = [{ mix: 0, r: r0, geo: startGeo }, ...lineup.map(geo => ({ mix: 0, r: r0, geo }))];
  // each frame is also held softly near the one before (a spring of 1 on every atom), so that a
  // loose part (a proton between two oxygens, a leaving ion) moves a little each frame instead of
  // jumping a whole angstrom in one; the products' own shape comes in the settle at the end
  const all = node.atoms.map((_, i) => i);
  for (const f of frames) {
    const geo = still(relax(node, shapeOrders(f.r, o => Math.min(0.97, Math.max(o, 0.03))), { start: pos, anchor: { atoms: all, k: 1 } }));
    pos = geo.pos;
    steps.push({ mix: f.mix, r: f.r, geo });
  }
  const settled = still(relax(out.after, out.rAfter, { start: pos }));
  steps.push({ mix: 1, r: frames[frames.length - 1].r, geo: settled });
  return { node, steps, roles, endPos: settled.pos, lineup: LINEUP_FRAMES };
}

// The lineup (items 69 and 76): before the film proper the molecules come together and line up for
// the move, each as a whole, turning and sliding without changing its shape. The substrate (the
// molecule holding the film's core) stays where it is; every other molecule swings round the
// substrate's centre, turning with the swing as the moon turns with its orbit (a hydroxide keeps
// facing the carbon it is heading for), turning on itself for whatever the swing leaves, its distance
// from the centre changing evenly. Whatever a molecule changes in its own shape on the way (an O-H
// turning towards the base, butadiene turning s-cis) is done atom by atom along its bonds, each atom
// turning about the one before it, so no bond shortens. A loose ion goes round the atom it is tied
// to, turning with that atom's molecule, by the shorter way round unless the longer one, or a wider
// one, keeps it clear of the other atoms (a bromide must not pass through the water coming in). Where a swing
// brings two molecules closer than they are at either end, it goes wider in the middle, by as little
// as clears them. Moving each atom straight from start to end instead shrinks a molecule that turns
// (at the middle of a half turn its bonds are a third short) and can take one molecule through
// another. Returns LINEUP_FRAMES frames of positions, the last one the approach itself.
const LINEUP_FRAMES = 24;
const ION_CLEAR = 2.2;
// The pieces the lineup moves: molecules by their bonds; a loose ion tied to an atom that is not loose
// goes round that atom (partnerOf); loose ions tied only to each other (a sodium and its chloride)
// move together as one piece.
function lineupPieces(before) {
  const n = before.atoms.length;
  const loose = before.atoms.map((_, i) => !!before.fixed[i]);
  const parent = before.atoms.map((_, i) => i);
  const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  const join = (i, j) => { parent[find(i)] = find(j); };
  const nbr = Array.from({ length: n }, () => []);
  before.edges.forEach(([i, j]) => { if (!loose[i] && !loose[j]) { join(i, j); nbr[i].push(j); nbr[j].push(i); } });
  const partnerOf = new Array(n).fill(-1);
  before.edges.forEach(([i, j]) => {
    if (loose[i] && !loose[j] && partnerOf[i] < 0) partnerOf[i] = j;
    if (loose[j] && !loose[i] && partnerOf[j] < 0) partnerOf[j] = i;
  });
  before.edges.forEach(([i, j]) => { if (loose[i] && loose[j] && partnerOf[i] < 0 && partnerOf[j] < 0) join(i, j); });
  const groups = new Map();
  for (let i = 0; i < n; i++) {
    if (partnerOf[i] >= 0) continue;
    const g = find(i);
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push(i);
  }
  return { n, nbr, partnerOf, find, list: [...groups.values()] };
}

function lineupFrames(before, from, to, core, count) {
  const { n, nbr, partnerOf, find, list } = lineupPieces(before);
  const inCore = new Set(core);
  const coreCount = g => g.filter(i => inCore.has(i)).length;
  const held = core.length ? list.reduce((a, b) => (coreCount(b) > coreCount(a) ? b : a)) : null;
  const centre = (atoms, P) => scale(atoms.reduce((t, i) => add(t, P[i]), [0, 0, 0]), 1 / atoms.length);
  const pivot = held ? centre(held, from) : null;
  const lerp = (a, b, u) => add(a, scale(sub(b, a), u));
  // a turn from one direction and length to another, by the shorter way round or the longer (way -1)
  const arcOf = (d0, d1, way = 1) => {
    const r0 = len(d0), r1 = len(d1);
    if (r0 < 1e-6 || r1 < 1e-6) return null;
    const u0 = scale(d0, 1 / r0), u1 = scale(d1, 1 / r1);
    const phi = Math.acos(Math.max(-1, Math.min(1, dot(u0, u1))));
    const c = cross(u0, u1);
    return { u0, r0, r1, phi: way > 0 ? phi : phi - 2 * Math.PI, axis: len(c) > 1e-6 ? norm(c) : anyPerp(u0) };
  };
  const along = (arc, u, wider = 0) => scale(rotate(arc.u0, arc.axis, arc.phi * u), arc.r0 + (arc.r1 - arc.r0) * u + 4 * u * (1 - u) * wider);
  // the order in which a molecule's atoms are placed, each from the one before it along a bond
  const treeOf = g => {
    const inG = new Set(g), seen = new Set(), order = [];
    for (const root of g) {
      if (seen.has(root)) continue;
      seen.add(root); order.push({ i: root, p: -1 });
      for (let q = order.length - 1; q < order.length; q++) {
        for (const j of nbr[order[q].i]) if (inG.has(j) && !seen.has(j)) { seen.add(j); order.push({ i: j, p: order[q].i }); }
      }
    }
    return order;
  };
  // a molecule's own shape at u of the way from Y0 to Y1 (both in the same frame)
  const shapeAt = (order, Y0, Y1, u) => {
    const y = new Map();
    for (const { i, p } of order) {
      if (p < 0) { y.set(i, lerp(Y0[i], Y1[i], u)); continue; }
      const arc = arcOf(sub(Y0[i], Y0[p]), sub(Y1[i], Y1[p]));
      y.set(i, arc ? add(y.get(p), along(arc, u)) : lerp(Y0[i], Y1[i], u));
    }
    return y;
  };
  const motion = new Map(); // group -> { order, Y1, turnAt(u, v), place(u, wider) }
  for (const g of list) {
    const order = treeOf(g);
    if (g === held) {
      motion.set(g, { g, order, Y1: to, R: null, turnAt: (u, v) => v, centreAt: () => null });
      continue;
    }
    const move = bestFit(g.map(i => from[i]), g.map(i => to[i]));
    const R = matrixOf(move), Rt = transpose(R);
    const c0 = centre(g, from), c1 = centre(g, to);
    const Y1 = [];
    for (const i of g) Y1[i] = add(c0, Rt.map(row => dot(row, sub(to[i], c1))));
    const swing = pivot ? arcOf(sub(c0, pivot), sub(c1, pivot)) : null;
    const turn = axisAngle(R);
    let turnAt, centreAt;
    if (swing && swing.r0 > 0.3 && swing.r1 > 0.3) {
      const self = axisAngle(times(R, transpose(rotation(swing.axis, swing.phi))));
      turnAt = (u, v) => rotate(rotate(v, swing.axis, swing.phi * u), self.axis, self.angle * u);
      centreAt = (u, wider) => add(pivot, along(swing, u, wider));
    } else {
      turnAt = (u, v) => rotate(v, turn.axis, turn.angle * u);
      centreAt = u => lerp(c0, c1, u);
    }
    motion.set(g, { g, order, Y1, R, c0, turnAt, centreAt, swings: !!(swing && swing.r0 > 0.3 && swing.r1 > 0.3) });
  }
  const groupOfAtom = new Map();
  for (const g of list) for (const i of g) groupOfAtom.set(i, g);
  const placeGroups = (u, wider) => {
    const pos = from.map(p => [...p]);
    for (const m of motion.values()) {
      const y = shapeAt(m.order, from, m.Y1, u);
      if (m.g === held) { for (const i of m.g) pos[i] = y.get(i); continue; }
      const c = m.centreAt(u, wider);
      for (const i of m.g) pos[i] = add(c, m.turnAt(u, sub(y.get(i), m.c0)));
    }
    return pos;
  };
  // the loose ions, each about its partner in the frame of the partner's molecule
  const ions = [];
  for (let i = 0; i < n; i++) {
    const j = partnerOf[i];
    if (j < 0) continue;
    const m = motion.get(groupOfAtom.get(j));
    const o1 = m.R ? transpose(m.R).map(row => dot(row, sub(to[i], to[j]))) : sub(to[i], to[j]);
    ions.push({ i, j, m, o0: sub(from[i], from[j]), o1 });
  }
  const placeIons = (pos, u, choice) => {
    ions.forEach((ion, t) => {
      const { way, wider } = choice[t];
      const arc = arcOf(ion.o0, ion.o1, way);
      const o = arc ? along(arc, u, wider) : lerp(ion.o0, ion.o1, u);
      pos[ion.i] = add(pos[ion.j], ion.m.turnAt(u, o));
    });
    return pos;
  };
  // how much closer than at either end two atoms come on a path (in A): `pairs` between molecules
  // for the swings, and each ion against every atom but its partner for the ion's way round
  const gap = (i, j) => Math.min(len(sub(from[i], from[j])), len(sub(to[i], to[j])));
  const pairs = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    if (partnerOf[i] < 0 && partnerOf[j] < 0 && find(i) !== find(j)) pairs.push([i, j, gap(i, j)]);
  }
  const intrusion = (frames, list) => Math.max(0, ...frames.flatMap(P => list.map(([i, j, d]) => d - len(sub(P[i], P[j])))));
  // a loose ion only has to keep clear: 2.2 A from any atom but its partner, or as close as it is at
  // either end if that is closer
  const ionGap = (i, k) => Math.min(ION_CLEAR, gap(i, k));
  const us = Array.from({ length: count - 1 }, (_, k) => (k + 1) / count);
  let wider = 0;
  let base = us.map(u => placeGroups(u, 0));
  if ([...motion.values()].some(m => m.swings)) {
    for (let w = 0.5; w <= 6 && intrusion(base, pairs) > 0.2; w += 0.5) { wider = w; base = us.map(u => placeGroups(u, w)); }
  }
  const choice = ions.map(() => ({ way: 1, wider: 0 }));
  const options = [];
  for (const w of [0, 0.5, 1, 1.5]) for (const way of [1, -1]) options.push({ way, wider: w });
  ions.forEach((ion, t) => {
    const against = [];
    for (let k = 0; k < n; k++) if (k !== ion.i && k !== ion.j && (partnerOf[k] < 0 || ions.findIndex(x => x.i === k) < t)) against.push([ion.i, k, ionGap(ion.i, k)]);
    let best = null, bestScore = Infinity;
    for (const opt of options) {
      choice[t] = opt;
      const frames = us.map((u, q) => placeIons(base[q].map(p => [...p]), u, choice));
      const score = intrusion(frames, against);
      if (score <= 0.05) { best = opt; break; }
      if (score < bestScore) { bestScore = score; best = opt; }
    }
    choice[t] = best;
  });
  const frames = us.map((u, q) => placeIons(base[q], u, choice));
  frames.push(to.map(p => [...p]));
  return frames;
}

// turns as 3 by 3 matrices: the matrix of a rigid move (from bestFit), a turn about an axis, the
// product of two, the reverse of one, and the axis and angle of one
function matrixOf(move) {
  const o = move([0, 0, 0]);
  const col = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(e => sub(move(e), o));
  return [0, 1, 2].map(a => [0, 1, 2].map(b => col[b][a]));
}
function rotation(axis, angle) {
  const col = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(e => rotate(e, axis, angle));
  return [0, 1, 2].map(a => [0, 1, 2].map(b => col[b][a]));
}
function times(A, B) { return A.map(row => [0, 1, 2].map(b => row[0] * B[0][b] + row[1] * B[1][b] + row[2] * B[2][b])); }
function transpose(A) { return [0, 1, 2].map(a => [0, 1, 2].map(b => A[b][a])); }
function axisAngle(R) {
  const angle = Math.acos(Math.max(-1, Math.min(1, (R[0][0] + R[1][1] + R[2][2] - 1) / 2)));
  if (angle < 1e-6) return { axis: [1, 0, 0], angle: 0 };
  const axis = [R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]];
  if (len(axis) > 1e-6) return { axis: norm(axis), angle };
  // a half turn: the axis from the diagonal
  const k = [0, 1, 2].reduce((m, i) => (R[i][i] > R[m][m] ? i : m), 0);
  const nk = Math.sqrt(Math.max(1e-12, (R[k][k] + 1) / 2));
  return { axis: norm([0, 1, 2].map(i => (i === k ? nk : (R[i][k] + R[k][i]) / (4 * nk)))), angle };
}

// the molecules of a node, each as a list of its atoms, by its bonds (a loose ion is a molecule of
// its own, not counted with its partner)
function pieces(node) {
  const parent = node.atoms.map((_, i) => i);
  const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  node.edges.forEach(([i, j]) => { if (!node.fixed[i] && !node.fixed[j]) parent[find(i)] = find(j); });
  const groups = new Map();
  node.atoms.forEach((_, i) => { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); });
  return [...groups.values()];
}

// the atoms of the molecule with the most atoms; none when there is only one molecule
export function biggestPiece(node) {
  const all = pieces(node);
  if (all.length < 2) return [];
  return all.reduce((a, b) => (b.length > a.length ? b : a));
}

// The atoms that hold the picture still during a film: those of the biggest molecule before the move
// (the only one, if there is one) that stay together after it, so a leaving group or a proton that
// goes is left out: CH3 of chloromethane, the tert-butyl of SN1 step 1.
export function filmCore(before, after) {
  const big = pieces(before).reduce((a, b) => (b.length > a.length ? b : a), []);
  const where = new Map();
  pieces(after).forEach((p, q) => p.forEach(i => where.set(i, q)));
  const count = new Map();
  for (const i of big) count.set(where.get(i), (count.get(where.get(i)) || 0) + 1);
  let keep = -1, most = 0;
  for (const [q, c] of count) if (c > most) { most = c; keep = q; }
  return big.filter(i => where.get(i) === keep);
}
