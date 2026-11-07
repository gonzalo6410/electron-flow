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
