// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Whose this is: kept with the engine wherever it goes, and shown on the page.
export const ENGINE = Object.freeze({
  name: 'Electron Flow',
  author: 'Gonzalo Fernandez de la Mora',
  notice: 'Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved.',
});
// Lewis engine. Ledger items 7, 8, 13, 23, 24.
//
// A "node" is a set of atoms with an adjacency list (which pairs of atoms may share electrons)
// and an overall charge. The engine enumerates every assignment of bond orders (0..3 per pair)
// and lone electrons (per atom) that
//   - uses exactly the electrons the elements bring minus the overall charge,
//   - never lets any atom count more than its ceiling (2 for H and He, 8 otherwise),
//   - keeps every atom closed shell unless the electron count is odd (then exactly one atom is a radical),
// then keeps only the assignments with the fewest holes (missing electron pairs),
// then among those the ones with the fewest adjacent pairs cut to order 0 (both rules are built into
// the search itself, see searchOrders). His rule (5 Nov 2026):
// sigma bonds do not break inside a node; a pair is cut only when the octet ceiling leaves no room
// for the bond, as the two spare F in SF6. This one rule removes three acetylenes from benzene,
// O=C=O plus a free CH3- from acetate, the free Cl- structures from POCl3, and N+ with one O2- from
// nitrate. Then no formal charge beyond +1 or -1 unless every structure has one. No other charge
// filter: structures that differ only in charge separation (amide O=C-N vs O(-)-C=N(+), CO2's
// O(-)-C#O(+)) all contribute with equal weight; the Lewis layer is pure formalism. Ledger items 8 to 10.
// Survivors are the contributing structures; they get equal weight; averages are reported.
//
// Group 1 and 2 metals are "ionic": they give their valence electrons to the pool, form no bonds,
// carry charge +valence, and are not octet counted.
//
// Pieces. Atoms joined by adjacency (bonds and the loose ties of spectator ions) form a piece. When
// the atoms carry the charges they were written with (a molecule read from SMILES), each piece keeps
// its own electrons: a lone pair cannot jump from an amide in one piece to a carbocation in another,
// and two methyl radicals written as two pieces stay two radicals. Without written charges (the
// known-answer cases, built by hand with one overall charge) every piece shares one pool, as before.

import { element } from './elements.js';

export function buildNode({ atoms, adj, charge = 0 }) {
  // atoms: [{ el: 'C', pos: [x,y,z]?, ion: -1?, written: 0? }], adj: [[i,j], ...], charge: integer
  // ion: -1 marks a spectator anion (a free Cl-, Br-, I-, F-): it keeps a full shell, carries its
  // charge, forms no bond. Group 1 and 2 metals are spectator cations by their element.
  // written: the formal charge the atom was written with (SMILES); only its sum over a piece is used.
  const els = atoms.map(a => element(a.el));
  const edges = adj.map(([i, j]) => [Math.min(i, j), Math.max(i, j)]);
  const fixed = atoms.map((a, i) => {
    const e = els[i];
    if (e.ionic) return { lone: 0, charge: e.valence };
    if (a.ion === -1) return { lone: e.ceiling, charge: e.valence - e.ceiling };
    return null;
  });
  const out = atoms.map((a, k) => ({ el: a.el, pos: a.pos ? [...a.pos] : null, idx: k, ion: a.ion || 0 }));
  if (atoms.every(a => Number.isInteger(a.written))) out.forEach((a, k) => { a.written = atoms[k].written; });
  return { atoms: out, els, edges, charge, fixed };
}

// The pieces of a node: connected parts of the adjacency, each with the electrons it owns.
export function pieces(node) {
  const { atoms, els, edges, charge, fixed } = node;
  const n = atoms.length;
  const parent = atoms.map((_, i) => i);
  const find = i => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  edges.forEach(([i, j]) => { parent[find(i)] = find(j); });
  const byRoot = new Map();
  for (let i = 0; i < n; i++) { const r = find(i); if (!byRoot.has(r)) byRoot.set(r, []); byRoot.get(r).push(i); }
  const written = atoms.every(a => Number.isInteger(a.written));
  let parts = [...byRoot.values()];
  if (!written) parts = [atoms.map((_, i) => i)]; // one pool, the overall charge is all we know
  return parts.map(members => {
    const valence = members.reduce((s, i) => s + (fixed[i] ? 0 : els[i].valence), 0);
    const fixedCharge = members.reduce((s, i) => s + (fixed[i] ? fixed[i].charge : 0), 0);
    const q = written ? members.reduce((s, i) => s + atoms[i].written, 0) : charge;
    return { members, pool: valence - (q - fixedCharge) };
  });
}

export function enumerate(node) {
  const { els, edges, fixed } = node;
  const n = els.length;
  const free = i => !fixed[i];
  // spectator ions keep their electrons and charge; the pool of a piece is what the rest share
  const parts = pieces(node);
  const pool = parts.reduce((s, p) => s + p.pool, 0);
  if (parts.some(p => p.pool < 0)) throw new Error('more positive charge than valence electrons');
  const sumCeil = els.reduce((s, e, i) => s + (free(i) ? e.ceiling : 0), 0);
  const radicals = parts.reduce((s, p) => s + (p.pool % 2), 0); // an odd piece has one unpaired electron
  const pieceOf = new Array(n);
  parts.forEach((p, q) => p.members.forEach(i => { pieceOf[i] = q; }));

  // bonding edges: both ends free
  const bondable = edges.map(([i, j]) => free(i) && free(j));
  const maxOrder = edges.map(([i, j], k) => bondable[k]
    ? Math.min(3, Math.floor(els[i].ceiling / 2), Math.floor(els[j].ceiling / 2))
    : 0);
  const capHalf = els.map((e, i) => (free(i) ? Math.floor(e.ceiling / 2) : 0)); // most bond order at an atom

  // S = total bond order. Holes (in electrons) = sumCeil - pool - 2S, so more bonds means fewer holes.
  // Search from the largest S that could possibly work downward until some assignment exists.
  let sBound = Math.min(
    Math.floor((sumCeil - pool) / 2),
    maxOrder.reduce((a, b) => a + b, 0),
    Math.floor(capHalf.reduce((a, b) => a + b, 0) / 2),
  );
  if (sBound < 0) sBound = 0;

  // for each bond-order assignment, place the lone electrons piece by piece and combine
  const structuresFor = orders => {
    const bondsAt = new Array(n).fill(0);
    edges.forEach(([i, j], k) => { bondsAt[i] += orders[k]; bondsAt[j] += orders[k]; });
    const cap = els.map((e, i) => free(i) ? e.ceiling - 2 * bondsAt[i] : 0);
    const cuts = orders.reduce((s, o, k) => s + (bondable[k] && o === 0 ? 1 : 0), 0);
    // each piece places what its pool leaves after its bonds; a piece with an odd pool has a radical
    const perPiece = parts.map((p, q) => {
      const inBonds = edges.reduce((s, [i], k) => s + (pieceOf[i] === q ? 2 * orders[k] : 0), 0);
      const L = p.pool - inBonds;
      if (L < 0) return [];
      const local = p.members.filter(free);
      return placeLone(local.map(i => cap[i]), L, p.pool % 2).map(v => local.map((i, t) => [i, v[t]]));
    });
    if (perPiece.some(list => !list.length)) return [];
    return combine(perPiece).map(choice => {
      const lone = els.map((e, i) => free(i) ? 0 : fixed[i].lone);
      choice.forEach(pairs => pairs.forEach(([i, v]) => { lone[i] = v; }));
      const fc = els.map((e, i) => free(i) ? e.valence - lone[i] - bondsAt[i] : fixed[i].charge);
      const count = els.map((e, i) => free(i) ? 2 * bondsAt[i] + lone[i] : 0);
      return { orders, lone, charge: fc, count, cuts };
    });
  };

  // Largest S first (fewest holes). At each S the search first allows no cut pair at all, then one,
  // then two, and stops at the first count that gives structures: the second rule (fewest cut pairs)
  // is built into the search, so the many assignments with cuts that it would throw away are never
  // listed (on caffeine they ran to hundreds of thousands and overflowed the stack).
  let structures = [];
  let S = sBound;
  const bondableCount = bondable.filter(Boolean).length;
  outer: for (; S >= 0; S--) {
    for (let budget = 0; budget <= bondableCount; budget++) {
      const found = searchOrders(n, edges, maxOrder, capHalf, S, budget);
      for (const orders of found) for (const st of structuresFor(orders)) structures.push(st);
      if (structures.length) break outer;
    }
  }
  if (!structures.length) throw new Error('no Lewis structure found');
  const holes = sumCeil - pool - 2 * S; // in electrons

  // Ring rule (7 Nov 2026, replacing the rule that refused benzyne). A ring of seven or fewer atoms
  // cannot straighten an atom, so a ring atom keeps one orbital above and below the ring and, if it
  // has only two neighbours, one orbital in the plane of the ring pointing outward. A second pi part
  // at a ring atom has to use that in-plane orbital, and it bonds best with the in-plane orbital of a
  // neighbour doing the same: two neighbouring ring atoms share one bond in the plane. That is
  // benzyne's third bond. So (1) the ring atoms with two pi parts pair off along bonds of order two or
  // more between them, unless no structure manages it (then the bond in the plane goes to a neighbour
  // that has to twist to take it, as in a ring allene such as 1,2-cyclohexadiene: drawn, not precise,
  // item 67). Like the octet target and the charge-size rule, it keeps the structures that pass and
  // never refuses the molecule. Pyridine's C=N(+)=C has no partner for N and is dropped, as before;
  // so is a phenyl anion whose ring would take the lone pair's place in the plane.
  // (2) Electrons in the plane of the ring do not move between structures; only those above and
  // below it do (the two sets of orbitals are at right angles). Each two-neighbour ring atom's
  // in-plane orbital holds the in-plane bond, a lone pair, one electron or nothing, and every
  // structure must agree with the one with the fewest charges. This refuses pyridazine's
  // N(+)#N(+) form, which would empty both nitrogens' in-plane lone pairs into the ring above and
  // below. Resonance locks nothing in (item 11): every contributor fits the same ring.
  const smallRing = smallRingAtoms(n, edges, bondable, 7);
  if (smallRing.size) {
    const sigma = new Array(n).fill(0);
    edges.forEach(([i, j], k) => { if (bondable[k]) { sigma[i]++; sigma[j]++; } });
    const piAt = s => {
      const pi = new Array(n).fill(0);
      edges.forEach(([i, j], k) => { if (s.orders[k] > 1) { pi[i] += s.orders[k] - 1; pi[j] += s.orders[k] - 1; } });
      return pi;
    };
    const paired = s => {
      const pi = piAt(s);
      const two = [...smallRing].filter(i => pi[i] >= 2);
      if (two.some(i => pi[i] > 2)) return false;
      const partners = new Map(two.map(i => [i, []]));
      edges.forEach(([i, j], k) => {
        if (s.orders[k] >= 2 && partners.has(i) && partners.has(j)) { partners.get(i).push(j); partners.get(j).push(i); }
      });
      const used = new Set();
      const match = q => {
        while (q < two.length && used.has(two[q])) q++;
        if (q === two.length) return true;
        const i = two[q];
        used.add(i);
        for (const j of partners.get(i)) {
          if (used.has(j)) continue;
          used.add(j);
          if (match(q + 1)) return true;
          used.delete(j);
        }
        used.delete(i);
        return false;
      };
      return match(0);
    };
    if (structures.some(paired)) structures = structures.filter(paired);
    const inPlane = s => {
      const pi = piAt(s);
      return [...smallRing].filter(i => free(i) && sigma[i] === 2)
        .map(i => (pi[i] >= 2 ? 'b' : s.lone[i] >= 2 ? 'p' : s.lone[i] === 1 ? 'r' : 'e')).join('');
    };
    const charges = s => s.charge.reduce((t, q, i) => t + (free(i) ? Math.abs(q) : 0), 0);
    const reference = structures.reduce((best, s) => (charges(s) < charges(best) ? s : best));
    const want = inPlane(reference);
    structures = structures.filter(s => inPlane(s) === want);
  }

  // third rule (his, 5 Nov 2026): no formal charge beyond +1 or -1 unless every structure has one
  // (SF6 keeps S2+; nitrous oxide loses N2-). Ionic cations are exempt (Mg2+, Ca2+). Judged piece by
  // piece, so an SF6 in one piece does not excuse a 2- somewhere else.
  const rejectedByMagnitude = [];
  for (const p of parts) {
    const within = s => p.members.every(i => !free(i) || Math.abs(s.charge[i]) <= 1);
    if (!structures.some(within)) continue;
    for (const s of structures) if (!within(s)) rejectedByMagnitude.push(s);
    structures = structures.filter(within);
  }

  const w = 1 / structures.length;
  structures.forEach(s => { s.weight = w; });

  return { pool, holes, structures, avg: average(structures, edges.length, n), rejectedByMagnitude, edges, radicals };
}

// every way of taking one item from each list
function combine(lists) {
  let out = [[]];
  for (const list of lists) {
    const next = [];
    for (const partial of out) for (const item of list) next.push([...partial, item]);
    out = next;
  }
  return out;
}

function average(structures, m, n) {
  const sum = (f) => structures.reduce((s, st) => s + st.weight * f(st), 0);
  return {
    orders: Array.from({ length: m }, (_, k) => sum(st => st.orders[k])),
    charge: Array.from({ length: n }, (_, i) => sum(st => st.charge[i])),
    lone: Array.from({ length: n }, (_, i) => sum(st => st.lone[i])),
    count: Array.from({ length: n }, (_, i) => sum(st => st.count[i])),
  };
}

// Transition hybrid between two nodes with the same atoms (ledger items 7, 9, 12).
// Each node is enumerated on its own; the hybrid is the union of their edges, each node's structures
// carried over with order 0 on edges it lacks. Weight: each node gets half, split equally among its
// structures (default, item 37). Enumerating freely over the union adjacency is wrong: for SN2 it admits
// a hydride leaving in place of chloride.
export function hybrid(nodeA, nodeB, mix = 0.5) {
  if (nodeA.atoms.length !== nodeB.atoms.length) throw new Error('hybrid needs the same atoms in both nodes');
  const key = ([i, j]) => `${i}-${j}`;
  const edgeMap = new Map();
  for (const e of [...nodeA.edges, ...nodeB.edges]) if (!edgeMap.has(key(e))) edgeMap.set(key(e), e);
  const edges = [...edgeMap.values()];
  const index = new Map(edges.map((e, k) => [key(e), k]));
  const lift = (node, result, share) => result.structures.map(st => {
    const orders = new Array(edges.length).fill(0);
    node.edges.forEach((e, k) => { orders[index.get(key(e))] = st.orders[k]; });
    return { ...st, orders, weight: share / result.structures.length };
  });
  const rA = enumerate(nodeA), rB = enumerate(nodeB);
  const structures = [...lift(nodeA, rA, 1 - mix), ...lift(nodeB, rB, mix)];
  const n = nodeA.atoms.length;
  const inA = new Set(nodeA.edges.map(key)), inB = new Set(nodeB.edges.map(key));
  const transitionEdges = edges.map(e => !(inA.has(key(e)) && inB.has(key(e)))); // made or broken
  return {
    pool: rA.pool, holes: Math.max(rA.holes, rB.holes), structures, edges, transitionEdges,
    avg: average(structures, edges.length, n), rejectedByMagnitude: [], radicals: rA.radicals,
    before: rA, after: rB,
  };
}

// Atoms lying in a ring of at most `maxSize` atoms (shortest ring through each bondable edge).
function smallRingAtoms(n, edges, bondable, maxSize) {
  const nbrs = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { if (bondable[k]) { nbrs[i].push(j); nbrs[j].push(i); } });
  const inRing = new Set();
  edges.forEach(([i, j], k) => {
    if (!bondable[k]) return;
    // shortest path from j back to i avoiding the edge itself
    const prev = new Map([[j, null]]);
    const queue = [j];
    while (queue.length) {
      const a = queue.shift();
      for (const b of nbrs[a]) {
        if (a === j && b === i) continue;
        if (prev.has(b)) continue;
        prev.set(b, a);
        if (b === i) {
          const ring = [];
          for (let x = i; x !== null; x = prev.get(x)) ring.push(x);
          if (ring.length <= maxSize) ring.forEach(x => inRing.add(x));
          queue.length = 0;
          break;
        }
        queue.push(b);
      }
    }
  });
  return inRing;
}

// All bond-order assignments with total exactly S and at most `budget` pairs cut to order 0,
// respecting per-atom ceilings (capHalf: the most bond order an atom may carry). Pairs are visited
// atom by atom, each atom's pairs to later atoms in one block, so an atom is closed as soon as its
// block is done. Bond order an atom could still have carried but no remaining pair can use is slack
// spent; the total slack is the sum of the ceilings less 2S, and a branch that spends more than that
// cannot reach S, so it is cut. Without this the search on a trityl cation runs for minutes; with
// it, milliseconds.
function searchOrders(n, edges, maxOrder, capHalf, S, budget) {
  const m = edges.length;
  const totalSlack = capHalf.reduce((a, b) => a + b, 0) - 2 * S;
  if (totalSlack < 0) return [];
  // atoms in breadth-first order over the pairs that can bond, then pairs ordered by the earlier atom
  const rank = new Array(n).fill(-1);
  const nbrs = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { if (maxOrder[k] > 0) { nbrs[i].push(j); nbrs[j].push(i); } });
  let next = 0;
  for (let s = 0; s < n; s++) {
    if (rank[s] >= 0) continue;
    rank[s] = next++;
    const queue = [s];
    while (queue.length) { const a = queue.shift(); for (const b of nbrs[a]) if (rank[b] < 0) { rank[b] = next++; queue.push(b); } }
  }
  const order = edges.map((_, k) => k).filter(k => maxOrder[k] > 0)
    .sort((a, b) => {
      const [a1, a2] = [Math.min(rank[edges[a][0]], rank[edges[a][1]]), Math.max(rank[edges[a][0]], rank[edges[a][1]])];
      const [b1, b2] = [Math.min(rank[edges[b][0]], rank[edges[b][1]]), Math.max(rank[edges[b][0]], rank[edges[b][1]])];
      return a1 - b1 || a2 - b2;
    });
  const last = new Array(n).fill(-1); // position in `order` of the last pair touching each atom
  order.forEach((k, q) => { last[edges[k][0]] = q; last[edges[k][1]] = q; });
  const suffixMax = new Array(order.length + 1).fill(0);
  for (let q = order.length - 1; q >= 0; q--) suffixMax[q] = suffixMax[q + 1] + maxOrder[order[q]];

  if (order.length - budget > S) return []; // every pair not cut carries at least one bond
  const results = [];
  const orders = new Array(m).fill(0);
  const used = new Array(n).fill(0); // bond order already placed at each atom
  function rec(q, sum, slack, cuts) {
    if (sum > S || slack > totalSlack) return;
    if (sum + suffixMax[q] < S) return;
    if (sum + (order.length - q) - (budget - cuts) > S) return; // the pairs left, less the cuts left, each need a bond
    if (q === order.length) { if (sum === S) results.push([...orders]); return; }
    const k = order[q];
    const [i, j] = edges[k];
    const top = Math.min(maxOrder[k], capHalf[i] - used[i], capHalf[j] - used[j], S - sum);
    const bottom = cuts < budget ? 0 : 1;
    for (let o = top; o >= bottom; o--) {
      orders[k] = o; used[i] += o; used[j] += o;
      let spent = slack;
      if (last[i] === q) spent += capHalf[i] - used[i];
      if (last[j] === q) spent += capHalf[j] - used[j];
      rec(q + 1, sum + o, spent, cuts + (o === 0 ? 1 : 0));
      used[i] -= o; used[j] -= o;
    }
    orders[k] = 0;
  }
  rec(0, 0, 0, 0);
  return results;
}

// All ways to place L lone electrons given per-atom capacity, even counts except `radicals` odd atoms.
function placeLone(cap, L, radicals) {
  const n = cap.length;
  const out = [];
  const lone = new Array(n).fill(0);
  const suffixCap = new Array(n + 1).fill(0);
  for (let i = n - 1; i >= 0; i--) suffixCap[i] = suffixCap[i + 1] + cap[i];
  function rec(i, left, rad) {
    if (left < 0) return;
    if (left > suffixCap[i]) return;
    if (i === n) { if (left === 0 && rad === 0) out.push([...lone]); return; }
    for (let v = Math.min(cap[i], left); v >= 0; v--) {
      const odd = v % 2;
      if (odd && rad === 0) continue;
      lone[i] = v;
      rec(i + 1, left - v, rad - odd);
    }
    lone[i] = 0;
  }
  rec(0, L, radicals);
  return out;
}

// Pretty printer for tests and the console.
export function describe(node, result, opts = {}) {
  const sym = node.atoms.map(a => a.el);
  const lines = [];
  lines.push(`electrons ${result.pool}, holes ${result.holes} electrons, ${result.structures.length} contributing structure(s)` +
    (result.rejectedByMagnitude.length ? `, ${result.rejectedByMagnitude.length} dropped by the magnitude rule` : ''));
  const fmt = x => (Number.isInteger(x) ? String(x) : x.toFixed(3));
  lines.push('average bond orders: ' + result.edges.map(([i, j], k) => `${sym[i]}${i}-${sym[j]}${j} ${fmt(result.avg.orders[k])}`).join(', '));
  lines.push('average charges: ' + sym.map((s, i) => `${s}${i} ${fmt(result.avg.charge[i])}`).join(', '));
  lines.push('average counts: ' + sym.map((s, i) => `${s}${i} ${fmt(result.avg.count[i])}`).join(', '));
  if (opts.all) {
    result.structures.forEach((st, q) => {
      lines.push(`  structure ${q + 1}: orders [${st.orders.join(' ')}] lone [${st.lone.join(' ')}] charge [${st.charge.join(' ')}]`);
    });
  }
  return lines.join('\n');
}
