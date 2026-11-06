// Lewis engine. Ledger items 7, 8, 13, 23, 24.
//
// A "node" is a set of atoms with an adjacency list (which pairs of atoms may share electrons)
// and an overall charge. The engine enumerates every assignment of bond orders (0..3 per pair)
// and lone electrons (per atom) that
//   - uses exactly the electrons the elements bring minus the overall charge,
//   - never lets any atom count more than its ceiling (2 for H and He, 8 otherwise),
//   - keeps every atom closed shell unless the electron count is odd (then exactly one atom is a radical),
// then keeps only the assignments with the fewest holes (missing electron pairs),
// then among those the ones with the fewest adjacent pairs cut to order 0. His rule (5 Nov 2026):
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

import { element } from './elements.js';

export function buildNode({ atoms, adj, charge = 0 }) {
  // atoms: [{ el: 'C', pos: [x,y,z]?, ion: -1? }], adj: [[i,j], ...], charge: integer
  // ion: -1 marks a spectator anion (a free Cl-, Br-, I-, F-): it keeps a full shell, carries its
  // charge, forms no bond. Group 1 and 2 metals are spectator cations by their element.
  const els = atoms.map(a => element(a.el));
  const edges = adj.map(([i, j]) => [Math.min(i, j), Math.max(i, j)]);
  const fixed = atoms.map((a, i) => {
    const e = els[i];
    if (e.ionic) return { lone: 0, charge: e.valence };
    if (a.ion === -1) return { lone: e.ceiling, charge: e.valence - e.ceiling };
    return null;
  });
  return { atoms: atoms.map((a, k) => ({ el: a.el, pos: a.pos ? [...a.pos] : null, idx: k, ion: a.ion || 0 })), els, edges, charge, fixed };
}

export const isIon = (node, i) => !!node.fixed[i];

export function enumerate(node) {
  const { els, edges, charge, fixed } = node;
  const n = els.length;
  const free = i => !fixed[i];
  // spectator ions keep their electrons and charge; the pool is what the rest share
  const fixedCharge = fixed.reduce((s, f) => s + (f ? f.charge : 0), 0);
  const pool = els.reduce((s, e, i) => s + (free(i) ? e.valence : 0), 0) - (charge - fixedCharge);
  if (pool < 0) throw new Error('more positive charge than valence electrons');
  const sumCeil = els.reduce((s, e, i) => s + (free(i) ? e.ceiling : 0), 0);
  const radicals = pool % 2; // odd electron count means one unpaired electron somewhere

  // bonding edges: both ends free
  const bondable = edges.map(([i, j]) => free(i) && free(j));
  const maxOrder = edges.map(([i, j], k) => bondable[k]
    ? Math.min(3, Math.floor(els[i].ceiling / 2), Math.floor(els[j].ceiling / 2))
    : 0);

  // S = total bond order. Holes (in electrons) = sumCeil - pool - 2S, so more bonds means fewer holes.
  // Search from the largest S that could possibly work downward until some assignment exists.
  let sBound = Math.min(
    Math.floor((sumCeil - pool) / 2),
    maxOrder.reduce((a, b) => a + b, 0),
  );
  if (sBound < 0) sBound = 0;

  // Fewest holes first (largest S), and among those the fewest cut pairs (edges at order 0): the search
  // asks for assignments with total S and at most `cuts` zeros, cuts rising from 0, so cuts are only
  // tried when the octet leaves no other way. Trying zeros everywhere, as before, made fused rings
  // take seconds to minutes.
  const bondableCount = bondable.filter(Boolean).length;
  let found = [];
  let S = sBound;
  outer: for (; S >= 0; S--) {
    for (let cuts = 0; cuts <= bondableCount; cuts++) {
      found = searchOrders(n, els, edges, maxOrder, bondable, S, cuts);
      if (found.length) break outer;
    }
  }
  const holes = sumCeil - pool - 2 * S; // in electrons
  const L = pool - 2 * S;                // lone electrons to place

  // place lone electrons for each bond-order assignment
  let structures = [];
  for (const orders of found) {
    const bondsAt = new Array(n).fill(0);
    edges.forEach(([i, j], k) => { bondsAt[i] += orders[k]; bondsAt[j] += orders[k]; });
    const cap = els.map((e, i) => free(i) ? e.ceiling - 2 * bondsAt[i] : 0);
    for (const lone0 of placeLone(cap, L, radicals)) {
      const lone = lone0.map((v, i) => free(i) ? v : fixed[i].lone);
      const fc = els.map((e, i) => free(i) ? e.valence - lone[i] - bondsAt[i] : fixed[i].charge);
      const count = els.map((e, i) => free(i) ? 2 * bondsAt[i] + lone[i] : 0);
      const cuts = orders.reduce((s, o, k) => s + (bondable[k] && o === 0 ? 1 : 0), 0);
      structures.push({ orders, lone, charge: fc, count, cuts });
    }
  }
  if (!structures.length) throw new Error('no Lewis structure found');

  // second rule, fewest cut pairs, is built into the search above; nothing to drop here
  const rejectedByCuts = [];

  // ring rule: no atom in a ring of seven or fewer atoms carries two pi parts (two double bonds, or a
  // triple bond), because such an atom would be linear and a small ring cannot hold it. This removes
  // the pyridine structures with C=N(+)=C in the ring. Checked against the geometry of the node, not
  // the structure: resonance locks nothing in (item 11), so every contributor must fit the same ring.
  const smallRing = smallRingAtoms(n, edges, bondable, 7);
  if (smallRing.size) {
    const fits = s => {
      const pi = new Array(n).fill(0);
      edges.forEach(([i, j], k) => { if (s.orders[k] > 1) { pi[i] += s.orders[k] - 1; pi[j] += s.orders[k] - 1; } });
      for (const i of smallRing) if (pi[i] >= 2) return false;
      return true;
    };
    structures = structures.filter(fits);
    if (!structures.length) throw new Error('no Lewis structure fits the ring');
  }

  // third rule (his, 5 Nov 2026): no formal charge beyond +1 or -1 unless every structure has one
  // (SF6 keeps S2+; nitrous oxide loses N2-). Ionic cations are exempt (Mg2+, Ca2+).
  const within = s => s.charge.every((q, i) => !free(i) || Math.abs(q) <= 1);
  const rejectedByMagnitude = structures.some(within) ? structures.filter(s => !within(s)) : [];
  if (rejectedByMagnitude.length) structures = structures.filter(within);

  const w = 1 / structures.length;
  structures.forEach(s => { s.weight = w; });

  return { pool, holes, structures, avg: average(structures, edges.length, n), rejectedByCuts, rejectedByMagnitude, edges, radicals };
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
    avg: average(structures, edges.length, n), rejectedByCuts: [], rejectedByMagnitude: [], radicals: rA.radicals,
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

// All bond-order assignments with total exactly S, at most `maxZeros` bondable edges at order 0,
// respecting per-atom ceilings. Edges are visited atom by atom so an atom's limit bites early, and two
// bounds prune: the most the remaining edges could still add, and the least they must add.
function searchOrders(n, els, edges, maxOrder, bondable, S, maxZeros) {
  const m = edges.length;
  // visit order: breadth first from the atom with most edges, each atom's edges together
  const edgesAt = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { edgesAt[i].push(k); edgesAt[j].push(k); });
  const order = [];
  const listed = new Array(m).fill(false);
  const seen = new Array(n).fill(false);
  const starts = [...Array(n).keys()].sort((a, b) => edgesAt[b].length - edgesAt[a].length);
  for (const s0 of starts) {
    if (seen[s0]) continue;
    const queue = [s0]; seen[s0] = true;
    while (queue.length) {
      const i = queue.shift();
      for (const k of edgesAt[i]) {
        if (!listed[k]) { listed[k] = true; order.push(k); }
        const j = edges[k][0] === i ? edges[k][1] : edges[k][0];
        if (!seen[j]) { seen[j] = true; queue.push(j); }
      }
    }
  }
  const capOf = i => Math.floor(els[i].ceiling / 2);
  const used = new Array(n).fill(0);
  const orders = new Array(m).fill(0);
  const results = [];
  const LIMIT = 20000; // more structures than this is not a teaching molecule

  function rec(q, sum, zeros) {
    if (results.length >= LIMIT) return;
    if (sum > S) return;
    // bounds from the edges still to assign
    let most = 0, mustHave = 0;
    for (let t = q; t < order.length; t++) {
      const k = order[t];
      const [i, j] = edges[k];
      const top = Math.min(maxOrder[k], capOf(i) - used[i], capOf(j) - used[j]);
      most += Math.max(0, top);
      if (bondable[k]) mustHave++;
    }
    if (sum + most < S) return;
    if (sum + Math.max(0, mustHave - (maxZeros - zeros)) > S) return;
    if (q === order.length) { if (sum === S) results.push([...orders]); return; }
    const k = order[q];
    const [i, j] = edges[k];
    const top = Math.min(maxOrder[k], capOf(i) - used[i], capOf(j) - used[j], S - sum);
    const low = bondable[k] && zeros >= maxZeros ? 1 : 0;
    for (let o = top; o >= low; o--) {
      orders[k] = o; used[i] += o; used[j] += o;
      rec(q + 1, sum + o, zeros + (bondable[k] && o === 0 ? 1 : 0));
      used[i] -= o; used[j] -= o;
    }
    orders[k] = 0;
  }
  rec(0, 0, 0);
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
    (result.rejectedByCuts.length ? `, ${result.rejectedByCuts.length} dropped by the cut-pair rule` : '') +
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
