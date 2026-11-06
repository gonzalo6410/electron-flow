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
// nitrate. No charge filter: structures that differ only in charge separation (amide O=C-N vs
// O(-)-C=N(+)) all contribute with equal weight; the Lewis layer is pure formalism. Ledger item 23.
// Survivors are the contributing structures; they get equal weight; averages are reported.
//
// Group 1 and 2 metals are "ionic": they give their valence electrons to the pool, form no bonds,
// carry charge +valence, and are not octet counted.

import { element } from './elements.js';

export function buildNode({ atoms, adj, charge = 0 }) {
  // atoms: [{ el: 'C', pos: [x,y,z]? }], adj: [[i,j], ...], charge: integer
  const els = atoms.map(a => element(a.el));
  const edges = adj.map(([i, j]) => [Math.min(i, j), Math.max(i, j)]);
  return { atoms: atoms.map((a, k) => ({ el: a.el, pos: a.pos ? [...a.pos] : null, idx: k })), els, edges, charge };
}

export function enumerate(node) {
  const { els, edges, charge } = node;
  const n = els.length;
  const pool = els.reduce((s, e) => s + e.valence, 0) - charge;
  if (pool < 0) throw new Error('more positive charge than valence electrons');
  const sumCeil = els.reduce((s, e) => s + (e.ionic ? 0 : e.ceiling), 0);
  const radicals = pool % 2; // odd electron count means one unpaired electron somewhere

  // bonding edges: both ends non-ionic
  const bondable = edges.map(([i, j]) => !els[i].ionic && !els[j].ionic);
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

  let found = [];
  let S = sBound;
  for (; S >= 0; S--) {
    found = searchOrders(n, els, edges, maxOrder, S);
    if (found.length) break;
  }
  const holes = sumCeil - pool - 2 * S; // in electrons
  const L = pool - 2 * S;                // lone electrons to place

  // place lone electrons for each bond-order assignment
  let structures = [];
  for (const orders of found) {
    const bondsAt = new Array(n).fill(0);
    edges.forEach(([i, j], k) => { bondsAt[i] += orders[k]; bondsAt[j] += orders[k]; });
    const cap = els.map((e, i) => e.ionic ? 0 : e.ceiling - 2 * bondsAt[i]);
    for (const lone of placeLone(cap, L, radicals)) {
      const fc = els.map((e, i) => e.ionic ? e.valence : e.valence - lone[i] - bondsAt[i]);
      const count = els.map((e, i) => e.ionic ? 0 : 2 * bondsAt[i] + lone[i]);
      const cuts = orders.reduce((s, o, k) => s + (bondable[k] && o === 0 ? 1 : 0), 0);
      structures.push({ orders, lone, charge: fc, count, cuts });
    }
  }
  if (!structures.length) throw new Error('no Lewis structure found');

  // second rule: fewest cut pairs (sigma bonds break only where the octet forces it)
  const minCuts = Math.min(...structures.map(s => s.cuts));
  const rejectedByCuts = structures.filter(s => s.cuts > minCuts);
  structures = structures.filter(s => s.cuts === minCuts);

  const w = 1 / structures.length;
  structures.forEach(s => { s.weight = w; });

  return { pool, holes, structures, avg: average(structures, edges.length, n), rejectedByCuts, edges, radicals };
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
    avg: average(structures, edges.length, n), rejectedByCuts: [], radicals: rA.radicals,
    before: rA, after: rB,
  };
}

// All bond-order assignments with total exactly S, respecting per-atom ceilings.
function searchOrders(n, els, edges, maxOrder, S) {
  const m = edges.length;
  const results = [];
  const orders = new Array(m).fill(0);
  const used = new Array(n).fill(0); // bond order already placed at each atom
  // suffix bound: most order the remaining edges could add (ignores atom caps, cheap)
  const suffixMax = new Array(m + 1).fill(0);
  for (let k = m - 1; k >= 0; k--) suffixMax[k] = suffixMax[k + 1] + maxOrder[k];

  function rec(k, sum) {
    if (sum > S) return;
    if (sum + suffixMax[k] < S) return;
    if (k === m) { if (sum === S) results.push([...orders]); return; }
    const [i, j] = edges[k];
    const capI = Math.floor((els[i].ceiling) / 2) - used[i];
    const capJ = Math.floor((els[j].ceiling) / 2) - used[j];
    const top = Math.min(maxOrder[k], capI, capJ, S - sum);
    for (let o = top; o >= 0; o--) {
      orders[k] = o; used[i] += o; used[j] += o;
      rec(k + 1, sum + o);
      used[i] -= o; used[j] -= o;
    }
    orders[k] = 0;
  }
  rec(0, 0);
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
    (result.rejectedByCuts.length ? `, ${result.rejectedByCuts.length} dropped by the cut-pair rule` : ''));
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
