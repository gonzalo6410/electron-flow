// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Symmetry of the sigma skeleton, and the grouping of contributing structures into families of
// equivalent structures (the fifteen of SF6 are one family; formamide's two structures are two).
// A relabelling of atoms that keeps every element and every adjacency is a symmetry of the skeleton.
// Two structures are equivalent when some symmetry carries one onto the other (bond orders, lone
// electrons and charges included). Equivalent structures are the same picture with atoms renamed, so
// the viewer cycles through them on its own; families are stepped through by hand (ledger item 13).
//
// Hydrogens that are plain (one single bond, no lone electrons, no charge, in every structure) are
// folded into their heavy atom as a count, which keeps the graph small.
//
// Two structures are compared directly: a search for a relabelling that carries one onto the other,
// with the bond orders, lone electrons and charges of the two pictures narrowing the choices at every
// step. The symmetries of the skeleton are never listed, so a molecule with tens of thousands of them
// (S2F10 has 28800) is grouped exactly; an earlier version listed the first 5000 and split a family
// in such cases.

export function structureFamilies(node, result) {
  const graph = reducedGraph(node, result);
  const sts = result.structures;
  const { keep, label, pairs } = graph;
  // a fingerprint no relabelling can change: sorted atom tokens and sorted pair tokens
  const fingerprint = st => {
    const atomPart = keep.map(i => `${label(i)}:${st.lone[i]}:${st.charge[i]}`).sort().join('|');
    const pairPart = pairs.map(([i, j, k]) => [label(i), label(j)].sort().join('~') + ':' + st.orders[k]).sort().join('|');
    return atomPart + '#' + pairPart;
  };
  const reps = new Map(); // fingerprint -> representatives of the families found so far
  const familyOf = new Array(sts.length);
  sts.forEach((st, q) => {
    const fp = fingerprint(st);
    if (!reps.has(fp)) reps.set(fp, []);
    const list = reps.get(fp);
    const match = list.find(rep => carries(graph, sts[rep], st));
    familyOf[q] = match === undefined ? q : match;
    if (match === undefined) list.push(q);
  });
  const groups = new Map();
  familyOf.forEach((r, q) => { if (!groups.has(r)) groups.set(r, []); groups.get(r).push(q); });
  const families = [...groups.values()].sort((a, b) => a[0] - b[0]);
  return { families };
}

// The skeleton with plain hydrogens folded in and spectator ions left out. Loosely associated ions
// must not break the symmetry of the ion they sit by.
function reducedGraph(node, result) {
  const n = node.atoms.length;
  const sym = node.atoms.map(a => a.el);
  const edges = result.edges;
  const sts = result.structures;
  const edgesAt = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { edgesAt[i].push(k); edgesAt[j].push(k); });
  const plainH = sym.map((el, i) => el === 'H' && edgesAt[i].length === 1 &&
    sts.every(st => st.orders[edgesAt[i][0]] === 1 && st.lone[i] === 0 && st.charge[i] === 0));
  const skip = i => plainH[i] || !!node.fixed[i];
  const keep = [];
  for (let i = 0; i < n; i++) if (!skip(i)) keep.push(i);
  const hCount = new Array(n).fill(0);
  edges.forEach(([i, j]) => { if (plainH[i]) hCount[j]++; if (plainH[j]) hCount[i]++; });
  const label = i => `${sym[i]}/${hCount[i]}`;
  const adj = new Map(keep.map(i => [i, new Set()]));
  const edgeIndex = new Map();
  const pairs = []; // [i, j, k] for every pair of kept atoms
  edges.forEach(([i, j], k) => {
    if (skip(i) || skip(j)) return;
    adj.get(i).add(j); adj.get(j).add(i);
    edgeIndex.set(`${i}-${j}`, k); edgeIndex.set(`${j}-${i}`, k);
    pairs.push([i, j, k]);
  });
  // atoms in breadth-first order from the busiest atom, so each new atom is tied to placed ones
  const order = [];
  const placed = new Set();
  const byDegree = [...keep].sort((a, b) => adj.get(b).size - adj.get(a).size);
  for (const s of byDegree) {
    if (placed.has(s)) continue;
    const queue = [s];
    placed.add(s);
    while (queue.length) {
      const a = queue.shift();
      order.push(a);
      for (const b of adj.get(a)) if (!placed.has(b)) { placed.add(b); queue.push(b); }
    }
  }
  return { keep, label, adj, edgeIndex, pairs, order };
}

// Is there a relabelling of the skeleton that carries structure a onto structure b?
function carries(graph, a, b) {
  const { keep, label, adj, edgeIndex, order } = graph;
  const token = (st, i) => `${label(i)}:${st.lone[i]}:${st.charge[i]}`;
  const perm = new Map(), used = new Set();
  function rec(q) {
    if (q === order.length) return true;
    const i = order[q];
    for (const j of keep) {
      if (used.has(j) || token(a, i) !== token(b, j) || adj.get(j).size !== adj.get(i).size) continue;
      let fine = true;
      for (const [x, y] of perm) {
        const joined = adj.get(i).has(x);
        if (joined !== adj.get(j).has(y)) { fine = false; break; }
        if (joined && a.orders[edgeIndex.get(`${i}-${x}`)] !== b.orders[edgeIndex.get(`${j}-${y}`)]) { fine = false; break; }
      }
      if (!fine) continue;
      perm.set(i, j); used.add(j);
      if (rec(q + 1)) return true;
      used.delete(j); perm.delete(i);
    }
    return false;
  }
  return rec(0);
}

// How many relabellings keep the reduced skeleton (elements, hydrogen counts, adjacency). Listed by
// backtracking up to `limit`; for the tests and the console, the grouping above does not need it.
export function skeletonSymmetries(node, result, limit = 5000) {
  const { keep, label, adj, order } = reducedGraph(node, result);
  const perm = new Map(), used = new Set();
  let count = 0;
  function rec(q) {
    if (count >= limit) return;
    if (q === order.length) { count++; return; }
    const i = order[q];
    for (const j of keep) {
      if (used.has(j) || label(j) !== label(i) || adj.get(j).size !== adj.get(i).size) continue;
      let fine = true;
      for (const [x, y] of perm) { if (adj.get(i).has(x) !== adj.get(j).has(y)) { fine = false; break; } }
      if (!fine) continue;
      perm.set(i, j); used.add(j);
      rec(q + 1);
      used.delete(j); perm.delete(i);
    }
  }
  rec(0);
  return { symmetries: count, truncated: count >= limit };
}
