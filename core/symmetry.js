// Symmetry of the sigma skeleton, and the grouping of contributing structures into families of
// equivalent structures (the fifteen of SF6 are one family; formamide's two structures are two).
// A relabelling of atoms that keeps every element and every adjacency is a symmetry of the skeleton.
// Two structures are equivalent when some symmetry carries one onto the other (bond orders, lone
// electrons and charges included). Equivalent structures are the same picture with atoms renamed, so
// the viewer cycles through them on its own; families are stepped through by hand (ledger item 13).
//
// Hydrogens that are plain (one single bond, no lone electrons, no charge, in every structure) are
// folded into their heavy atom as a count, which keeps the number of symmetries small.

export function structureFamilies(node, result) {
  const n = node.atoms.length;
  const sym = node.atoms.map(a => a.el);
  const edges = result.edges;
  const sts = result.structures;

  // which atoms are plain hydrogens
  const edgesAt = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { edgesAt[i].push(k); edgesAt[j].push(k); });
  const plainH = sym.map((el, i) => el === 'H' && edgesAt[i].length === 1 &&
    sts.every(st => st.orders[edgesAt[i][0]] === 1 && st.lone[i] === 0 && st.charge[i] === 0));

  // reduced graph: non-plain atoms, labelled by element and plain-H count. Loosely associated
  // cations (group 1 and 2) are left out: they must not break the symmetry of the ion they sit by.
  const keep = [];
  for (let i = 0; i < n; i++) if (!plainH[i] && !node.fixed[i]) keep.push(i);
  const hCount = new Array(n).fill(0);
  edges.forEach(([i, j]) => { if (plainH[i]) hCount[j]++; if (plainH[j]) hCount[i]++; });
  const label = i => `${sym[i]}/${hCount[i]}`;
  const adj = new Map(keep.map(i => [i, new Set()]));
  const edgeIndex = new Map();
  const skip = i => plainH[i] || node.fixed[i];
  edges.forEach(([i, j], k) => {
    if (skip(i) || skip(j)) return;
    adj.get(i).add(j); adj.get(j).add(i);
    edgeIndex.set(`${i}-${j}`, k); edgeIndex.set(`${j}-${i}`, k);
  });

  // all symmetries of the reduced graph, by backtracking in degree order
  const order = [...keep].sort((a, b) => adj.get(b).size - adj.get(a).size);
  const perm = new Map(), used = new Set();
  const autos = [];
  const LIMIT = 5000;
  function rec(q) {
    if (autos.length >= LIMIT) return;
    if (q === order.length) { autos.push(new Map(perm)); return; }
    const i = order[q];
    for (const j of keep) {
      if (used.has(j) || label(j) !== label(i) || adj.get(j).size !== adj.get(i).size) continue;
      let fine = true;
      for (const [a, b] of perm) { if (adj.get(i).has(a) !== adj.get(j).has(b)) { fine = false; break; } }
      if (!fine) continue;
      perm.set(i, j); used.add(j);
      rec(q + 1);
      used.delete(j); perm.delete(i);
    }
  }
  rec(0);

  // key of a structure after relabelling by p
  const keyUnder = (st, p) => {
    const atomPart = keep.map(i => `${p.get(i)}:${st.lone[i]}:${st.charge[i]}`).sort().join('|');
    const edgePart = edges.map(([i, j], k) => {
      if (skip(i) || skip(j)) return null;
      const kk = edgeIndex.get(`${p.get(i)}-${p.get(j)}`);
      return `${kk}:${st.orders[k]}`;
    }).filter(Boolean).sort().join('|');
    return atomPart + '#' + edgePart;
  };
  const identity = new Map(keep.map(i => [i, i]));
  const ownKey = sts.map(st => keyUnder(st, identity));
  const keyToIdx = new Map(ownKey.map((k, q) => [k, q]));

  // union-find over images
  const parent = sts.map((_, q) => q);
  const find = q => (parent[q] === q ? q : (parent[q] = find(parent[q])));
  sts.forEach((st, q) => {
    for (const p of autos) {
      const target = keyToIdx.get(keyUnder(st, p));
      if (target !== undefined && target !== q) parent[find(target)] = find(q);
    }
  });
  const groups = new Map();
  sts.forEach((_, q) => { const r = find(q); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(q); });
  const families = [...groups.values()].sort((a, b) => a[0] - b[0]);
  return { families, symmetries: autos.length, truncated: autos.length >= LIMIT };
}
