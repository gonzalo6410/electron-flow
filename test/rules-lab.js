// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Candidate filters tried after the engine's own two rules (fewest holes, fewest cut pairs).
// Prints how many structures survive each candidate, so a rule can be judged on known molecules.
// node test/rules-lab.js
import { buildNode, enumerate } from '../core/lewis.js';
import { element } from '../core/elements.js';

const M = (els, adj, charge) => ({ atoms: els.split(' ').map(el => ({ el })), adj, charge });
const CASES = {
  'carbon dioxide': M('C O O', [[0, 1], [0, 2]], 0),
  'formamide': M('C O N H H H', [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5]], 0),
  'methyl formate (ester)': M('C O O H C H H H', [[0, 1], [0, 2], [0, 3], [2, 4], [4, 5], [4, 6], [4, 7]], 0),
  'methyl vinyl ether': M('C C O C H H H H H H', [[0, 1], [1, 2], [2, 3], [0, 4], [0, 5], [1, 6], [3, 7], [3, 8], [3, 9]], 0),
  'nitrous oxide': M('N N O', [[0, 1], [1, 2]], 0),
  'ozone': M('O O O', [[0, 1], [1, 2]], 0),
  'hydronium': M('O H H H', [[0, 1], [0, 2], [0, 3]], 1),
  'protonated formaldehyde': M('C O H H H', [[0, 1], [0, 2], [0, 3], [1, 4]], 1),
  'diazomethane': M('C N N H H', [[0, 1], [1, 2], [0, 3], [0, 4]], 0),
  'carbon monoxide': M('C O', [[0, 1]], 0),
  'phenol': M('C C C C C C O H H H H H H', [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [6, 7], [1, 8], [2, 9], [3, 10], [4, 11], [5, 12]], 0),
  'aniline': M('C C C C C C N H H H H H H H', [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [6, 7], [6, 8], [1, 9], [2, 10], [3, 11], [4, 12], [5, 13]], 0),
  'nitrate': M('N O O O', [[0, 1], [0, 2], [0, 3]], -1),
  'SF6': M('S F F F F F F', [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6]], 0),
  'acetate': M('C O O C H H H', [[0, 1], [0, 2], [0, 3], [3, 4], [3, 5], [3, 6]], -1),
  'benzene': M('C C C C C C H H H H H H', [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11]], 0),
  'enolate': M('C C O H H H', [[0, 1], [1, 2], [0, 3], [0, 4], [1, 5]], -1),
  'acetone': M('C O C C H H H H H H', [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [2, 6], [3, 7], [3, 8], [3, 9]], 0),
  'pyridine N-oxide': M('N C C C C C O H H H H H', [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11]], 0),
  'sulfate': M('S O O O O', [[0, 1], [0, 2], [0, 3], [0, 4]], -2),
};

// candidate 1: no formal charge of magnitude above 1 unless every structure has one
function magnitude(sts) {
  const ok = sts.filter(s => s.charge.every(q => Math.abs(q) <= 1));
  return ok.length ? ok : sts;
}
// candidate 2: a positive charge may sit only on an atom less electronegative than every atom
// carrying a negative charge, unless no structure satisfies that
function enOrder(sts, els) {
  const ok = sts.filter(s => {
    const plus = s.charge.map((q, i) => q > 0 ? els[i].en : null).filter(x => x !== null);
    const minus = s.charge.map((q, i) => q < 0 ? els[i].en : null).filter(x => x !== null);
    if (!plus.length || !minus.length) return true;
    return Math.max(...plus) < Math.min(...minus);
  });
  return ok.length ? ok : sts;
}
// candidate 3: symmetry. keep only structures unchanged by every relabelling that preserves the
// sigma skeleton and the elements, when at least one such structure exists
function symmetric(sts, node) {
  const n = node.atoms.length;
  const sym = node.atoms.map(a => a.el);
  const nbrs = Array.from({ length: n }, () => new Set());
  node.edges.forEach(([i, j]) => { nbrs[i].add(j); nbrs[j].add(i); });
  // all automorphisms by backtracking (molecules are small)
  const autos = [];
  const perm = new Array(n).fill(-1), used = new Array(n).fill(false);
  function rec(i) {
    if (i === n) { autos.push([...perm]); return; }
    for (let j = 0; j < n; j++) {
      if (used[j] || sym[j] !== sym[i] || nbrs[j].size !== nbrs[i].size) continue;
      // edges to already-mapped atoms must map to edges
      let fine = true;
      for (let k = 0; k < i; k++) { if (nbrs[i].has(k) !== nbrs[j].has(perm[k])) { fine = false; break; } }
      if (!fine) continue;
      perm[i] = j; used[j] = true; rec(i + 1); used[j] = false; perm[i] = -1;
    }
  }
  rec(0);
  const edgeIndex = new Map(node.edges.map(([i, j], k) => [`${Math.min(i, j)}-${Math.max(i, j)}`, k]));
  const invariant = st => autos.every(p => {
    for (let i = 0; i < n; i++) if (st.charge[i] !== st.charge[p[i]] || st.lone[i] !== st.lone[p[i]]) return false;
    return node.edges.every(([i, j], k) => st.orders[k] === st.orders[edgeIndex.get(`${Math.min(p[i], p[j])}-${Math.max(p[i], p[j])}`)]);
  });
  const ok = sts.filter(invariant);
  return ok.length ? ok : sts;
}

const show = (node, sts) => sts.map(s => {
  const sym = node.atoms.map(a => a.el);
  const ch = s.charge.map((q, i) => q ? `${sym[i]}${i}${q > 0 ? '+' : ''}${q}` : null).filter(Boolean).join(' ');
  return ch || 'no charges';
}).join(' | ');

for (const [name, c] of Object.entries(CASES)) {
  const node = buildNode(c);
  const r = enumerate(node);
  const base = r.structures;
  const m = magnitude(base);
  const e = enOrder(m, node.els);
  const s = symmetric(m, node);
  console.log(`\n${name}: engine ${base.length} [${show(node, base)}]`);
  console.log(`  + magnitude rule: ${m.length}`);
  console.log(`  + magnitude + electronegativity order: ${e.length} [${show(node, e)}]`);
  console.log(`  + magnitude + symmetry: ${s.length} [${show(node, s)}]`);
}
