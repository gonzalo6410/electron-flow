// Rule-based 3D geometry (ledger item 22). Nothing stored per molecule.
//
// 1. Each atom's electron groups = bonded neighbours (average bond order above zero, ionic partners
//    excluded) + lone pairs (average lone electrons / 2, rounded). The group count fixes the shape:
//    2 linear, 3 flat triangle, 4 tetrahedron, 5 trigonal bipyramid, 6 octahedron.
//    In a trigonal bipyramid lone pairs go equatorial and the weakest bonds go axial (so the SN2
//    transition state puts the two half bonds on the axis). In an octahedron two lone pairs go trans.
// 2. Atoms are first placed by walking the bond graph and copying template directions.
// 3. Then positions are relaxed as springs: bond lengths from covalent radii scaled by bond order,
//    angles toward the template angles, atoms sharing a double or partial-double bond held flat
//    together, non-bonded atoms pushed apart. Numerical gradient descent; molecules are small.
// 4. Group 1 and 2 cations are placed beside their adjacency partner at a loose distance.

import { element } from './elements.js';

const DEG = Math.PI / 180;

// template unit vectors by group count; slot roles mark axial/equatorial for 5 and 6
const TEMPLATES = {
  1: [[1, 0, 0]],
  2: [[1, 0, 0], [-1, 0, 0]],
  3: [[1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]],
  4: [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(norm),
  5: [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]],
  6: [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
};
const AXIAL = { 5: [0, 1], 6: [0, 1] };

function norm(v) { const l = Math.hypot(...v); return v.map(x => x / l); }
function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
function len(a) { return Math.hypot(a[0], a[1], a[2]); }

// Bond length rule: sum of covalent radii, times a factor by bond order.
// Order 1 -> 1.00, 2 -> 0.87, 3 -> 0.78 (C-C 1.54, C=C 1.34, C#C 1.20). Below 1 the bond lengthens.
// A partial bond inside one node (the 2/3 bonds of SF6) keeps the single-bond length: the octet
// limit shares the bonds out, it does not stretch them. A partial bond of a transition hybrid
// (the half bonds of the SN2 state) is a bond being made or broken and is drawn long.
export function bondLength(elI, elJ, order, transition = false) {
  const base = element(elI).radius + element(elJ).radius;
  let f;
  if (order >= 2) f = 0.87 - 0.09 * (order - 2);
  else if (order >= 1) f = 1 - 0.13 * (order - 1);
  else f = transition ? 1 + 0.5 * (1 - order) : 1;
  return base * f;
}

// Electron groups per atom: neighbours bonded in the average (any order above zero; the carbon of
// the SN2 state has five) plus whole lone pairs in the average (a fractional lone pair, like the
// oxygen pair that is becoming the O-C bond, is already represented by its partial bond, so it is
// not counted again; XeF4 gets 4 + 2 = 6 and comes out square planar).
export function electronGroups(node, result) {
  const n = node.atoms.length;
  const nbrs = Array.from({ length: n }, () => []); // {j, order}
  result.edges.forEach(([i, j], k) => {
    const o = result.avg.orders[k];
    if (o > 1e-9) { nbrs[i].push({ j, order: o }); nbrs[j].push({ j: i, order: o }); }
  });
  const assoc = Array.from({ length: n }, () => []); // ionic partners
  result.edges.forEach(([i, j], k) => {
    if (result.avg.orders[k] <= 1e-9 && (node.fixed[i] || node.fixed[j])) { assoc[i].push(j); assoc[j].push(i); }
  });
  return node.atoms.map((a, i) => {
    if (node.fixed[i]) return { nbrs: nbrs[i], lonePairs: 0, groups: 1, assoc: assoc[i] };
    const lonePairs = Math.floor(result.avg.lone[i] / 2 + 1e-9);
    const groups = Math.max(1, Math.min(6, nbrs[i].length + lonePairs));
    return { nbrs: nbrs[i], lonePairs: Math.max(0, groups - nbrs[i].length), groups, assoc: assoc[i] };
  });
}

// Shortest ring through the edge i-j, as a set of atoms, or null.
function ringThrough(groups, i, j) {
  const prev = new Map([[j, null]]);
  const queue = [j];
  while (queue.length) {
    const a = queue.shift();
    for (const nb of groups[a].nbrs) {
      const b = nb.j;
      if (a === j && b === i) continue;
      if (prev.has(b)) continue;
      prev.set(b, a);
      if (b === i) {
        const ring = new Set();
        for (let x = i; x !== null; x = prev.get(x)) ring.add(x);
        return ring;
      }
      queue.push(b);
    }
  }
  return null;
}

// Decide which template slot each neighbour and lone pair takes.
function assignSlots(g) {
  const T = TEMPLATES[g.groups];
  const slots = new Array(T.length).fill(null); // 'lp' or neighbour index j
  const axial = AXIAL[g.groups] || [];
  const eq = T.map((_, s) => s).filter(s => !axial.includes(s));
  const order = [...g.nbrs].sort((a, b) => a.order - b.order); // weakest first
  if (g.groups === 5) {
    const lpSlots = eq.slice(0, g.lonePairs);
    lpSlots.forEach(s => { slots[s] = 'lp'; });
    const free = [...axial, ...eq.filter(s => !lpSlots.includes(s))];
    order.forEach((nb, q) => { slots[free[q]] = nb.j; });
  } else if (g.groups === 6) {
    const lpSlots = axial.slice(0, g.lonePairs);
    lpSlots.forEach(s => { slots[s] = 'lp'; });
    const free = T.map((_, s) => s).filter(s => !lpSlots.includes(s));
    order.forEach((nb, q) => { slots[free[q]] = nb.j; });
  } else {
    let s = 0;
    order.forEach(nb => { slots[s++] = nb.j; });
    for (let q = 0; q < g.lonePairs && s < T.length; q++) slots[s++] = 'lp';
  }
  return slots;
}

// rotation taking unit vector a onto unit vector b (Rodrigues)
function rotationAtoB(a, b) {
  const v = cross(a, b), c = dot(a, b);
  if (c < -0.999999) {
    // opposite: rotate 180 degrees about any axis perpendicular to a
    let p = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    p = norm(cross(a, p));
    return u => add(scale(p, 2 * dot(p, u)), scale(u, -1));
  }
  const k = 1 / (1 + c);
  return u => add(add(scale(u, c), scale(cross(v, u), 1)), scale(v, dot(v, u) * k));
}
function rotateAbout(axis, ang) {
  const c = Math.cos(ang), s = Math.sin(ang);
  return u => add(add(scale(u, c), scale(cross(axis, u), s)), scale(axis, dot(axis, u) * (1 - c)));
}

function isTransitionEdge(result, i, j) {
  if (!result.transitionEdges) return false;
  const k = result.edges.findIndex(([a, b]) => (a === i && b === j) || (a === j && b === i));
  return k >= 0 && !!result.transitionEdges[k];
}

// Initial placement by walking the graph.
function initialPlace(node, result, groups) {
  const n = node.atoms.length;
  const pos = new Array(n).fill(null);
  const slotOf = groups.map(assignSlots);
  const placedDir = Array.from({ length: n }, () => ({})); // atom -> {neighbour j: unit dir}
  const sym = node.atoms.map(a => a.el);

  // choose a start per connected piece: most groups
  const visited = new Array(n).fill(false);
  const pieces = [];
  for (let s = 0; s < n; s++) {
    if (visited[s]) continue;
    const piece = [];
    const stack = [s];
    visited[s] = true;
    while (stack.length) {
      const i = stack.pop(); piece.push(i);
      for (const nb of groups[i].nbrs) if (!visited[nb.j]) { visited[nb.j] = true; stack.push(nb.j); }
      for (const j of groups[i].assoc) if (!visited[j]) { visited[j] = true; stack.push(j); }
    }
    pieces.push(piece);
  }

  let offsetX = 0;
  for (const piece of pieces) {
    const start = piece.reduce((best, i) => (groups[i].groups > groups[best].groups ? i : best), piece[0]);
    pos[start] = [offsetX, 0, 0];
    const queue = [{ i: start, parent: null, grand: null, ringSign: 1 }];
    const seen = new Set([start]);
    while (queue.length) {
      const { i, parent, grand, ringSign } = queue.pop(); // depth first, so a ring is walked in one direction
      const T = TEMPLATES[groups[i].groups];
      const slots = slotOf[i];
      // orient template: slot holding parent points to parent
      let rot = u => u;
      if (parent !== null) {
        const ps = slots.indexOf(parent);
        const toParent = norm(sub(pos[parent], pos[i]));
        const r1 = rotationAtoB(ps >= 0 ? T[ps] : T[0], toParent);
        // spin about the parent axis: a child in the same ring as the grandparent goes to the same side
        // (rings close), any other child goes to the opposite side (zigzag chains, staggered), and a
        // child already placed (ring closure) pulls its slot toward where it is.
        let best = r1, bestScore = -Infinity;
        // reference atom: the grandparent, or failing that any other placed neighbour of the parent
        let refAtom = grand;
        if (refAtom === null) {
          const other = groups[parent].nbrs.find(x => x.j !== i && pos[x.j] !== null);
          if (other) refAtom = other.j;
        }
        if (refAtom !== null) {
          const grandRef = refAtom;
          const ref = sub(pos[grandRef], pos[parent]);
          const refPerp = sub(ref, scale(toParent, dot(ref, toParent)));
          if (len(refPerp) > 1e-6) {
            const rp = norm(refPerp);
            const sameRing = {};
            T.forEach((t, s) => {
              const j = slots[s];
              if (j === null || j === 'lp' || j === parent) return;
              const ring = ringThrough(groups, i, j);
              sameRing[j] = !!(ring && ring.has(parent) && ring.has(grandRef));
            });
            // signed dihedral of a direction about the parent axis, 0 = same side as the reference atom
            const dihedral = d => {
              const dp = sub(d, scale(toParent, dot(d, toParent)));
              return Math.atan2(dot(cross(rp, dp), toParent), dot(rp, dp));
            };
            // target dihedral per child: ring neighbour at a tetrahedral atom gauche (chair), ring neighbour
            // at a flat atom cis (flat ring), anything else anti (zigzag chains)
            const target = {};
            T.forEach((t, s) => {
              const j = slots[s];
              if (j === null || j === 'lp' || j === parent) return;
              target[j] = sameRing[j] ? (groups[i].groups === 4 ? ringSign * 60 * DEG : 0) : Math.PI;
            });
            for (let a = 0; a < 360; a += 5) {
              const r2 = rotateAbout(toParent, a * DEG);
              const rr = u => r2(r1(u));
              let score = 0, maxDot = -1;
              T.forEach((t, s) => {
                const j = slots[s];
                if (s === ps || j === 'lp' || j === null) return;
                const d = rr(t);
                const w = sym[j] === 'H' ? 0.3 : 1; // heavy atoms decide, hydrogens follow
                if (pos[j] !== null) score += 3 * dot(d, norm(sub(pos[j], pos[i])));
                else score += w * Math.cos(dihedral(d) - target[j]);
                maxDot = Math.max(maxDot, dot(d, rp));
              });
              score -= 0.5 * maxDot; // among ties, avoid eclipsing (staggered ethane)
              T.forEach((t, s) => { if (slots[s] === 'lp') score += 0.1 * dot(rr(t), rp); });
              if (score > bestScore) { bestScore = score; best = rr; }
            }
          }
        }
        rot = best;
      }
      T.forEach((t, s) => {
        const j = slots[s];
        if (j === null || j === 'lp' || j === parent) return;
        const dir = rot(t);
        placedDir[i][j] = dir;
        if (pos[j] !== null) return; // ring closure, leave where it is
        const nb = groups[i].nbrs.find(x => x.j === j);
        const d = bondLength(sym[i], sym[j], nb.order, isTransitionEdge(result, i, j));
        pos[j] = add(pos[i], scale(dir, d));
        if (!seen.has(j)) { seen.add(j); queue.push({ i: j, parent: i, grand: parent, ringSign: -ringSign }); }
      });
      // ionic partners: beside the atom, away from its neighbours
      for (const j of groups[i].assoc) {
        if (pos[j] !== null) continue;
        let away = [0, 0, 0];
        for (const nb of groups[i].nbrs) away = add(away, norm(sub(pos[i], pos[nb.j])));
        if (len(away) < 1e-6) away = [0, 1, 0];
        const d = element(sym[i]).radius + element(sym[j]).radius + 0.2;
        pos[j] = add(pos[i], scale(norm(away), d));
        if (!seen.has(j)) { seen.add(j); queue.push({ i: j, parent: i, grand: parent, ringSign }); }
      }
    }
    // next piece a few angstrom to the right
    const maxX = Math.max(...piece.map(i => pos[i][0]));
    offsetX = maxX + 4;
  }
  return pos;
}

// Energy of a set of positions under the spring rules.
function makeEnergy(node, result, groups) {
  const n = node.atoms.length;
  const sym = node.atoms.map(a => a.el);
  const bonds = []; // {i, j, d0}
  result.edges.forEach(([i, j], k) => {
    const o = result.avg.orders[k];
    if (o > 1e-9) bonds.push({ i, j, d0: bondLength(sym[i], sym[j], o, !!(result.transitionEdges && result.transitionEdges[k])), order: o });
    else if (node.fixed[i] || node.fixed[j]) bonds.push({ i, j, d0: element(sym[i]).radius + element(sym[j]).radius + 0.2, order: 0, loose: true });
  });
  const bonded = new Set(bonds.map(b => `${b.i}-${b.j}`));
  const isBonded = (i, j) => bonded.has(`${Math.min(i, j)}-${Math.max(i, j)}`);

  // angle terms from template slots
  const angles = []; // {i, a, b, theta0}
  const slotOf = groups.map(assignSlots);
  groups.forEach((g, i) => {
    const T = TEMPLATES[g.groups];
    const slots = slotOf[i];
    for (let s = 0; s < T.length; s++) for (let t = s + 1; t < T.length; t++) {
      const a = slots[s], b = slots[t];
      if (a === null || b === null || a === 'lp' || b === 'lp') continue;
      angles.push({ i, a, b, theta0: Math.acos(Math.max(-1, Math.min(1, dot(T[s], T[t])))) });
    }
  });
  const oneThree = new Set(angles.map(x => `${Math.min(x.a, x.b)}-${Math.max(x.a, x.b)}`));

  // flatness across bonds between two flat (three-group) centres: strong for double or partial-double
  // bonds, weak for a single bond joining two flat centres (conjugation, butadiene's middle bond)
  const flats = []; // {i, j, a, b, k}
  bonds.forEach(({ i, j, order }) => {
    if (groups[i].groups !== 3 || groups[j].groups !== 3) return;
    const k = order > 1.05 ? 1 : 0.3;
    for (const na of groups[i].nbrs) for (const nb of groups[j].nbrs) {
      if (na.j === j || nb.j === i) continue;
      flats.push({ i, j, a: na.j, b: nb.j, k });
    }
  });

  // at a linear atom with two multiple bonds (allene, ketene) the planes of its two neighbours are at
  // right angles, because the two pi parts are
  const twists = []; // {x, z, a, b}
  groups.forEach((g, y) => {
    if (g.groups !== 2 || g.nbrs.length !== 2) return;
    const [x, z] = g.nbrs.map(nb => nb.j);
    const ox = g.nbrs[0].order, oz = g.nbrs[1].order;
    if (ox <= 1.05 || oz <= 1.05) return;
    for (const na of groups[x].nbrs) for (const nb of groups[z].nbrs) {
      if (na.j === y || nb.j === y) continue;
      twists.push({ x, z, a: na.j, b: nb.j });
    }
  });

  const kb = 100, ka = 30, kf = 10, kr = 5;
  function energy(pos) {
    let E = 0;
    for (const b of bonds) {
      const d = len(sub(pos[b.i], pos[b.j]));
      E += (b.loose ? 0.2 : 1) * kb * (d - b.d0) ** 2;
    }
    for (const x of angles) {
      const u = sub(pos[x.a], pos[x.i]), v = sub(pos[x.b], pos[x.i]);
      const c = dot(u, v) / (len(u) * len(v) + 1e-12);
      const th = Math.acos(Math.max(-1, Math.min(1, c)));
      E += ka * (th - x.theta0) ** 2;
    }
    for (const f of flats) {
      // dihedral a-i-j-b; flat means sin = 0
      const b1 = sub(pos[f.i], pos[f.a]), b2 = sub(pos[f.j], pos[f.i]), b3 = sub(pos[f.b], pos[f.j]);
      const n1 = cross(b1, b2), n2 = cross(b2, b3);
      const l1 = len(n1), l2 = len(n2);
      if (l1 < 1e-9 || l2 < 1e-9) continue;
      const s = len(cross(n1, n2)) / (l1 * l2);
      E += kf * f.k * s * s;
    }
    for (const t of twists) {
      // dihedral a-x-z-b about the x...z axis should be 90 degrees: cos = 0
      const b1 = sub(pos[t.x], pos[t.a]), b2 = sub(pos[t.z], pos[t.x]), b3 = sub(pos[t.b], pos[t.z]);
      const n1 = cross(b1, b2), n2 = cross(b2, b3);
      const l1 = len(n1), l2 = len(n2);
      if (l1 < 1e-9 || l2 < 1e-9) continue;
      const c = dot(n1, n2) / (l1 * l2);
      E += kf * c * c;
    }
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (isBonded(i, j) || oneThree.has(`${i}-${j}`)) continue;
      const d = len(sub(pos[i], pos[j]));
      const d0 = 1.6 * (element(sym[i]).radius + element(sym[j]).radius);
      if (d < d0) E += kr * (d0 - d) ** 2;
    }
    return E;
  }
  return energy;
}

export function relax(node, result, opts = {}) {
  const groups = electronGroups(node, result);
  let pos = opts.start ? opts.start.map(p => [...p]) : initialPlace(node, result, groups);
  // a small fixed jitter so a symmetric start (a flat ring of tetrahedral atoms) can fall off its ridge
  if (!opts.start) {
    let seed = 12345;
    const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };
    pos = pos.map(p => [p[0] + 0.05 * rnd(), p[1] + 0.05 * rnd(), p[2] + 0.05 * rnd()]);
  }
  const energy = makeEnergy(node, result, groups);
  const n = pos.length;
  const h = 1e-4;
  let step = 0.01;
  let E = energy(pos);
  const iters = opts.iters ?? 400;
  for (let it = 0; it < iters; it++) {
    const grad = pos.map(() => [0, 0, 0]);
    for (let i = 0; i < n; i++) for (let c = 0; c < 3; c++) {
      pos[i][c] += h; const ep = energy(pos);
      pos[i][c] -= 2 * h; const em = energy(pos);
      pos[i][c] += h;
      grad[i][c] = (ep - em) / (2 * h);
    }
    const gnorm = Math.sqrt(grad.reduce((s, g) => s + dot(g, g), 0));
    if (gnorm < 1e-3) break;
    const trial = pos.map((p, i) => sub(p, scale(grad[i], step)));
    const Et = energy(trial);
    if (Et < E) { pos = trial; E = Et; step *= 1.1; } else { step *= 0.5; }
    if (step < 1e-7) break;
  }
  // centre on the heavy atoms
  const heavy = node.atoms.map((a, i) => i).filter(i => node.atoms[i].el !== 'H');
  const ref = heavy.length ? heavy : node.atoms.map((_, i) => i);
  const c = scale(ref.reduce((s, i) => add(s, pos[i]), [0, 0, 0]), 1 / ref.length);
  pos = pos.map(p => sub(p, c));
  return { pos, energy: E, groups };
}

// Directions for drawing lone pairs: unit vectors that keep away from the bonded neighbours and from
// each other (points repelling on a sphere). Returns an array per atom of `count` unit vectors.
export function lonePairDirections(pos, groups, counts) {
  return pos.map((p, i) => {
    const want = counts[i];
    if (!want) return [];
    const nbrDirs = groups[i].nbrs.map(nb => norm(sub(pos[nb.j], p)));
    // start: opposite the mean neighbour direction, spread by small fixed offsets
    let mean = nbrDirs.reduce((s, d) => add(s, d), [0, 0, 0]);
    let away = len(mean) > 1e-6 ? scale(norm(mean), -1) : [0, 0, 1];
    const perp = norm(cross(away, Math.abs(away[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
    const perp2 = cross(away, perp);
    let dirs = [];
    for (let q = 0; q < want; q++) {
      const ang = (2 * Math.PI * q) / want;
      dirs.push(norm(add(away, add(scale(perp, 0.6 * Math.cos(ang)), scale(perp2, 0.6 * Math.sin(ang))))));
    }
    for (let it = 0; it < 60; it++) {
      dirs = dirs.map((d, q) => {
        let f = [0, 0, 0];
        for (const nd of nbrDirs) { const diff = sub(d, nd); const r = len(diff) + 1e-6; f = add(f, scale(diff, 1 / (r * r * r))); }
        dirs.forEach((o, r) => { if (r !== q) { const diff = sub(d, o); const rr = len(diff) + 1e-6; f = add(f, scale(diff, 1 / (rr * rr * rr))); } });
        return norm(add(d, scale(f, 0.1)));
      });
    }
    return dirs;
  });
}

export function geometryReport(node, result, geo) {
  const sym = node.atoms.map(a => a.el);
  const lines = [];
  lines.push(`energy ${geo.energy.toFixed(4)}`);
  lines.push('groups: ' + geo.groups.map((g, i) => `${sym[i]}${i}:${g.groups}${g.lonePairs ? '(' + g.lonePairs + 'lp)' : ''}`).join(' '));
  const pos = geo.pos;
  const bl = result.edges.map(([i, j], k) => `${sym[i]}${i}-${sym[j]}${j} ${len(sub(pos[i], pos[j])).toFixed(2)}`);
  lines.push('distances: ' + bl.join(', '));
  // a few angles at the first heavy atom with more than one neighbour
  for (let i = 0; i < sym.length; i++) {
    const nb = geo.groups[i].nbrs;
    if (nb.length < 2) continue;
    const angs = [];
    for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++) {
      const u = sub(pos[nb[a].j], pos[i]), v = sub(pos[nb[b].j], pos[i]);
      angs.push(`${sym[nb[a].j]}${nb[a].j}-${sym[i]}${i}-${sym[nb[b].j]}${nb[b].j} ${(Math.acos(dot(u, v) / (len(u) * len(v))) / DEG).toFixed(1)}`);
    }
    lines.push('angles: ' + angs.join(', '));
    break;
  }
  return lines.join('\n');
}
