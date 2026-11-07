// Rule-based 3D geometry (ledger item 22). Nothing stored per molecule.
//
// 1. Each atom's electron groups = bonded neighbours (average bond order above zero, ionic partners
//    excluded) + lone pairs (average lone electrons / 2, rounded). The group count fixes the shape:
//    2 linear, 3 flat triangle, 4 tetrahedron, 5 trigonal bipyramid, 6 octahedron.
//    In a trigonal bipyramid lone pairs go equatorial and the weakest bonds go axial (so the SN2
//    transition state puts the two half bonds on the axis). In an octahedron two lone pairs go trans.
// 2. Atoms are first placed by walking the bond graph and copying template directions (a flat atom
//    in a small ring uses the ring's polygon angle, see atomTemplate).
// 3. Then positions are relaxed as springs: bond lengths from covalent radii scaled by bond order,
//    angles toward the template angles as adjusted by the ring and lone pair rules (see
//    angleTargets), atoms sharing a double or partial-double bond held flat together, non-bonded
//    atoms pushed apart. The forces are exact derivatives of the spring energy and the relaxation
//    is a quasi-Newton descent (L-BFGS), so twenty atoms take a few milliseconds.
// 4. Group 1 and 2 cations are placed beside their adjacency partner at a loose distance, held by
//    a tether too weak to bend the partner.

import { element } from './elements.js';
import { sub, add, scale, dot, cross, len, norm, rotate } from './vectors.js';

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

// Bond length rule: sum of covalent radii, times a factor by bond order.
// Order 1 -> 1.00, 2 -> 0.87, 3 -> 0.78 (C-C 1.54, C=C 1.34, C#C 1.20). Below 1 the bond lengthens.
// A partial bond inside one node (the 2/3 bonds of SF6) keeps the single-bond length: the octet
// limit shares the bonds out, it does not stretch them. A partial bond of a transition hybrid
// (the half bonds of the SN2 state) is a bond being made or broken and is drawn long: 1.25 times
// the single length at order 1/2, and below that longer still, up to twice the single length as
// the order goes to zero, so that a move played as a film shows the two atoms drifting together
// before the bond takes hold (the spring of such a bond is also weak at low order, springTerms).
export function bondLength(elI, elJ, order, transition = false) {
  const base = element(elI).radius + element(elJ).radius;
  let f;
  if (order >= 2) f = 0.87 - 0.09 * (order - 2);
  else if (order >= 1) f = 1 - 0.13 * (order - 1);
  else if (!transition) f = 1;
  else if (order >= 0.5) f = 1 + 0.5 * (1 - order);
  else f = 1.25 + 1.5 * (0.5 - order);
  return base * f;
}

// Length of the bond on edge k of a result, by the rule above plus one more rule: a single bond
// from a positively charged atom of the third row or beyond to a negatively charged partner (the
// S2+ to O- bonds of sulfate, P+ to O- of phosphate, Cl3+ to O- of perchlorate, S+ to O- of DMSO)
// is drawn at the double-bond length. The charge separation that the octet ceiling forces on these
// atoms is a short, strongly polar bond: sulfate's S-O is 1.49 angstrom, not the 1.71 of a plain
// single bond. Charges are the averages over the structures.
export function edgeLength(node, result, k) {
  const [i, j] = result.edges[k];
  const sym = node.atoms.map(a => a.el);
  const order = result.avg.orders[k];
  const transition = !!(result.transitionEdges && result.transitionEdges[k]);
  if (Math.abs(order - 1) < 1e-6 && !transition) {
    const qi = result.avg.charge[i], qj = result.avg.charge[j];
    const pos = qi >= 0.5 - 1e-9 ? i : qj >= 0.5 - 1e-9 ? j : -1;
    const neg = qi <= -0.5 + 1e-9 ? i : qj <= -0.5 + 1e-9 ? j : -1;
    if (pos >= 0 && neg >= 0 && element(sym[pos]).Z > 10) return bondLength(sym[i], sym[j], 2);
  }
  return bondLength(sym[i], sym[j], order, transition);
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
  // pairs: the average lone pairs, fractions included, used by the angle rules
  return node.atoms.map((a, i) => {
    if (node.fixed[i]) return { nbrs: nbrs[i], lonePairs: 0, pairs: 0, groups: 1, assoc: assoc[i] };
    const lonePairs = Math.floor(result.avg.lone[i] / 2 + 1e-9);
    const groups = Math.max(1, Math.min(6, nbrs[i].length + lonePairs));
    return { nbrs: nbrs[i], lonePairs: Math.max(0, groups - nbrs[i].length), pairs: result.avg.lone[i] / 2, groups, assoc: assoc[i] };
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

// Shortest ring holding the path a-i-b (the two ring bonds of atom i), as a set of atoms, or null.
function ringThroughAngle(groups, a, i, b) {
  const prev = new Map([[a, null]]);
  const queue = [a];
  while (queue.length) {
    const x = queue.shift();
    for (const nb of groups[x].nbrs) {
      const y = nb.j;
      if (y === i || prev.has(y)) continue;
      prev.set(y, x);
      if (y === b) {
        const ring = new Set([i]);
        for (let z = b; z !== null; z = prev.get(z)) ring.add(z);
        return ring;
      }
      queue.push(y);
    }
  }
  return null;
}

// A two-group atom whose two neighbours close a ring of at most seven atoms through it.
function inSmallRing(groups, i) {
  const g = groups[i];
  if (g.groups !== 2 || g.nbrs.length !== 2) return false;
  const ring = ringThroughAngle(groups, g.nbrs[0].j, i, g.nbrs[1].j);
  return !!ring && ring.size <= 7;
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
function rotateAbout(axis, ang) { return u => rotate(u, axis, ang); }

function edgeIndex(result, i, j) {
  return result.edges.findIndex(([a, b]) => (a === i && b === j) || (a === j && b === i));
}

// Initial placement by walking the graph.
function initialPlace(node, result, groups) {
  const n = node.atoms.length;
  const pos = new Array(n).fill(null);
  const templates = groups.map((g, i) => atomTemplate(groups, i));
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
      const { T, slots } = templates[i];
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
        const d = edgeLength(node, result, edgeIndex(result, i, j));
        pos[j] = add(pos[i], scale(dir, d));
        if (!seen.has(j)) { seen.add(j); queue.push({ i: j, parent: i, grand: parent, ringSign: -ringSign }); }
      });
      // ionic partners: beside the atom, in the most open direction between its neighbours
      for (const j of groups[i].assoc) {
        if (pos[j] !== null) continue;
        const away = ionDirection(groups[i].nbrs.filter(nb => pos[nb.j] !== null).map(nb => norm(sub(pos[nb.j], pos[i]))));
        pos[j] = add(pos[i], scale(away, ionDistance(node, i, j)));
        if (!seen.has(j)) { seen.add(j); queue.push({ i: j, parent: i, grand: parent, ringSign }); }
      }
    }
    // next piece a few angstrom to the right
    const maxX = Math.max(...piece.map(i => pos[i][0]));
    offsetX = maxX + 4;
  }
  return pos;
}

// Ring rule. A flat ring of n atoms has interior angles of 180 - 360/n, so the angle between the two
// ring bonds of a ring atom aims at that value when it is below the template angle (108 in a
// five-ring, 90 in a four-ring, 60 in a three-ring), and also when it is above it but every atom of
// the ring is flat (128.6 in tropylium's seven-ring; a seven-ring with a tetrahedral atom keeps 120
// and folds instead). Rings of up to seven atoms count, the shortest ring through the two bonds.
// Returns the angle in radians, or null when the pair a-i-b is not such a ring pair.
function ringAngle(groups, i, a, b, template) {
  const ring = ringThroughAngle(groups, a, i, b);
  if (!ring || ring.size > 7) return null;
  const polygon = Math.PI - 2 * Math.PI / ring.size;
  if (polygon < template - 1e-9) return polygon;
  if (groups[i].groups === 3 && [...ring].every(x => groups[x].groups === 3)) return polygon;
  return null;
}

// Template directions and slot assignment for one atom. A flat atom whose two ring bonds fall under
// the ring rule gets its own flat template with the ring angle between those two slots and the third
// slot opening the rest of the plane (126 beside the 108 of a five-ring); a bridgehead between two
// such rings gets both ring angles. Everything else uses the shared templates. The placement walks
// with these directions, so a five- or seven-ring starts out as a regular polygon.
function atomTemplate(groups, i) {
  const g = groups[i];
  const plain = { T: TEMPLATES[g.groups], slots: assignSlots(g) };
  const flat = ang => [Math.cos(ang), Math.sin(ang), 0];
  if (inSmallRing(groups, i)) {
    // a linear atom held bent by a small ring (benzyne): its two bonds at the ring's angle
    const [a, b] = g.nbrs.map(x => x.j);
    const ang = ringAngle(groups, i, a, b, Math.PI);
    if (ang !== null) return { T: [flat(0), flat(ang)], slots: [a, b] };
  }
  if (g.groups !== 3) return plain;
  const nb = g.nbrs.map(x => x.j);
  const ringPairs = [];
  for (let x = 0; x < nb.length; x++) for (let y = x + 1; y < nb.length; y++) {
    const ang = ringAngle(groups, i, nb[x], nb[y], 2 * Math.PI / 3);
    if (ang !== null) ringPairs.push({ a: nb[x], b: nb[y], ang });
  }
  if (ringPairs.length === 1) {
    const { a, b, ang } = ringPairs[0];
    const third = nb.find(j => j !== a && j !== b);
    return { T: [flat(0), flat(ang), flat(Math.PI + ang / 2)], slots: [a, b, third === undefined ? 'lp' : third] };
  }
  if (ringPairs.length === 2) {
    const [p, q] = ringPairs;
    const shared = [p.a, p.b].find(j => j === q.a || j === q.b);
    if (shared === undefined) return plain;
    const pa = p.a === shared ? p.b : p.a, qa = q.a === shared ? q.b : q.a;
    return { T: [flat(0), flat(p.ang), flat(-q.ang)], slots: [shared, pa, qa] };
  }
  return plain;
}

// Angle targets at atom i, one per pair of neighbours, from the template angles and two rules:
// the ring rule above (for a tetrahedral atom it is applied here; for a flat atom it is already in
// the template), and the lone pair rule. A lone pair that sits in one of the atom's slots takes more
// room than a bond and narrows every angle between the atom's bonds by 2.5 degrees per pair, counted
// with the weights of the structures: water 104.5, ammonia 107, pyridine's C-N-C 117.5, furan's
// C-O-C 105 (108 minus 3 for its 1.2 pairs). A pair with no slot (the amide nitrogen's pair, which
// is in the p orbital above the plane) narrows nothing. Atoms with five or six groups keep their
// templates.
const NARROWING = 2.5 * DEG;
function angleTargets(groups, i, T, slots) {
  const g = groups[i];
  const pairs = []; // {a, b, theta0}
  for (let s = 0; s < T.length; s++) for (let t = s + 1; t < T.length; t++) {
    const a = slots[s], b = slots[t];
    if (a === null || b === null || a === 'lp' || b === 'lp') continue;
    pairs.push({ a, b, theta0: Math.acos(Math.max(-1, Math.min(1, dot(T[s], T[t])))) });
  }
  // a linear atom in a small ring (benzyne's triple-bond carbons) takes the ring's angle too: the
  // ring cannot straighten it (the Lewis ring rule, core/lewis.js)
  if (g.groups === 4 || g.groups === 2) {
    for (const p of pairs) {
      const ring = ringAngle(groups, i, p.a, p.b, p.theta0);
      if (ring !== null) p.theta0 = ring;
    }
  }
  if ((g.groups === 3 || g.groups === 4) && g.lonePairs > 0) pairs.forEach(p => { p.theta0 -= NARROWING * g.pairs; });
  return pairs;
}

// A spectator ion is held beside its partner by a tether this fraction as stiff as a bond: enough to
// keep it there, never enough to bend the partner (the Na+ of the SN2 state is tied to both the O
// and the Cl and must not pull the O-C-Cl axis off straight; the pushes keep it clear of the rest).
const LOOSE = 0.01;

// How far a spectator ion sits from its partner: the ion at its ionic radius, the partner at its own
// ionic radius if it is an ion too (Na+ to Cl- 2.83) and otherwise at its covalent radius plus 0.7,
// which is about what the first-row atoms gain on becoming anions (O 0.66 to O2- 1.40). Na+ sits
// 2.4 from an O-, Cl- 3.3 from a C+ and 3.2 from the N of ammonium, so a loose ion never comes as
// close as a bond would. Covalent radii alone put Cl- 2.0 from a carbocation, which read as a bond.
function ionDistance(node, i, j) {
  const size = k => {
    const e = element(node.atoms[k].el);
    return node.fixed[k] && e.ionRadius ? e.ionRadius : e.radius + 0.7;
  };
  return size(i) + size(j);
}

// Where a spectator ion starts beside an atom with bonds in the given unit directions: away from a
// single bond, away from the pair of a bent atom, and opposite one of the bonds of an atom with three
// or more, which puts it between the others: in the plane of a flat CH3+ (so the page can look at the
// face of the cation with its Cl- beside it) and through a face of a tetrahedral NH4+. It never starts
// on top of a hydrogen, and the pushes settle where it ends.
function ionDirection(dirs) {
  if (!dirs.length) return [0, 1, 0];
  if (dirs.length >= 3) return scale(dirs[0], -1);
  const away = scale(dirs.reduce((s, d) => add(s, d), [0, 0, 0]), -1);
  if (len(away) > 1e-6) return norm(away);
  const d = dirs[0]; // two bonds on one line: any direction across it
  return norm(cross(d, Math.abs(d[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
}

// The spring terms, built once per molecule from the rules above. Positions are handled as one flat
// array x = [x0, y0, z0, x1, y1, z1, ...] so the relaxation can work on it directly.
function springTerms(node, result, groups) {
  const n = node.atoms.length;
  const sym = node.atoms.map(a => a.el);
  const kb = 100, ka = 30, kf = 10, kr = 5;

  const bonds = []; // {i, j, d0, k}
  result.edges.forEach(([i, j], k) => {
    const o = result.avg.orders[k];
    const transition = !!(result.transitionEdges && result.transitionEdges[k]);
    // a bond being made or broken pulls weakly while its order is small (full strength from 1/2 up)
    if (o > 1e-9) bonds.push({ i, j, d0: edgeLength(node, result, k), order: o, k: kb * (transition ? Math.min(1, Math.max(0.1, 2 * o)) : 1) });
    else if (node.fixed[i] || node.fixed[j]) bonds.push({ i, j, d0: ionDistance(node, i, j), order: 0, k: LOOSE * kb });
  });

  // angle terms from template slots; the target can be narrowed by lone pairs or set by a ring
  const angles = []; // {i, a, b, theta0, k}
  groups.forEach((g, i) => {
    const { T, slots } = atomTemplate(groups, i);
    for (const p of angleTargets(groups, i, T, slots)) angles.push({ i, a: p.a, b: p.b, theta0: p.theta0, k: ka });
  });

  // flatness across bonds between two flat (three-group) centres: strong for double or partial-double
  // bonds, weak for a single bond joining two flat centres (conjugation, butadiene's middle bond)
  // A two-group atom in a small ring counts as flat here: its two pi parts are one above and below
  // the ring and one in its plane, so it keeps its neighbours in one plane (benzyne).
  const flats = []; // {a, i, j, b, k}: dihedral a-i-j-b, flat means sin = 0
  const flatLike = i => groups[i].groups === 3 || inSmallRing(groups, i);
  bonds.forEach(({ i, j, order }) => {
    if (!flatLike(i) || !flatLike(j)) return;
    const k = kf * (order > 1.05 ? 1 : 0.3);
    for (const na of groups[i].nbrs) for (const nb of groups[j].nbrs) {
      if (na.j === j || nb.j === i) continue;
      flats.push({ a: na.j, i, j, b: nb.j, k });
    }
  });

  // at a linear atom with two multiple bonds (allene, ketene) the planes of its two neighbours are at
  // right angles, because the two pi parts are; not in a small ring, where one pi part lies in the
  // ring's plane and the other above and below it
  const twists = []; // {a, i, j, b, k}: dihedral a-i-j-b about the i...j axis, right angle means cos = 0
  groups.forEach((g, y) => {
    if (g.groups !== 2 || g.nbrs.length !== 2 || inSmallRing(groups, y)) return;
    const [x, z] = g.nbrs.map(nb => nb.j);
    const ox = g.nbrs[0].order, oz = g.nbrs[1].order;
    if (ox <= 1.05 || oz <= 1.05) return;
    for (const na of groups[x].nbrs) for (const nb of groups[z].nbrs) {
      if (na.j === y || nb.j === y) continue;
      twists.push({ a: na.j, i: x, j: z, b: nb.j, k: kf });
    }
  });

  // non-bonded push between atoms that share no bond and no angle term
  const excluded = new Uint8Array(n * n);
  for (const b of bonds) { excluded[b.i * n + b.j] = 1; excluded[b.j * n + b.i] = 1; }
  for (const x of angles) { excluded[x.a * n + x.b] = 1; excluded[x.b * n + x.a] = 1; }
  const pushes = []; // {i, j, d0, k}
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    if (excluded[i * n + j]) continue;
    pushes.push({ i, j, d0: 1.6 * (element(sym[i]).radius + element(sym[j]).radius), k: kr });
  }
  return { n, bonds, angles, flats, twists, pushes };
}

// Energy and its exact derivative for the spring terms, on a flat position array.
// energy(x) returns E; gradient(x, g) fills g with dE/dx and returns E.
export function energyModel(node, result, groups = electronGroups(node, result)) {
  const terms = springTerms(node, result, groups);
  const { bonds, angles, flats, twists, pushes } = terms;

  // cosine of the dihedral a-i-j-b and, when asked, its derivative with respect to the four atoms
  function dihedralCos(x, a, i, j, b, want) {
    const b1x = x[3 * i] - x[3 * a], b1y = x[3 * i + 1] - x[3 * a + 1], b1z = x[3 * i + 2] - x[3 * a + 2];
    const b2x = x[3 * j] - x[3 * i], b2y = x[3 * j + 1] - x[3 * i + 1], b2z = x[3 * j + 2] - x[3 * i + 2];
    const b3x = x[3 * b] - x[3 * j], b3y = x[3 * b + 1] - x[3 * j + 1], b3z = x[3 * b + 2] - x[3 * j + 2];
    const n1 = [b1y * b2z - b1z * b2y, b1z * b2x - b1x * b2z, b1x * b2y - b1y * b2x];
    const n2 = [b2y * b3z - b2z * b3y, b2z * b3x - b2x * b3z, b2x * b3y - b2y * b3x];
    const l1 = len(n1), l2 = len(n2);
    if (l1 < 1e-9 || l2 < 1e-9) return null;
    const c = dot(n1, n2) / (l1 * l2);
    if (!want) return { c };
    // dc/dn1, dc/dn2, then back through the cross products to the bond vectors
    const g1 = sub(scale(n2, 1 / (l1 * l2)), scale(n1, c / (l1 * l1)));
    const g2 = sub(scale(n1, 1 / (l1 * l2)), scale(n2, c / (l2 * l2)));
    const b1 = [b1x, b1y, b1z], b2 = [b2x, b2y, b2z], b3 = [b3x, b3y, b3z];
    const d1 = cross(b2, g1);
    const d2 = add(cross(g1, b1), cross(b3, g2));
    const d3 = cross(g2, b2);
    return { c, da: scale(d1, -1), di: sub(d1, d2), dj: sub(d2, d3), db: d3 };
  }

  function evaluate(x, g) {
    let E = 0;
    if (g) g.fill(0);
    for (const t of bonds) {
      const dx = x[3 * t.i] - x[3 * t.j], dy = x[3 * t.i + 1] - x[3 * t.j + 1], dz = x[3 * t.i + 2] - x[3 * t.j + 2];
      const d = Math.hypot(dx, dy, dz);
      E += t.k * (d - t.d0) ** 2;
      if (!g || d < 1e-12) continue;
      const f = 2 * t.k * (d - t.d0) / d;
      g[3 * t.i] += f * dx; g[3 * t.i + 1] += f * dy; g[3 * t.i + 2] += f * dz;
      g[3 * t.j] -= f * dx; g[3 * t.j + 1] -= f * dy; g[3 * t.j + 2] -= f * dz;
    }
    for (const t of angles) {
      const u = [x[3 * t.a] - x[3 * t.i], x[3 * t.a + 1] - x[3 * t.i + 1], x[3 * t.a + 2] - x[3 * t.i + 2]];
      const v = [x[3 * t.b] - x[3 * t.i], x[3 * t.b + 1] - x[3 * t.i + 1], x[3 * t.b + 2] - x[3 * t.i + 2]];
      const lu = len(u), lv = len(v);
      const c = Math.max(-1, Math.min(1, dot(u, v) / (lu * lv + 1e-12)));
      const th = Math.acos(c);
      E += t.k * (th - t.theta0) ** 2;
      if (!g) continue;
      const s = Math.sqrt(1 - c * c);
      if (s < 1e-8) continue; // at a straight angle the direction of steepest change is undefined and the slope is zero
      const f = -2 * t.k * (th - t.theta0) / s; // dE/dc
      for (let q = 0; q < 3; q++) {
        const ga = f * (v[q] / (lu * lv) - c * u[q] / (lu * lu));
        const gb = f * (u[q] / (lu * lv) - c * v[q] / (lv * lv));
        g[3 * t.a + q] += ga; g[3 * t.b + q] += gb; g[3 * t.i + q] -= ga + gb;
      }
    }
    for (const t of flats) {
      const r = dihedralCos(x, t.a, t.i, t.j, t.b, !!g);
      if (!r) continue;
      E += t.k * (1 - r.c * r.c); // sin squared of the dihedral
      if (!g) continue;
      const f = -2 * t.k * r.c;
      for (let q = 0; q < 3; q++) { g[3 * t.a + q] += f * r.da[q]; g[3 * t.i + q] += f * r.di[q]; g[3 * t.j + q] += f * r.dj[q]; g[3 * t.b + q] += f * r.db[q]; }
    }
    for (const t of twists) {
      const r = dihedralCos(x, t.a, t.i, t.j, t.b, !!g);
      if (!r) continue;
      E += t.k * r.c * r.c;
      if (!g) continue;
      const f = 2 * t.k * r.c;
      for (let q = 0; q < 3; q++) { g[3 * t.a + q] += f * r.da[q]; g[3 * t.i + q] += f * r.di[q]; g[3 * t.j + q] += f * r.dj[q]; g[3 * t.b + q] += f * r.db[q]; }
    }
    for (const t of pushes) {
      const dx = x[3 * t.i] - x[3 * t.j], dy = x[3 * t.i + 1] - x[3 * t.j + 1], dz = x[3 * t.i + 2] - x[3 * t.j + 2];
      const d = Math.hypot(dx, dy, dz);
      if (d >= t.d0) continue;
      E += t.k * (t.d0 - d) ** 2;
      if (!g || d < 1e-12) continue;
      const f = -2 * t.k * (t.d0 - d) / d;
      g[3 * t.i] += f * dx; g[3 * t.i + 1] += f * dy; g[3 * t.i + 2] += f * dz;
      g[3 * t.j] -= f * dx; g[3 * t.j + 1] -= f * dy; g[3 * t.j + 2] -= f * dz;
    }
    return E;
  }
  return { terms, energy: x => evaluate(x, null), gradient: (x, g) => evaluate(x, g) };
}

// Minimise the model energy from x by the limited-memory quasi-Newton method (L-BFGS): each step
// goes along a direction built from the last few steps and gradient changes, and a backtracking
// line search keeps every step a real decrease. No step moves any coordinate more than maxMove.
function minimize(model, x, { maxIter = 500, gradTol = 1e-3, history = 8, maxMove = 0.3 } = {}) {
  const n = x.length;
  const dotv = (a, b) => { let s = 0; for (let q = 0; q < n; q++) s += a[q] * b[q]; return s; };
  let g = new Float64Array(n), gNew = new Float64Array(n);
  let xNew = new Float64Array(n);
  const d = new Float64Array(n);
  let E = model.gradient(x, g);
  const S = [], Y = [], rho = [];
  let it = 0;
  for (; it < maxIter; it++) {
    const gnorm = Math.sqrt(dotv(g, g));
    if (gnorm < gradTol) break;
    // two-loop recursion: d = H g, where H approximates the inverse curvature
    d.set(g);
    const alpha = new Array(S.length);
    for (let k = S.length - 1; k >= 0; k--) {
      alpha[k] = rho[k] * dotv(S[k], d);
      for (let q = 0; q < n; q++) d[q] -= alpha[k] * Y[k][q];
    }
    if (S.length) {
      const k = S.length - 1;
      const gamma = dotv(S[k], Y[k]) / dotv(Y[k], Y[k]);
      for (let q = 0; q < n; q++) d[q] *= gamma;
    } else {
      for (let q = 0; q < n; q++) d[q] *= 0.02 / gnorm; // first step: a small move along the gradient
    }
    for (let k = 0; k < S.length; k++) {
      const beta = rho[k] * dotv(Y[k], d);
      for (let q = 0; q < n; q++) d[q] += (alpha[k] - beta) * S[k][q];
    }
    let slope = -dotv(g, d);
    if (slope >= 0) { // the history gave an uphill direction: forget it and go down the gradient
      S.length = 0; Y.length = 0; rho.length = 0;
      for (let q = 0; q < n; q++) d[q] = g[q] * 0.02 / gnorm;
      slope = -dotv(g, d);
    }
    let biggest = 0;
    for (let q = 0; q < n; q++) biggest = Math.max(biggest, Math.abs(d[q]));
    let step = biggest > maxMove ? maxMove / biggest : 1;
    let ENew = Infinity, accepted = false;
    for (let tries = 0; tries < 40; tries++) {
      for (let q = 0; q < n; q++) xNew[q] = x[q] - step * d[q];
      ENew = model.gradient(xNew, gNew);
      if (ENew <= E + 1e-4 * step * slope) { accepted = true; break; }
      step *= 0.5;
    }
    if (!accepted) break;
    const s = new Float64Array(n), y = new Float64Array(n);
    for (let q = 0; q < n; q++) { s[q] = xNew[q] - x[q]; y[q] = gNew[q] - g[q]; }
    const sy = dotv(s, y);
    if (sy > 1e-12) {
      S.push(s); Y.push(y); rho.push(1 / sy);
      if (S.length > history) { S.shift(); Y.shift(); rho.shift(); }
    }
    [x, xNew] = [xNew, x];
    [g, gNew] = [gNew, g];
    E = ENew;
  }
  return { x, energy: E, iterations: it };
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
  const model = energyModel(node, result, groups);
  const x = new Float64Array(pos.flat());
  const done = minimize(model, x, { maxIter: opts.iters ?? 500 });
  pos = pos.map((_, i) => [done.x[3 * i], done.x[3 * i + 1], done.x[3 * i + 2]]);
  // centre on the heavy atoms; a molecule relaxed from given positions (an edit) stays where it was,
  // so the picture does not shift with every change
  if (!opts.start) {
    const heavy = node.atoms.map((a, i) => i).filter(i => node.atoms[i].el !== 'H');
    const ref = heavy.length ? heavy : node.atoms.map((_, i) => i);
    const c = scale(ref.reduce((s, i) => add(s, pos[i]), [0, 0, 0]), 1 / ref.length);
    pos = pos.map(p => sub(p, c));
  }
  return { pos, energy: done.energy, groups, iterations: done.iterations };
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
