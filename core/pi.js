// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Where the pi bonds and the lone pairs point. Rules, nothing stored per molecule.
//
// Pi bonds. A bond of order above 1 has a pi part whose p orbitals point along a unit vector n at
// right angles to the bond. n comes from the flat plane of a three-group atom at either end (the
// normal of that plane). At a linear atom carrying two multiple bonds (allene, CO2, N=N=O) the two
// pi parts are at right angles to each other, so a known n on one side fixes n on the other. A
// triple bond has two pi parts at right angles. With nothing to go on, n is chosen at random at
// right angles to the bond.
//
// Lone pairs. An atom has slots from its electron-group count (item 29): a three-group atom with one
// neighbour has two slots in its plane at 120 degrees from the bond; a four-group atom has the
// tetrahedral slots left free by its neighbours; a two-group atom with one neighbour has the slot
// opposite the bond; six-group with four neighbours has the two slots above and below. Lone pairs
// beyond the slots are the ones taking part in bonding: a pair toward a partner whose bond is
// absent or partial in this picture (the F- of SF6 points its fourth pair at S; the nucleophile of
// the SN2 state points its half pair at carbon), otherwise a pair in the p orbital (the amide N pair,
// the third pair of a carboxylate O), otherwise placed by repulsion.

import { sub, add, scale, dot, cross, len, norm, perpTo, anyPerp, rotate } from './vectors.js';

// pic: { orders, lone } for the picture being drawn (one structure, or the average)
export function piFrame(node, result, pos, groups, pic) {
  const n = node.atoms.length;
  const edges = result.edges;
  const ionic = i => node.fixed[i];

  // neighbours in the average (geometry) and in this picture
  const avgNbrs = groups.map(g => g.nbrs.map(x => x.j));
  const picOrder = (i, j) => { const k = edges.findIndex(([a, b]) => (a === i && b === j) || (a === j && b === i)); return k < 0 ? 0 : pic.orders[k]; };
  const dirTo = (i, j) => norm(sub(pos[j], pos[i]));

  // plane normal of three-group atoms with two or more neighbours
  const planeN = new Array(n).fill(null);
  for (let i = 0; i < n; i++) {
    if (groups[i].groups !== 3 || avgNbrs[i].length < 2) continue;
    const c = cross(dirTo(i, avgNbrs[i][0]), dirTo(i, avgNbrs[i][1]));
    if (len(c) > 1e-6) planeN[i] = norm(c);
  }

  // pi normals per edge: [] for no pi, [n] for one, [n1, n2] for two
  const piCount = edges.map((_, k) => {
    const o = pic.orders[k];
    if (o <= 1 + 1e-6) return 0;
    return o <= 2 + 1e-6 ? 1 : 2;
  });
  const normals = edges.map(() => []);
  const axisOf = k => norm(sub(pos[edges[k][1]], pos[edges[k][0]]));
  const edgesAt = Array.from({ length: n }, () => []);
  edges.forEach(([i, j], k) => { edgesAt[i].push(k); edgesAt[j].push(k); });

  // seed from planes, then borrow from next door, then random
  let pending = edges.map((_, k) => k).filter(k => piCount[k] > 0);
  const setFrom = (k, v) => {
    const axis = axisOf(k);
    const n1 = perpTo(v, axis);
    if (len(n1) < 1e-6) return false;
    normals[k] = piCount[k] === 2 ? [n1, norm(cross(axis, n1))] : [n1];
    return true;
  };
  for (const k of pending) {
    const [i, j] = edges[k];
    const v = planeN[i] || planeN[j];
    if (v) setFrom(k, v);
  }
  // A pi part with no flat atom at either end takes its direction from next door. At a linear atom
  // it is at right angles to the other pi part there (allene, CO2, N#N-O). Across a single bond from
  // a linear atom it lies parallel to a flat atom's p orbitals or to the pi part on the far side
  // (the triple bond of vinylacetylene or benzonitrile, the two of a diyne), so the lines of a
  // triple bond are drawn in the plane of the molecule, not edge-on at random.
  const borrow = k => {
    const [i, j] = edges[k];
    for (const atom of [i, j]) {
      if (groups[atom].groups !== 2) continue;
      const other = edgesAt[atom].find(q => q !== k && normals[q].length);
      if (other !== undefined && setFrom(k, cross(axisOf(k), normals[other][0]))) return true;
    }
    for (const atom of [i, j]) {
      if (groups[atom].groups !== 2) continue;
      for (const q of edgesAt[atom]) {
        if (q === k || piCount[q] > 0) continue;
        const b = edges[q][0] === atom ? edges[q][1] : edges[q][0];
        if (planeN[b] && setFrom(k, planeN[b])) return true;
        const far = edgesAt[b].find(r => r !== q && normals[r].length);
        if (far !== undefined && setFrom(k, normals[far][0])) return true;
      }
    }
    return false;
  };
  const propagate = () => {
    let progress = true;
    while (pending.length && progress) {
      progress = false;
      for (const k of pending) if (!normals[k].length && borrow(k)) progress = true;
      pending = pending.filter(k => !normals[k].length);
    }
  };
  pending = pending.filter(k => !normals[k].length);
  propagate();
  // nothing to go on: a fixed perpendicular, then propagate again
  while (pending.length) {
    setFrom(pending[0], anyPerp(axisOf(pending[0])));
    pending = pending.filter(q => !normals[q].length);
    propagate();
  }

  // plane normals for terminal three-group atoms come from their pi bond; failing that from a
  // neighbour's plane (the conjugated p orbital of a carboxylate O or an amide N)
  const pOrbital = new Array(n).fill(null); // direction of the p orbital available at the atom
  for (let i = 0; i < n; i++) {
    if (planeN[i]) { pOrbital[i] = planeN[i]; continue; }
    const withPi = edgesAt[i].find(k => normals[k].length);
    if (withPi !== undefined) { pOrbital[i] = normals[withPi][0]; continue; }
    const j = avgNbrs[i].find(j => planeN[j]);
    if (j !== undefined) { pOrbital[i] = perpTo(planeN[j], dirTo(i, j)); continue; }
    // beside a linear atom carrying a pi bond on its far side (the O of N#N-O, the O- of O-C#O+):
    // the conjugated p orbital lines up with that pi part
    for (const jj of avgNbrs[i]) {
      if (groups[jj].groups !== 2) continue;
      const far = edgesAt[jj].find(k => normals[k].length && !edges[k].includes(i));
      if (far !== undefined) { pOrbital[i] = perpTo(normals[far][0], dirTo(i, jj)); break; }
    }
  }

  // lone pairs
  const lonePairs = []; // per atom: [{dir, kind}] kind: slot | partner | p | free
  for (let i = 0; i < n; i++) {
    const out = [];
    lonePairs.push(out);
    if (ionic(i)) continue;
    const pairs = pic.lone[i] / 2;
    const want = Math.ceil(pairs - 1e-9);
    if (!want) continue;
    const G = groups[i].groups;
    const nb = avgNbrs[i];
    const d = nb.map(j => dirTo(i, j));
    const slotDirs = [];
    const slots = Math.max(0, G - nb.length);
    if (slots > 0) {
      if (nb.length === 1 && want === 1) slotDirs.push(scale(d[0], -1)); // one neighbour, one pair: opposite (sp)
      else if (G === 2 && nb.length === 1) slotDirs.push(scale(d[0], -1));
      else if (G === 3 && nb.length === 1 && pOrbital[i]) {
        const nrm = pOrbital[i];
        slotDirs.push(rotate(d[0], nrm, 2 * Math.PI / 3), rotate(d[0], nrm, -2 * Math.PI / 3));
      } else if (G === 3 && nb.length === 2) slotDirs.push(norm(scale(add(d[0], d[1]), -1)));
      else if (G === 4 && nb.length === 3) slotDirs.push(norm(scale(add(add(d[0], d[1]), d[2]), -1)));
      else if (G === 4 && nb.length === 2) {
        const bis = norm(scale(add(d[0], d[1]), -1));
        const nrm = norm(cross(d[0], d[1]));
        const a = 54.75 * Math.PI / 180;
        slotDirs.push(add(scale(bis, Math.cos(a)), scale(nrm, Math.sin(a))), add(scale(bis, Math.cos(a)), scale(nrm, -Math.sin(a))));
      } else if (G === 6 && nb.length === 4) {
        let nrm = null;
        for (let a = 0; a < 4 && !nrm; a++) for (let b = a + 1; b < 4 && !nrm; b++) { const c = cross(d[a], d[b]); if (len(c) > 0.5) nrm = norm(c); }
        if (nrm) slotDirs.push(nrm, scale(nrm, -1));
      } else if (G === 5 && nb.length === 4) {
        // seesaw: the lone pair sits in the equatorial plane opposite the two equatorial bonds
        slotDirs.push(norm(scale(add(add(add(d[0], d[1]), d[2]), d[3]), -1)));
      }
      // fill any slots the simple rules did not give by repulsion
      while (slotDirs.length < slots) slotDirs.push(...repulsion(d, slotDirs, slots - slotDirs.length));
    }
    slotDirs.slice(0, Math.min(want, slots)).forEach(dir => out.push({ dir, kind: 'slot' }));
    // pairs beyond the slots
    let extra = want - Math.min(want, slots);
    const partners = nb.filter(j => picOrder(i, j) < 1 - 1e-6); // bond absent or partial in this picture
    for (const j of partners) { if (extra <= 0) break; out.push({ dir: dirTo(i, j), kind: 'partner', partner: j }); extra--; }
    if (extra > 0 && pOrbital[i]) { out.push({ dir: pOrbital[i], kind: 'p' }); extra--; }
    if (extra > 0 && pOrbital[i]) { out.push({ dir: scale(pOrbital[i], -1), kind: 'p' }); extra--; }
    if (extra > 0) repulsion(d, out.map(o => o.dir), extra).forEach(dir => out.push({ dir, kind: 'free' }));
  }
  return { normals, lonePairs, planeN, pOrbital };
}

// unit vectors keeping away from the given directions and each other
export function repulsion(avoid, existing, count) {
  let dirs = [];
  let mean = avoid.concat(existing).reduce((s, v) => add(s, v), [0, 0, 0]);
  const away = len(mean) > 1e-6 ? scale(norm(mean), -1) : [0, 0, 1];
  const p1 = anyPerp(away), p2 = cross(away, p1);
  for (let q = 0; q < count; q++) {
    const ang = (2 * Math.PI * q) / count + 0.3;
    dirs.push(norm(add(away, add(scale(p1, 0.6 * Math.cos(ang)), scale(p2, 0.6 * Math.sin(ang))))));
  }
  for (let it = 0; it < 60; it++) {
    dirs = dirs.map((dv, q) => {
      let f = [0, 0, 0];
      for (const o of avoid.concat(existing)) { const diff = sub(dv, o); const r = len(diff) + 1e-6; f = add(f, scale(diff, 1 / (r * r * r))); }
      dirs.forEach((o, r) => { if (r !== q) { const diff = sub(dv, o); const rr = len(diff) + 1e-6; f = add(f, scale(diff, 1 / (rr * rr * rr))); } });
      return norm(add(dv, scale(f, 0.1)));
    });
  }
  return dirs;
}
