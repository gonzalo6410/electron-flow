// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Small helpers for vectors written as arrays of three numbers, shared by the geometry and the
// pi-frame rules. norm of a vector too short to have a direction is the zero vector.

const EPS = 1e-9;
export function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
export function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
export function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
export function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
export function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
export function len(a) { return Math.hypot(a[0], a[1], a[2]); }
export function norm(a) { const l = len(a); return l < EPS ? [0, 0, 0] : scale(a, 1 / l); }
// the part of v at right angles to a unit axis, as a unit vector
export function perpTo(v, axis) { return norm(sub(v, scale(axis, dot(v, axis)))); }
// some unit vector at right angles to a unit axis
export function anyPerp(axis) { return perpTo(Math.abs(axis[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0], axis); }
// v turned by ang about a unit axis (Rodrigues)
export function rotate(v, axis, ang) {
  const c = Math.cos(ang), s = Math.sin(ang);
  return add(add(scale(v, c), scale(cross(axis, v), s)), scale(axis, dot(axis, v) * (1 - c)));
}

// The turn and shift that best lay `points` onto `targets` (each pair counted with its weight), as a
// function that moves any point the same way. Horn's method: the best turn is the eigenvector, with
// the largest eigenvalue, of a 4 by 4 matrix built from the pairs, read as a quaternion. Used to keep
// a molecule still on the screen while the atoms around it move.
export function bestFit(points, targets, weights = points.map(() => 1)) {
  const W = weights.reduce((s, w) => s + w, 0) || 1;
  const centre = ps => [0, 1, 2].map(a => ps.reduce((s, p, k) => s + weights[k] * p[a], 0) / W);
  const cp = centre(points), ct = centre(targets);
  const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  points.forEach((p, k) => {
    const u = sub(p, cp), v = sub(targets[k], ct);
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += weights[k] * u[a] * v[b];
  });
  const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
  const N = [
    [xx + yy + zz, yz - zy, zx - xz, xy - yx],
    [yz - zy, xx - yy - zz, xy + yx, zx + xz],
    [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
    [xy - yx, zx + xz, yz + zy, -xx - yy + zz],
  ];
  const [w, x, y, z] = largestEigenvector(N);
  const R = [
    [w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
    [2 * (y * x + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
    [2 * (z * x - w * y), 2 * (z * y + w * x), w * w - x * x - y * y + z * z],
  ];
  return p => { const u = sub(p, cp); return add(R.map(row => dot(row, u)), ct); };
}

// the unit eigenvector with the largest eigenvalue of a small symmetric matrix, by Jacobi rotations
function largestEigenvector(M) {
  const n = M.length;
  const m = M.map(row => [...row]);
  const v = m.map((_, i) => m.map((__, j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += m[p][q] * m[p][q];
    if (off < 1e-20) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(m[p][q]) < 1e-15) continue;
      const th = 0.5 * Math.atan2(2 * m[p][q], m[q][q] - m[p][p]);
      const c = Math.cos(th), s = Math.sin(th);
      for (let k = 0; k < n; k++) { const a = m[k][p], b = m[k][q]; m[k][p] = c * a - s * b; m[k][q] = s * a + c * b; }
      for (let k = 0; k < n; k++) { const a = m[p][k], b = m[q][k]; m[p][k] = c * a - s * b; m[q][k] = s * a + c * b; }
      for (let k = 0; k < n; k++) { const a = v[k][p], b = v[k][q]; v[k][p] = c * a - s * b; v[k][q] = s * a + c * b; }
    }
  }
  let best = 0;
  for (let k = 1; k < n; k++) if (m[k][k] > m[best][best]) best = k;
  const col = v.map(row => row[best]);
  const l = Math.hypot(...col) || 1;
  return col.map(c => c / l);
}
