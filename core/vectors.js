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
