// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// Element table for the Lewis engine.
// valence: electrons the neutral atom brings to the pool.
// ceiling: most electrons the atom may count (bonds count two each, lone electrons one each).
// ionic: group 1 and 2 metals; they give their electrons away, form no bonds, are not octet counted (ledger item 24).
// radius: single-bond covalent radius in angstrom (Cordero 2008 values, rounded), used only to set spring lengths.
// ionRadius: radius of the element as a spectator ion (Na+, Cl-, O2-; Shannon crystal radii for six
//   neighbours, rounded; the hydride value is a rough one), used only to set how far a loose ion sits
//   from its partner. Only elements that are written as single-atom ions have one.
// en: Pauling electronegativity, used only for display ordering and later for weighting rules if he wants them.
// color: display colour (CPK-like), renderer only.

export const ELEMENTS = {
  H:  { Z: 1,  valence: 1, ceiling: 2, radius: 0.31, ionRadius: 1.40, en: 2.20, color: 0xffffff },
  He: { Z: 2,  valence: 2, ceiling: 2, radius: 0.28, en: 0,    color: 0xd9ffff },
  Li: { Z: 3,  valence: 1, ceiling: 0, ionic: true, radius: 1.28, ionRadius: 0.76, en: 0.98, color: 0xcc80ff },
  Be: { Z: 4,  valence: 2, ceiling: 0, ionic: true, radius: 0.96, ionRadius: 0.45, en: 1.57, color: 0xc2ff00 },
  B:  { Z: 5,  valence: 3, ceiling: 8, radius: 0.84, en: 2.04, color: 0xffb5b5 },
  C:  { Z: 6,  valence: 4, ceiling: 8, radius: 0.76, en: 2.55, color: 0x909090 },
  N:  { Z: 7,  valence: 5, ceiling: 8, radius: 0.71, ionRadius: 1.46, en: 3.04, color: 0x3050f8 },
  O:  { Z: 8,  valence: 6, ceiling: 8, radius: 0.66, ionRadius: 1.40, en: 3.44, color: 0xff0d0d },
  F:  { Z: 9,  valence: 7, ceiling: 8, radius: 0.57, ionRadius: 1.33, en: 3.98, color: 0x90e050 },
  Ne: { Z: 10, valence: 8, ceiling: 8, radius: 0.58, en: 0,    color: 0xb3e3f5 },
  Na: { Z: 11, valence: 1, ceiling: 0, ionic: true, radius: 1.66, ionRadius: 1.02, en: 0.93, color: 0xab5cf2 },
  Mg: { Z: 12, valence: 2, ceiling: 0, ionic: true, radius: 1.41, ionRadius: 0.72, en: 1.31, color: 0x8aff00 },
  Al: { Z: 13, valence: 3, ceiling: 8, radius: 1.21, en: 1.61, color: 0xbfa6a6 },
  Si: { Z: 14, valence: 4, ceiling: 8, radius: 1.11, en: 1.90, color: 0xf0c8a0 },
  P:  { Z: 15, valence: 5, ceiling: 8, radius: 1.07, en: 2.19, color: 0xff8000 },
  S:  { Z: 16, valence: 6, ceiling: 8, radius: 1.05, ionRadius: 1.84, en: 2.58, color: 0xffff30 },
  Cl: { Z: 17, valence: 7, ceiling: 8, radius: 1.02, ionRadius: 1.81, en: 3.16, color: 0x1ff01f },
  Ar: { Z: 18, valence: 8, ceiling: 8, radius: 1.06, en: 0,    color: 0x80d1e3 },
  K:  { Z: 19, valence: 1, ceiling: 0, ionic: true, radius: 2.03, ionRadius: 1.38, en: 0.82, color: 0x8f40d4 },
  Ca: { Z: 20, valence: 2, ceiling: 0, ionic: true, radius: 1.76, ionRadius: 1.00, en: 1.00, color: 0x3dff00 },
  Se: { Z: 34, valence: 6, ceiling: 8, radius: 1.20, ionRadius: 1.98, en: 2.55, color: 0xffa100 },
  Br: { Z: 35, valence: 7, ceiling: 8, radius: 1.20, ionRadius: 1.96, en: 2.96, color: 0xa62929 },
  Kr: { Z: 36, valence: 8, ceiling: 8, radius: 1.16, en: 3.00, color: 0x5cb8d1 },
  I:  { Z: 53, valence: 7, ceiling: 8, radius: 1.39, ionRadius: 2.20, en: 2.66, color: 0x940094 },
  Xe: { Z: 54, valence: 8, ceiling: 8, radius: 1.40, en: 2.60, color: 0x429eb0 },
};

export function element(symbol) {
  const e = ELEMENTS[symbol];
  if (!e) throw new Error(`unknown element ${symbol}`);
  return e;
}
