// node test/edit.js   builds molecules by editing operations and checks the engine's answer
import { emptySpec, addAtom, removeAtom, toggleBond, setCharge, fillHydrogens, placeNear } from '../core/edit.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { relax } from '../core/geometry.js';

let bad = 0;
const check = (name, cond) => { console.log(`${cond ? 'ok  ' : 'FAIL'} ${name}`); if (!cond) bad++; };

// ethane from scratch: C, then C on it, then fill
let s = emptySpec();
let r = addAtom(s, 'C'); s = r.spec; const c0 = r.index;
r = addAtom(s, 'C', c0); s = r.spec; const c1 = r.index;
s = fillHydrogens(s).spec;
check('ethane has 8 atoms', s.atoms.length === 8);
let res = enumerate(buildNode(s));
check('ethane: one structure, no holes, C-C single', res.structures.length === 1 && res.holes === 0 && res.structures[0].orders[0] === 1);

// take one hydrogen off each carbon: the engine should find the double bond
const hOnC0 = s.adj.filter(([a, b]) => a === c0 || b === c0).map(([a, b]) => (a === c0 ? b : a)).find(j => s.atoms[j].el === 'H');
s = removeAtom(s, hOnC0);
const hOnC1 = s.adj.filter(([a, b]) => a === c1 || b === c1).map(([a, b]) => (a === c1 ? b : a)).find(j => s.atoms[j].el === 'H');
s = removeAtom(s, hOnC1);
res = enumerate(buildNode(s));
const cc = s.adj.findIndex(([a, b]) => s.atoms[a].el === 'C' && s.atoms[b].el === 'C');
check('ethylene: C=C after removing two hydrogens', res.structures.length === 1 && res.structures[0].orders[cc] === 2 && res.holes === 0);

// take one more hydrogen off: vinyl cation needs a charge; with charge +1 one hole
const hAny = s.atoms.findIndex(a => a.el === 'H');
let t = removeAtom(s, hAny);
let rt = enumerate(buildNode(t));
check('C2H3 neutral is a radical', rt.radicals === 1);
t = setCharge(t, 1);
rt = enumerate(buildNode(t));
check('C2H3+ has one hole and no radical', rt.radicals === 0 && rt.holes === 2);

// ring closure: build a six-carbon chain, close it, fill, and get cyclohexane; remove six H, benzene
let ring = emptySpec();
let prev = null;
const cs = [];
for (let k = 0; k < 6; k++) { const q = addAtom(ring, 'C', prev); ring = q.spec; prev = q.index; cs.push(q.index); }
ring = toggleBond(ring, cs[0], cs[5]);
check('ring closed: six C-C bonds', ring.adj.length === 6);
ring = fillHydrogens(ring).spec;
check('cyclohexane: 18 atoms', ring.atoms.length === 18);
let rr = enumerate(buildNode(ring));
check('cyclohexane: one structure, all single', rr.structures.length === 1 && rr.structures[0].orders.every(o => o === 1));
for (const c of cs) {
  const h = ring.adj.filter(([a, b]) => a === c || b === c).map(([a, b]) => (a === c ? b : a)).find(j => ring.atoms[j].el === 'H');
  // removing shifts indices; recompute carbon indices by element each time is simpler: carbons are the first six
  ring = removeAtom(ring, h);
}
rr = enumerate(buildNode(ring));
check('benzene from cyclohexane minus six H: two structures, no holes', rr.structures.length === 2 && rr.holes === 0);

// toggling the same bond twice restores the spec
const once = toggleBond(ring, 0, 3);
const twice = toggleBond(once, 0, 3);
check('toggle twice restores', twice.adj.length === ring.adj.length);

// placing a new atom near a parent keeps a plausible distance and avoids existing neighbours
const node = buildNode(ring); const geo = relax(node, rr);
const pNew = placeNear(geo.pos, 0, 'Cl', ring);
const dist = Math.hypot(pNew[0] - geo.pos[0][0], pNew[1] - geo.pos[0][1], pNew[2] - geo.pos[0][2]);
check('new atom placed about a bond length away', dist > 1.5 && dist < 2.0);
const freePos = placeNear(geo.pos, null, 'Na', ring);
check('free atom placed to the right of everything', freePos[0] > Math.max(...geo.pos.map(p => p[0])));

console.log(bad ? `${bad} problem(s)` : 'edit tests pass');
process.exit(bad ? 1 : 0);
