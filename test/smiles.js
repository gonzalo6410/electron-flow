// node test/smiles.js   reading SMILES: aromatic rings with heteroatoms, hydrogens in brackets, charges,
// isotopes (dropped), stereo marks (dropped), several pieces, bad strings (a plain error, not a
// different molecule), elements not in the table; and that pieces written apart keep their own
// electrons (a lone pair does not jump from an amide to a carbocation in another piece).
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';

let bad = 0;
const check = (what, ok, detail = '') => { console.log(`${ok ? 'ok  ' : 'FAIL'} ${what}${detail ? ': ' + detail : ''}`); if (!ok) bad++; };
const formula = spec => { const c = {}; spec.atoms.forEach(a => { c[a.el] = (c[a.el] || 0) + 1; }); return Object.entries(c).sort().map(([e, n]) => e + n).join(' '); };
const run = s => { const node = buildNode(nodeFromSmiles(s)); return { node, r: enumerate(node) }; };
const charges = r => r.avg.charge.map(q => Math.round(q * 1000) / 1000);

// reading
check('pyridine: aromatic ring with N, one H per carbon', formula(nodeFromSmiles('c1ccncc1')) === 'C5 H5 N1');
check('pyrrole: [nH] keeps its hydrogen', formula(nodeFromSmiles('c1cc[nH]c1')) === 'C4 H5 N1');
check('furan and thiophene read', formula(nodeFromSmiles('c1ccoc1')) === 'C4 H4 O1' && formula(nodeFromSmiles('c1ccsc1')) === 'C4 H4 S1');
check('imidazole reads', formula(nodeFromSmiles('c1cnc[nH]1')) === 'C3 H4 N2');
check('[NH4+] keeps four hydrogens and the charge', formula(nodeFromSmiles('[NH4+]')) === 'H4 N1' && nodeFromSmiles('[NH4+]').charge === 1);
check('[13CH4] drops the isotope mark', formula(nodeFromSmiles('[13CH4]')) === 'C1 H4');
check('[2H]C reads deuterium as hydrogen', formula(nodeFromSmiles('[2H]C')) === 'C1 H4');
check('stereo marks are dropped', formula(nodeFromSmiles('F/C=C\\F')) === 'C2 F2 H2' && formula(nodeFromSmiles('N[C@@H](C)C(=O)O')) === 'C3 H7 N1 O2');
{
  const spec = nodeFromSmiles('[Na+].[O-]C(=O)C');
  check('two pieces: Na+ tied loosely to the O written negative', spec.adj.some(([i, j]) => (spec.atoms[i].el === 'Na' && spec.atoms[j].el === 'O') || (spec.atoms[j].el === 'Na' && spec.atoms[i].el === 'O')));
  check('a single-atom anion piece is a spectator ion', nodeFromSmiles('[CH3+].[Cl-]').atoms.find(a => a.el === 'Cl').ion === -1);
  const two = nodeFromSmiles('[O-]C(=O)[O-].[Na+].[Na+]');
  const partners = two.adj.filter(([i, j]) => two.atoms[i].el === 'Na' || two.atoms[j].el === 'Na').map(([i, j]) => (two.atoms[i].el === 'Na' ? j : i));
  check('the two Na+ of carbonate sit by two different O', partners.length === 2 && partners[0] !== partners[1], `partners ${partners}`);
}
for (const [s, why] of [['C(', 'unmatched parenthesis'], ['C)', 'closing parenthesis'], ['C1CC', 'dangling ring closure'], ['hello', 'not an element'], ['*', 'wildcard atom'], ['=C', 'bond with nothing before it'], ['CC=', 'bond with nothing after it'], ['[CH3', 'unmatched bracket'], ['', 'empty']]) {
  let msg = null;
  try { nodeFromSmiles(s); } catch (e) { msg = e.message; }
  check(`bad SMILES ${JSON.stringify(s)} (${why}) gives a plain error`, msg !== null && !/undefined|null|TypeError/.test(msg), msg || 'no error');
}
{
  let msg = null;
  try { nodeFromSmiles('[Fe]Cl'); } catch (e) { msg = e.message; }
  check('an element outside the table is named in the error', msg === 'element Fe is not in the table yet', msg);
}

// pieces keep their electrons
{
  const { r } = run('[CH3+].[CH3-]');
  check('CH3+ and CH3- written apart stay a cation and an anion', r.structures.length === 1 && JSON.stringify(charges(r).slice(0, 2)) === '[1,-1]', JSON.stringify(charges(r)));
}
{
  const { r } = run('[CH3].[CH3]');
  check('two methyl radicals written apart stay two radicals', r.radicals === 2 && r.structures.length === 1 && r.structures[0].lone[0] === 1 && r.structures[0].lone[1] === 1);
}
{
  const { node, r } = run('C[C+](C)C.CC(=O)N');
  const cPlus = node.atoms.findIndex((a, i) => a.el === 'C' && node.edges.filter(([p, q]) => p === i || q === i).filter(([p, q]) => node.atoms[p === i ? q : p].el === 'C').length === 3);
  check('tert-butyl cation beside acetamide: the amide pair stays in the amide', r.structures.length === 2 && r.structures.every(st => st.charge[cPlus] === 1), `structures ${r.structures.length}, C+ charges ${r.structures.map(st => st.charge[cPlus])}`);
}
{
  const { r } = run('[NH4+].CC(=O)[O-]');
  check('ammonium acetate: ammonium +1, acetate two equal structures', r.structures.length === 2 && charges(r)[0] === 1);
}
{
  const { r } = run('C[Mg]Br');
  check('a Grignard written as one piece: Mg2+, C-, Br-', JSON.stringify(charges(r).slice(0, 3)) === '[-1,2,-1]', JSON.stringify(charges(r)));
}
{
  const { r } = run('[Na]Cl');
  check('NaCl written with a bond: Na+ and Cl-', JSON.stringify(charges(r)) === '[1,-1]', JSON.stringify(charges(r)));
}
{
  const { r } = run('O=C=O.O=C=O');
  check('two CO2 pieces: three structures each, nine together', r.structures.length === 9);
}

console.log(bad ? `${bad} problem(s)` : 'all SMILES checks pass');
process.exit(bad ? 1 : 0);
