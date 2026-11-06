// node test/rdkit-dump.js   prints, as JSON, what the engine finds for every named compound and for
// a set of textbook resonance cases, so test/rdkit-compare.py can hold it against RDKit. Each
// structure is given as its atoms (element, formal charge, lone electrons) and its bonds with orders;
// pairs cut to order 0 are left out of the bond list.
import { COMPOUNDS } from '../core/compounds.js';
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { structureFamilies } from '../core/symmetry.js';

// textbook cases beyond the library (the SMILES only says which atoms are joined)
export const TEXTBOOK = [
  ['acetate', 'CC(=O)[O-]'],
  ['formate', 'C(=O)[O-]'],
  ['enolate of acetaldehyde', 'C=C[O-]'],
  ['enolate of acetone', 'CC(=C)[O-]'],
  ['allyl cation', 'C=C[CH2+]'],
  ['allyl anion', 'C=C[CH2-]'],
  ['allyl radical', 'C=C[CH2]'],
  ['benzyl cation', '[CH2+]c1ccccc1'],
  ['benzyl anion', '[CH2-]c1ccccc1'],
  ['benzyl radical', '[CH2]c1ccccc1'],
  ['nitromethane', 'C[N+](=O)[O-]'],
  ['nitrate', '[O-][N+](=O)[O-]'],
  ['nitrite', '[O-]N=O'],
  ['methanesulfonate', 'CS(=O)(=O)[O-]'],
  ['phosphate', '[O-]P(=O)([O-])[O-]'],
  ['azide', '[N-]=[N+]=[N-]'],
  ['hydrazoic acid', 'N=[N+]=[N-]'],
  ['diazomethane', 'C=[N+]=[N-]'],
  ['ozone', '[O-][O+]=O'],
  ['sulfur dioxide', 'O=S=O'],
  ['sulfur trioxide', 'O=S(=O)=O'],
  ['perchlorate', '[O-]Cl(=O)(=O)=O'],
  ['pyridine', 'c1ccncc1'],
  ['pyridinium', 'c1cc[nH+]cc1'],
  ['pyrrole', 'c1cc[nH]c1'],
  ['furan', 'c1ccoc1'],
  ['thiophene', 'c1ccsc1'],
  ['imidazole', 'c1cnc[nH]1'],
  ['naphthalene', 'c1ccc2ccccc2c1'],
  ['cyclopentadienyl anion', '[CH-]1C=CC=C1'],
  ['tropylium', '[CH+]1C=CC=CC=C1'],
  ['guanidinium', 'NC(N)=[NH2+]'],
  ['formamidinium', 'NC=[NH2+]'],
  ['acetylacetonate', 'CC(=O)[CH-]C(C)=O'],
  ['nitrosonium', '[O+]#N'],
  ['nitronium', 'O=[N+]=O'],
  ['cyanate', '[O-]C#N'],
  ['methyl isocyanate', 'CN=C=O'],
  ['thiocyanate', '[S-]C#N'],
  ['fulminate', '[C-]#[N+][O-]'],
  ['ketene', 'C=C=O'],
  ['carbon suboxide', 'O=C=C=C=O'],
  ['butatriene', 'C=C=C=C'],
  ['vinylacetylene', 'C=CC#C'],
  ['carbon disulfide', 'S=C=S'],
  ['methyl isocyanide', 'C[N+]#[C-]'],
  ['urea', 'NC(N)=O'],
  ['carbonate', '[O-]C(=O)[O-]'],
  ['oxalate', '[O-]C(=O)C(=O)[O-]'],
  ['vinylamine (enamine)', 'C=CN'],
  ['methyl vinyl ether', 'C=COC'],
  ['chlorobenzene', 'Clc1ccccc1'],
  ['fluoroethylene', 'C=CF'],
  ['2-pyridone', 'O=c1cccc[nH]1'],
  ['acylium', 'CC#[O+]'],
  ['protonated formaldehyde', 'C=[OH+]'],
  ['methyleneiminium', 'C=[NH2+]'],
  ['hypochlorite', '[O-]Cl'],
  ['cyanamide', 'NC#N'],
  ['carbodiimide', 'N=C=N'],
  ['superoxide', '[O][O-]'],
  ['hydroxyl radical', '[OH]'],
  ['nitrogen dioxide', '[O]N=O'],
  ['nitric oxide', '[N]=O'],
  ['methyl radical', '[CH3]'],
];

export function dumpOne(name, smiles) {
  const out = { name, smiles };
  try {
    const spec = nodeFromSmiles(smiles);
    const node = buildNode(spec);
    const r = enumerate(node);
    const fam = structureFamilies(node, r);
    out.atoms = node.atoms.map(a => a.el);
    out.charge = spec.charge;
    out.pool = r.pool;
    out.holes = r.holes;
    out.radicals = r.radicals;
    out.families = fam.families;
    out.structures = r.structures.map(st => ({
      charge: st.charge,
      lone: st.lone,
      bonds: r.edges.map(([i, j], k) => [i, j, st.orders[k]]).filter(b => b[2] > 0),
    }));
  } catch (e) {
    out.error = e.message;
  }
  return out;
}

if (process.argv[1] && process.argv[1].endsWith('rdkit-dump.js')) {
  const all = [
    ...COMPOUNDS.map(c => ['library: ' + c.name, c.smiles]),
    ...TEXTBOOK.map(([n, s]) => ['textbook: ' + n, s]),
  ];
  console.log(JSON.stringify(all.map(([n, s]) => dumpOne(n, s))));
}
