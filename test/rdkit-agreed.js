// Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
// node test/rdkit-agreed.js   known answers for every named compound and the textbook cases of
// test/rdkit-dump.js, recorded on 6 Nov 2026 after holding the engine against RDKit's
// ResonanceMolSupplier (test/rdkit-compare.py) and reading every difference. Each row: name, SMILES,
// holes (electrons), radicals, family sizes, the pictures found (charges, radicals as El., multiple
// bonds, cut pairs; "N x" for N structures that are the same picture with atoms renamed), and what
// the comparison showed: "RDKit agrees", or the rule of ours that explains the difference. A change
// that alters any row is caught here; if the change is meant, re-record the row and say why.
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { structureFamilies } from '../core/symmetry.js';

export const AGREED = [
  ["water", "O", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["ammonia", "N", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["methane", "C", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["ethane", "CC", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["ethylene", "C=C", 0, 0, [1], ["1 x C=C"], "RDKit agrees"],
  ["acetylene", "C#C", 0, 0, [1], ["1 x C#C"], "RDKit agrees"],
  ["allene", "C=C=C", 0, 0, [1], ["1 x C=C C=C"], "RDKit agrees"],
  ["butadiene", "C=CC=C", 0, 0, [1], ["1 x C=C C=C"], "RDKit agrees"],
  ["carbon dioxide", "O=C=O", 0, 0, [2, 1], ["1 x C=O C=O","2 x O+1 O-1 C#O"], "no charge filter: O(-)-C#O(+) kept; RDKit has only O=C=O"],
  ["carbon monoxide", "[C-]#[O+]", 0, 0, [1], ["1 x C-1 O+1 C#O"], "RDKit agrees"],
  ["nitrous oxide", "N#[N+][O-]", 0, 0, [1, 1], ["1 x N+1 N-1 N=N N=O","1 x N+1 O-1 N#N"], "fewest holes: RDKit adds forms with a six-electron N"],
  ["ozone", "[O-][O+]=O", 0, 0, [2], ["2 x O+1 O-1 O=O"], "fewest holes: RDKit adds forms with a six-electron O"],
  ["sulfur dioxide", "O=S=O", 0, 0, [2], ["2 x O-1 S+1 O=S"], "octet ceiling: O=S(+)-O(-), not O=S=O"],
  ["hydrogen peroxide", "OO", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["hydrogen cyanide", "C#N", 0, 0, [1], ["1 x C#N"], "RDKit agrees"],
  ["formaldehyde", "C=O", 0, 0, [1], ["1 x C=O"], "RDKit agrees"],
  ["ketene", "C=C=O", 0, 0, [1, 1], ["1 x C-1 O+1 C#O","1 x C=C C=O"], "no charge filter: CH2(-)-C#O(+) kept (a textbook form RDKit lacks)"],
  ["diazomethane", "C=[N+]=[N-]", 0, 0, [1, 1], ["1 x C-1 N+1 N#N","1 x N+1 N-1 C=N N=N"], "no charge filter: CH2(-)-N(+)#N kept; fewest holes drops RDKit's C=N-N with a six-electron N"],
  ["hydroxide (Na+)", "[OH-].[Na+]", 0, 0, [1], ["1 x Na+1 O-1"], "RDKit agrees"],
  ["hydronium (Cl-)", "[OH3+].[Cl-]", 0, 0, [1], ["1 x Cl-1 O+1"], "RDKit agrees"],
  ["ammonium (Cl-)", "[NH4+].[Cl-]", 0, 0, [1], ["1 x Cl-1 N+1"], "RDKit agrees"],
  ["nitrate (Na+)", "[O-][N+](=O)[O-].[Na+]", 0, 0, [3], ["3 x N+1 Na+1 O-1 O-1 N=O"], "RDKit agrees"],
  ["nitrite (Na+)", "[O-]N=O.[Na+]", 0, 0, [2], ["2 x Na+1 O-1 N=O"], "fewest holes: RDKit adds O(-)-N(+)-O(-) with a six-electron N"],
  ["carbonate (2 Na+)", "[O-]C(=O)[O-].[Na+].[Na+]", 0, 0, [3], ["3 x Na+1 Na+1 O-1 O-1 C=O"], "RDKit agrees"],
  ["sulfate (2 Na+)", "[O-]S(=O)(=O)[O-].[Na+].[Na+]", 0, 0, [1], ["1 x Na+1 Na+1 O-1 O-1 O-1 O-1 S+2"], "octet ceiling: S(2+) with four O(-)"],
  ["phosphate (3 Na+)", "[O-]P(=O)([O-])[O-].[Na+].[Na+].[Na+]", 0, 0, [1], ["1 x Na+1 Na+1 Na+1 O-1 O-1 O-1 O-1 P+1"], "octet ceiling: P(+) with four O(-)"],
  ["perchlorate (Na+)", "[O-]Cl(=O)(=O)=O.[Na+]", 0, 0, [1], ["1 x Cl+3 Na+1 O-1 O-1 O-1 O-1"], "octet ceiling: Cl(3+) with four O(-); RDKit happens to give the same"],
  ["cyanide (Na+)", "[C-]#N.[Na+]", 0, 0, [1], ["1 x C-1 Na+1 C#N"], "RDKit agrees"],
  ["azide (Na+)", "[N-]=[N+]=[N-].[Na+]", 0, 0, [1], ["1 x N+1 N-1 N-1 Na+1 N=N N=N"], "magnitude: N#N(+)-N(2-) dropped; fewest holes drops RDKit's radical and six-electron forms"],
  ["acetate (Na+)", "CC(=O)[O-].[Na+]", 0, 0, [2], ["2 x Na+1 O-1 C=O"], "RDKit agrees"],
  ["formate (Na+)", "C(=O)[O-].[Na+]", 0, 0, [2], ["2 x Na+1 O-1 C=O"], "RDKit agrees"],
  ["enolate of acetaldehyde (Na+)", "C=C[O-].[Na+]", 0, 0, [1, 1], ["1 x C-1 Na+1 C=O","1 x Na+1 O-1 C=C"], "RDKit agrees"],
  ["methyl cation (Cl-)", "[CH3+].[Cl-]", 2, 0, [1], ["1 x C+1 Cl-1"], "RDKit agrees"],
  ["tert-butyl cation (Cl-)", "C[C+](C)C.[Cl-]", 2, 0, [1], ["1 x C+1 Cl-1"], "RDKit agrees"],
  ["allyl cation (Cl-)", "C=C[CH2+].[Cl-]", 2, 0, [2], ["2 x C+1 Cl-1 C=C"], "RDKit agrees"],
  ["allyl anion (Na+)", "C=C[CH2-].[Na+]", 0, 0, [2], ["2 x C-1 Na+1 C=C"], "RDKit agrees"],
  ["pentadienyl cation (Cl-)", "C=CC=C[CH2+].[Cl-]", 2, 0, [2, 1], ["3 x C+1 Cl-1 C=C C=C"], "RDKit agrees"],
  ["benzyl cation (Cl-)", "[CH2+]c1ccccc1.[Cl-]", 2, 0, [2, 2, 1], ["5 x C+1 Cl-1 C=C C=C C=C"], "RDKit agrees"],
  ["phenoxide (Na+)", "[O-]c1ccccc1.[Na+]", 0, 0, [2, 2, 1], ["3 x C-1 Na+1 C=C C=C C=O","2 x Na+1 O-1 C=C C=C C=C"], "RDKit agrees"],
  ["cyclopentadienyl anion (Na+)", "[CH-]1C=CC=C1.[Na+]", 0, 0, [5], ["5 x C-1 Na+1 C=C C=C"], "RDKit agrees"],
  ["tropylium (Cl-)", "[CH+]1C=CC=CC=C1.[Cl-]", 2, 0, [7], ["7 x C+1 Cl-1 C=C C=C C=C"], "RDKit agrees"],
  ["methanol", "CO", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["ethanol", "CCO", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["dimethyl ether", "COC", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["acetone", "CC(=O)C", 0, 0, [1], ["1 x C=O"], "RDKit agrees"],
  ["acetaldehyde", "CC=O", 0, 0, [1], ["1 x C=O"], "RDKit agrees"],
  ["acetic acid", "CC(=O)O", 0, 0, [1, 1], ["1 x C=O","1 x O+1 O-1 C=O"], "no charge filter: O(-)-C=O(+)H kept"],
  ["methyl acetate", "CC(=O)OC", 0, 0, [1, 1], ["1 x C=O","1 x O+1 O-1 C=O"], "no charge filter: O(-)-C=O(+)Me kept"],
  ["formamide", "C(=O)N", 0, 0, [1, 1], ["1 x C=O","1 x N+1 O-1 C=N"], "RDKit agrees"],
  ["acetamide", "CC(=O)N", 0, 0, [1, 1], ["1 x C=O","1 x N+1 O-1 C=N"], "RDKit agrees"],
  ["acetonitrile", "CC#N", 0, 0, [1], ["1 x C#N"], "RDKit agrees"],
  ["nitromethane", "C[N+](=O)[O-]", 0, 0, [2], ["2 x N+1 O-1 N=O"], "RDKit agrees"],
  ["methyl vinyl ether", "C=COC", 0, 0, [1, 1], ["1 x C-1 O+1 C=O","1 x C=C"], "no charge filter: CH2(-)-CH=O(+)Me kept"],
  ["acrolein", "C=CC=O", 0, 0, [1], ["1 x C=C C=O"], "RDKit agrees"],
  ["methyl vinyl ketone", "C=CC(=O)C", 0, 0, [1], ["1 x C=C C=O"], "RDKit agrees"],
  ["chloromethane", "CCl", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["cyclohexane", "C1CCCCC1", 0, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["cyclohexene", "C1CCC=CC1", 0, 0, [1], ["1 x C=C"], "RDKit agrees"],
  ["benzene", "c1ccccc1", 0, 0, [2], ["2 x C=C C=C C=C"], "RDKit agrees"],
  ["toluene", "Cc1ccccc1", 0, 0, [2], ["2 x C=C C=C C=C"], "RDKit agrees"],
  ["phenol", "Oc1ccccc1", 0, 0, [2, 2, 1], ["3 x C-1 O+1 C=C C=C C=O","2 x C=C C=C C=C"], "no charge filter: O(+)H=C with C(-) at ortho and para kept"],
  ["aniline", "Nc1ccccc1", 0, 0, [2, 2, 1], ["3 x C-1 N+1 C=C C=C C=N","2 x C=C C=C C=C"], "no charge filter: N(+)H2=C with C(-) at ortho and para kept"],
  ["nitrobenzene", "[O-][N+](=O)c1ccccc1", 0, 0, [4], ["4 x N+1 O-1 C=C C=C C=C N=O"], "RDKit agrees"],
  ["pyridine", "c1ccncc1", 0, 0, [2], ["2 x C=C C=C C=N"], "RDKit agrees"],
  ["pyridine N-oxide", "[O-][n+]1ccccc1", 0, 0, [2, 2, 1], ["3 x C-1 N+1 C=C C=C N=O","2 x N+1 O-1 C=C C=C C=N"], "no charge filter: O=N(+) with C(-) at ortho and para kept"],
  ["pyrrole", "c1cc[nH]c1", 0, 0, [2, 2, 1], ["4 x C-1 N+1 C=C C=N","1 x C=C C=C"], "no charge filter: N(+)=C with C(-) kept"],
  ["furan", "c1ccoc1", 0, 0, [2, 2, 1], ["4 x C-1 O+1 C=C C=O","1 x C=C C=C"], "no charge filter: O(+)=C with C(-) kept"],
  ["naphthalene", "c1ccc2ccccc2c1", 0, 0, [2, 1], ["3 x C=C C=C C=C C=C C=C"], "RDKit lists one Kekule form; the three are right"],
  ["DMSO", "CS(=O)C", 0, 0, [1], ["1 x O-1 S+1"], "octet ceiling: S(+)-O(-), not S=O"],
  ["SF6", "FS(F)(F)(F)(F)F", 0, 0, [15], ["15 x F-1 F-1 S+2 2 pairs cut"], "octet ceiling: S(2+) with two F(-) cut, fifteen ways"],
  ["SF4", "FS(F)(F)F", 0, 0, [4], ["4 x F-1 S+1 1 pair cut"], "octet ceiling: S(+) with one F(-) cut"],
  ["PCl5", "ClP(Cl)(Cl)(Cl)Cl", 0, 0, [5], ["5 x Cl-1 P+1 1 pair cut"], "octet ceiling: P(+) with one Cl(-) cut"],
  ["POCl3", "O=P(Cl)(Cl)Cl", 0, 0, [1], ["1 x O-1 P+1"], "octet ceiling: P(+)-O(-)"],
  ["XeF4", "F[Xe](F)(F)F", 0, 0, [6], ["6 x F-1 F-1 Xe+2 2 pairs cut"], "octet ceiling: Xe(2+) with two F(-) cut"],
  ["XeF2", "F[Xe]F", 0, 0, [2], ["2 x F-1 Xe+1 1 pair cut"], "octet ceiling: Xe(+) with one F(-) cut"],
  ["triiodide (Na+)", "I[I-]I.[Na+]", 0, 0, [2], ["2 x I-1 Na+1 1 pair cut"], "octet ceiling: I-I with I(-) cut, two ways"],
  ["borane", "B", 2, 0, [1], ["1 x single bonds only"], "RDKit agrees"],
  ["boron trifluoride", "FB(F)F", 0, 0, [3], ["3 x B-1 F+1 B=F"], "fewest holes: F(+)=B(-) forms only; the picture with six electrons on B is dropped (his call: see REPORT)"],
  ["ammonia borane", "[NH3+][BH3-]", 0, 0, [1], ["1 x B-1 N+1"], "RDKit agrees"],
  ["sulfuric acid", "OS(=O)(=O)O", 0, 0, [1], ["1 x O-1 O-1 S+2"], "octet ceiling: S(2+) with two O(-)"],
  ["nitric acid", "O[N+](=O)[O-]", 0, 0, [2, 1], ["1 x N+1 O+1 O-1 O-1 N=O","2 x N+1 O-1 N=O"], "no charge filter: O(-)N(+)(O(-))=O(+)H kept beside the two textbook forms"],
  ["nitrogen dioxide (radical)", "[O]N=O", 1, 1, [2, 2, 2], ["2 x N+1 N. O-1 N=O","2 x O+1 O-1 O. N=O","2 x O. N=O"], "no charge filter: radical on N (N(+) O(-)) and on O, both ways"],
  ["nitric oxide (radical)", "[N]=O", 1, 1, [1, 1], ["1 x N-1 O+1 O. N=O","1 x N. N=O"], "no charge filter: N(-)=O(+) radical kept beside N=O"],
  ["methyl radical", "[CH3]", 1, 1, [1], ["1 x C."], "RDKit agrees"],
  ["acetate", "CC(=O)[O-]", 0, 0, [2], ["2 x O-1 C=O"], "RDKit agrees"],
  ["formate", "C(=O)[O-]", 0, 0, [2], ["2 x O-1 C=O"], "RDKit agrees"],
  ["enolate of acetaldehyde", "C=C[O-]", 0, 0, [1, 1], ["1 x C-1 C=O","1 x O-1 C=C"], "RDKit agrees"],
  ["enolate of acetone", "CC(=C)[O-]", 0, 0, [1, 1], ["1 x C-1 C=O","1 x O-1 C=C"], "RDKit agrees"],
  ["allyl cation", "C=C[CH2+]", 2, 0, [2], ["2 x C+1 C=C"], "RDKit agrees"],
  ["allyl anion", "C=C[CH2-]", 0, 0, [2], ["2 x C-1 C=C"], "RDKit agrees"],
  ["allyl radical", "C=C[CH2]", 1, 1, [2], ["2 x C. C=C"], "RDKit agrees"],
  ["benzyl cation", "[CH2+]c1ccccc1", 2, 0, [2, 2, 1], ["5 x C+1 C=C C=C C=C"], "RDKit agrees"],
  ["benzyl anion", "[CH2-]c1ccccc1", 0, 0, [2, 2, 1], ["5 x C-1 C=C C=C C=C"], "RDKit agrees"],
  ["benzyl radical", "[CH2]c1ccccc1", 1, 1, [2, 2, 1], ["5 x C. C=C C=C C=C"], "RDKit lists one form; the ring forms are textbook"],
  ["nitrate", "[O-][N+](=O)[O-]", 0, 0, [3], ["3 x N+1 O-1 O-1 N=O"], "RDKit agrees"],
  ["nitrite", "[O-]N=O", 0, 0, [2], ["2 x O-1 N=O"], "fewest holes: RDKit adds O(-)-N(+)-O(-) with a six-electron N"],
  ["methanesulfonate", "CS(=O)(=O)[O-]", 0, 0, [1], ["1 x O-1 O-1 O-1 S+2"], "octet ceiling: S(2+) with three O(-)"],
  ["phosphate", "[O-]P(=O)([O-])[O-]", 0, 0, [1], ["1 x O-1 O-1 O-1 O-1 P+1"], "octet ceiling: P(+) with four O(-)"],
  ["azide", "[N-]=[N+]=[N-]", 0, 0, [1], ["1 x N+1 N-1 N-1 N=N N=N"], "magnitude: N#N(+)-N(2-) dropped; fewest holes drops RDKit's radical and six-electron forms"],
  ["hydrazoic acid", "N=[N+]=[N-]", 0, 0, [1, 1], ["1 x N+1 N-1 N#N","1 x N+1 N-1 N=N N=N"], "magnitude and fewest holes, as for azide"],
  ["sulfur trioxide", "O=S(=O)=O", 0, 0, [3], ["3 x O-1 O-1 S+2 O=S"], "octet ceiling: S(2+) with one double bond"],
  ["perchlorate", "[O-]Cl(=O)(=O)=O", 0, 0, [1], ["1 x Cl+3 O-1 O-1 O-1 O-1"], "octet ceiling: Cl(3+) with four O(-); RDKit happens to give the same"],
  ["pyridinium", "c1cc[nH+]cc1", 0, 0, [2], ["2 x N+1 C=C C=C C=N"], "RDKit agrees"],
  ["thiophene", "c1ccsc1", 0, 0, [2, 2, 1], ["4 x C-1 S+1 C=C C=S","1 x C=C C=C"], "no charge filter: S(+)=C with C(-) kept"],
  ["imidazole", "c1cnc[nH]1", 0, 0, [1, 1, 1, 1, 1], ["3 x C-1 N+1 C=N C=N","1 x C=C C=N","1 x N+1 N-1 C=C C=N"], "no charge filter: N(+)=C with C(-) or N(-) kept"],
  ["cyclopentadienyl anion", "[CH-]1C=CC=C1", 0, 0, [5], ["5 x C-1 C=C C=C"], "RDKit agrees"],
  ["tropylium", "[CH+]1C=CC=CC=C1", 2, 0, [7], ["7 x C+1 C=C C=C C=C"], "RDKit agrees"],
  ["guanidinium", "NC(N)=[NH2+]", 0, 0, [3], ["3 x N+1 C=N"], "RDKit agrees"],
  ["formamidinium", "NC=[NH2+]", 0, 0, [2], ["2 x N+1 C=N"], "RDKit agrees"],
  ["acetylacetonate", "CC(=O)[CH-]C(C)=O", 0, 0, [2, 1], ["1 x C-1 C=O C=O","2 x O-1 C=C C=O"], "RDKit agrees"],
  ["nitrosonium", "[O+]#N", 0, 0, [1], ["1 x O+1 N#O"], "RDKit agrees"],
  ["nitronium", "O=[N+]=O", 0, 0, [2, 1], ["1 x N+1 N=O N=O","2 x N+1 O+1 O-1 N#O"], "no charge filter: O(-)-N(+)#O(+) kept, as for CO2"],
  ["cyanate", "[O-]C#N", 0, 0, [1, 1], ["1 x N-1 C=N C=O","1 x O-1 C#N"], "RDKit agrees"],
  ["methyl isocyanate", "CN=C=O", 0, 0, [1, 1, 1], ["1 x C=N C=O","1 x N+1 O-1 C#N","1 x N-1 O+1 C#O"], "RDKit lacks CH3-N(-)-C#O(+), a textbook form"],
  ["thiocyanate", "[S-]C#N", 0, 0, [1, 1], ["1 x N-1 C=N C=S","1 x S-1 C#N"], "RDKit agrees"],
  ["fulminate", "[C-]#[N+][O-]", 0, 0, [1], ["1 x C-1 N+1 O-1 C#N"], "RDKit agrees"],
  ["carbon suboxide", "O=C=C=C=O", 0, 0, [2, 2, 1], ["2 x C-1 O+1 C#O C=C C=O","1 x C=C C=C C=O C=O","2 x O+1 O-1 C#C C#O"], "no charge filter: O(+)#C-C#C-O(-) and O=C=C(-)-C#O(+) kept"],
  ["butatriene", "C=C=C=C", 0, 0, [1], ["1 x C=C C=C C=C"], "RDKit agrees"],
  ["vinylacetylene", "C=CC#C", 0, 0, [1], ["1 x C#C C=C"], "RDKit agrees"],
  ["carbon disulfide", "S=C=S", 0, 0, [2, 1], ["1 x C=S C=S","2 x S+1 S-1 C#S"], "no charge filter: S(+)#C-S(-) kept, as for CO2"],
  ["methyl isocyanide", "C[N+]#[C-]", 0, 0, [1], ["1 x C-1 N+1 C#N"], "RDKit agrees"],
  ["urea", "NC(N)=O", 0, 0, [2, 1], ["1 x C=O","2 x N+1 O-1 C=N"], "RDKit agrees"],
  ["carbonate", "[O-]C(=O)[O-]", 0, 0, [3], ["3 x O-1 O-1 C=O"], "RDKit agrees"],
  ["oxalate", "[O-]C(=O)C(=O)[O-]", 0, 0, [4], ["4 x O-1 O-1 C=O C=O"], "RDKit agrees"],
  ["vinylamine (enamine)", "C=CN", 0, 0, [1, 1], ["1 x C-1 N+1 C=N","1 x C=C"], "no charge filter: CH2(-)-CH=N(+)H2 kept"],
  ["chlorobenzene", "Clc1ccccc1", 0, 0, [2, 2, 1], ["3 x C-1 Cl+1 C=C C=C C=Cl","2 x C=C C=C C=C"], "no charge filter: Cl(+)=C with C(-) at ortho and para kept"],
  ["fluoroethylene", "C=CF", 0, 0, [1, 1], ["1 x C-1 F+1 C=F","1 x C=C"], "no charge filter: CH2(-)-CH=F(+) kept"],
  ["2-pyridone", "O=c1cccc[nH]1", 0, 0, [1, 1, 1, 1, 1], ["2 x C-1 N+1 C=C C=N C=O","1 x C=C C=C C=O","2 x N+1 O-1 C=C C=C C=N"], "no charge filter; RDKit also misses the second Kekule form of the O(-) N(+)H picture"],
  ["acylium", "CC#[O+]", 0, 0, [1], ["1 x O+1 C#O"], "RDKit agrees"],
  ["protonated formaldehyde", "C=[OH+]", 0, 0, [1], ["1 x O+1 C=O"], "RDKit agrees"],
  ["methyleneiminium", "C=[NH2+]", 0, 0, [1], ["1 x N+1 C=N"], "RDKit agrees"],
  ["hypochlorite", "[O-]Cl", 0, 0, [1], ["1 x O-1"], "RDKit agrees"],
  ["cyanamide", "NC#N", 0, 0, [1, 1], ["1 x C#N","1 x N+1 N-1 C=N C=N"], "RDKit agrees"],
  ["carbodiimide", "N=C=N", 0, 0, [2, 1], ["1 x C=N C=N","2 x N+1 N-1 C#N"], "RDKit agrees"],
  ["superoxide", "[O][O-]", 1, 1, [2], ["2 x O-1 O."], "RDKit agrees"],
  ["hydroxyl radical", "[OH]", 1, 1, [1], ["1 x O."], "RDKit agrees"],
  ["nitrogen dioxide", "[O]N=O", 1, 1, [2, 2, 2], ["2 x N+1 N. O-1 N=O","2 x O+1 O-1 O. N=O","2 x O. N=O"], "no charge filter: radical on N (N(+) O(-)) and on O, both ways"],
  ["nitric oxide", "[N]=O", 1, 1, [1, 1], ["1 x N-1 O+1 O. N=O","1 x N. N=O"], "no charge filter: N(-)=O(+) radical kept beside N=O"],
];

// the same index-free picture of a structure the table uses
export function picture(atoms, st, edges, fixed) {
  const marks = [];
  st.charge.forEach((q, i) => { if (q) marks.push(`${atoms[i]}${q > 0 ? '+' : ''}${q}`); });
  st.lone.forEach((l, i) => { if (l % 2) marks.push(`${atoms[i]}.`); });
  const bonds = [];
  edges.forEach(([i, j], k) => { if (st.orders[k] > 1) bonds.push([atoms[i], atoms[j]].sort().join(st.orders[k] === 2 ? '=' : '#')); });
  const cuts = edges.reduce((s, [i, j], k) => s + (st.orders[k] === 0 && !fixed[i] && !fixed[j] ? 1 : 0), 0);
  return [...marks.sort(), ...bonds.sort(), ...(cuts ? [`${cuts} pair${cuts > 1 ? 's' : ''} cut`] : [])].join(' ') || 'single bonds only';
}

export function picturesOf(node, r) {
  const atoms = node.atoms.map(a => a.el);
  const tally = new Map();
  r.structures.forEach(st => { const k = picture(atoms, st, r.edges, node.fixed); tally.set(k, (tally.get(k) || 0) + 1); });
  return [...tally.entries()].sort().map(([k, c]) => `${c} x ${k}`);
}

if (process.argv[1] && process.argv[1].endsWith('rdkit-agreed.js')) {
  let bad = 0;
  for (const [name, smiles, holes, radicals, sizes, pics, note] of AGREED) {
    const problems = [];
    try {
      const node = buildNode(nodeFromSmiles(smiles));
      const r = enumerate(node);
      const fam = structureFamilies(node, r);
      const got = fam.families.map(f => f.length).sort((a, b) => b - a);
      if (r.holes !== holes) problems.push(`holes ${r.holes}, expected ${holes}`);
      if (r.radicals !== radicals) problems.push(`radicals ${r.radicals}, expected ${radicals}`);
      if (JSON.stringify(got) !== JSON.stringify(sizes)) problems.push(`families [${got}], expected [${sizes}]`);
      const gotPics = picturesOf(node, r);
      if (JSON.stringify(gotPics) !== JSON.stringify(pics)) problems.push(`pictures ${JSON.stringify(gotPics)}, expected ${JSON.stringify(pics)}`);
    } catch (e) { problems.push(e.message); }
    if (problems.length) bad++;
    console.log(`${problems.length ? 'FAIL' : 'ok  '} ${name}: ${note}`);
    problems.forEach(p => console.log('  PROBLEM: ' + p));
  }
  console.log(bad ? `${bad} problem(s)` : `all ${AGREED.length} recorded answers hold`);
  process.exit(bad ? 1 : 0);
}
