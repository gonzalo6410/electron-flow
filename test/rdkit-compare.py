#!/usr/bin/env python3
# Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
# python3 test/rdkit-compare.py [names...]   (run from the electron-flow folder; needs rdkit:
# pip install rdkit --break-system-packages)
#
# Holds the engine's contributing structures against RDKit's ResonanceMolSupplier, for every named
# compound and the textbook cases listed in test/rdkit-dump.js. RDKit is run three ways: its default
# flags, with ALLOW_CHARGE_SEPARATION and ALLOW_INCOMPLETE_OCTETS, and with KEKULE_ALL added so that
# aromatic rings give every Kekule structure. Every structure, ours and theirs, is turned into one
# canonical SMILES (RDKit's, with explicit charges and radicals) so sets can be compared without any
# atom mapping. Our rules differ from RDKit's on purpose, so a difference is not a bug by itself; the
# script only prints the differences, and test/rdkit-agreed.js records the judged outcome.
import json, subprocess, sys, os
from rdkit import Chem, RDLogger
from rdkit.Chem import rdchem

RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ELEMENT_VALENCE = {'H': 1, 'He': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8, 'Al': 3, 'Si': 4,
                   'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8, 'Se': 6, 'Br': 7, 'Kr': 8, 'I': 7, 'Xe': 8,
                   'Li': 1, 'Na': 1, 'K': 1, 'Be': 2, 'Mg': 2, 'Ca': 2}


def engine_dump():
    out = subprocess.run(['node', os.path.join(HERE, 'rdkit-dump.js')], capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0:
        sys.exit(out.stderr)
    return json.loads(out.stdout)


def smiles_of_ours(entry, st):
    # build an RDKit molecule from one of our structures and return its canonical Kekule SMILES
    m = Chem.RWMol()
    for el in entry['atoms']:
        m.AddAtom(Chem.Atom(el))
    for i, j, o in st['bonds']:
        m.AddBond(i, j, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[o])
    for i, a in enumerate(m.GetAtoms()):
        a.SetFormalCharge(int(st['charge'][i]))
        a.SetNumRadicalElectrons(int(st['lone'][i]) % 2)
        a.SetNoImplicit(True)
    return canonical(m)


def canonical(m):
    # one SMILES per picture: bond orders as given, no aromaticity, so the two Kekule forms of a ring
    # with a substituent stay two pictures (they are two pictures for us too) while the two of
    # benzene, which are the same picture with atoms renamed, come out identical
    try:
        m = Chem.RWMol(m)
        m.UpdatePropertyCache(strict=False)
        if any(b.GetIsAromatic() for b in m.GetBonds()):
            Chem.Kekulize(m, clearAromaticFlags=True)
        for a in m.GetAtoms():
            a.SetIsAromatic(False)
        ops = Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_KEKULIZE
        Chem.SanitizeMol(m, sanitizeOps=ops, catchErrors=True)
        m = Chem.RemoveHs(m, sanitize=False)
        return Chem.MolToSmiles(m, kekuleSmiles=True, canonical=True)
    except Exception as e:  # keep going: an unsanitizable picture still gets a label
        return 'unsanitizable:' + Chem.MolToSmiles(m, canonical=False)


def rdkit_structures(smiles, flags):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    hs = [a.GetTotalNumHs() for a in m.GetAtoms()]
    sup = rdchem.ResonanceMolSupplier(m, flags)
    out = []
    for r in sup:
        if r is None:
            continue
        r = Chem.RWMol(r)
        # keep the hydrogens the input had: a structure with fewer bonds must not grow new ones
        for a, h in zip(r.GetAtoms(), hs):
            a.SetNumExplicitHs(h)
            a.SetNoImplicit(True)
        out.append(canonical(r))
    return sorted(set(out))


def main():
    wanted = sys.argv[1:]
    data = engine_dump()
    FULL = Chem.ALLOW_CHARGE_SEPARATION | Chem.ALLOW_INCOMPLETE_OCTETS
    summary = []
    for entry in data:
        if wanted and not any(w in entry['name'] for w in wanted):
            continue
        print('\n== %s   %s' % (entry['name'], entry['smiles']))
        if 'error' in entry:
            print('  engine error: ' + entry['error'])
            continue
        ours = sorted(set(smiles_of_ours(entry, st) for st in entry['structures']))
        rd_default = rdkit_structures(entry['smiles'], 0)
        rd_full = rdkit_structures(entry['smiles'], FULL)
        rd_kek = rdkit_structures(entry['smiles'], FULL | Chem.KEKULE_ALL)
        if rd_default is None:
            print('  rdkit could not read the SMILES')
            continue
        print('  engine: %d structure(s) in %d famil%s, %d holes, %s' % (
            len(entry['structures']), len(entry['families']), 'y' if len(entry['families']) == 1 else 'ies',
            entry['holes'], 'radical' if entry['radicals'] else 'closed shell'))
        for s in ours:
            print('    ours   ' + s)
        print('  rdkit: default %d, charge separation + incomplete octets %d, + every Kekule form %d' % (
            len(rd_default), len(rd_full), len(rd_kek)))
        for s in rd_kek:
            print('    rdkit  ' + s)
        only_ours = [s for s in ours if s not in rd_kek]
        only_rd = [s for s in rd_kek if s not in ours]
        for s in only_ours:
            print('  ONLY OURS  ' + s)
        for s in only_rd:
            print('  ONLY RDKIT ' + s)
        summary.append((entry['name'], len(ours), len(rd_kek), len(only_ours), len(only_rd)))
    print('\n== summary (name, ours, rdkit, only ours, only rdkit)')
    for row in summary:
        print('  %-45s %3d %3d %3d %3d%s' % (row[0], row[1], row[2], row[3], row[4], '' if row[3] == 0 and row[4] == 0 else '  DIFFERS'))


if __name__ == '__main__':
    main()
