Electron Flow

Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved; see LICENSE.

A rules-based 3D tool for building and viewing molecules, resonance and electron movement. Design decisions live in the project doc 00 DESIGN; this file only says how the code is laid out.

core/      chemistry, no dependence on any renderer
  elements.js   valence electrons, octet ceiling, radii (covalent and ionic), colours
  lewis.js      Lewis structure enumeration, filters, equal weights, averages, transition hybrids
  geometry.js   electron groups, template shapes, placement, spring relaxation, where loose ions sit
  cases.js      known-answer cases and runCase (transition hybrids live here)
  symmetry.js   symmetries of the skeleton, families of equivalent structures
  pi.js         where pi bonds and lone pairs point
  vectors.js    the vector arithmetic shared by geometry.js and pi.js
  smiles.js     molecules from SMILES (openchemlib reads the text; the engine redoes the chemistry) and back
  edit.js       editing operations on a molecule (add, remove, bond, charge, fill with hydrogens)
  compounds.js  84 named compounds as SMILES with counter-ions: a test list (test/compounds.js, tools/all.py)
  names.js      the name search: comparing names, searching the index, counter-ions for a bare ion, the
                skeleton key that names what is on screen
  moves.js      curved arrows applied to a structure: resonance, a move to a new node, or a refusal with the
                reason; the node after a move; the transition hybrid at a run of mixes; the shapes of the film
  reactions.js  reaction scenes (reagents placed together) with their textbook arrows, and piece formulas
test/      node test/run.js (known answers), families.js, pi.js, compounds.js (timings), geometry.js, extra.js,
           gradient.js (exact forces against finite differences), shapes.js (known-answer shapes), smiles.js
           (reading and writing SMILES), edit.js (editing operations), rdkit-agreed.js (138 recorded answers;
           rdkit-dump.js and rdkit-compare.py made them), rules-lab.js (candidate rules, prints only),
           moves.js (arrows on the basic mechanisms, refusals, hybrids), reactions.js (every scene's arrows
           give the textbook products), rings.js (the ring rule: benzyne, pyridine, pyridazine, the shapes),
           film.js (the film's shapes: no swing, no jump), names.js (the name search on the shipped index)
web/       index.html, app.js (three.js viewer, the search box, the film of a move), edit.js (editing on the page),
           arrows.js (drawing curved arrows with the pointer), names.json.gz (the names, built by
           tools/build-names.mjs), vendor/ (three.js and openchemlib shipped with the page)
build.sh   assembles dist/ (web + core) for publishing
tools/     build-names.mjs builds web/names.json.gz from the chemicals Python package's PubChem lists (how is at
           its top). Driving the page in headless Chrome with playwright, for screenshots and checks: shot.py (one
           name or SMILES, one picture), drive.py (editing: clicks, keys, drags), drive-moves.py (arrows, moves, the film),
           all.py (the test list and every reaction scene on contact sheets); each says how to run it at the top. Serve dist/ on port
           8765 first: (cd dist && python3 -m http.server 8765). build.sh deletes dist/, so restart the server
           after every build or it serves nothing. reviews/ holds the briefs given to the reviewing agents.

Run tests from this folder: node test/run.js
