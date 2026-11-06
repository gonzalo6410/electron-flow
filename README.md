Electron Flow

A rules-based 3D tool for building and viewing molecules, resonance and electron movement. Design decisions live in the project doc 00 DESIGN; this file only says how the code is laid out.

core/      chemistry, no dependence on any renderer
  elements.js   valence electrons, octet ceiling, radii, colours
  lewis.js      Lewis structure enumeration, filters, equal weights, averages, transition hybrids
  geometry.js   electron groups, template shapes, placement, spring relaxation, lone pair directions
  cases.js      known-answer cases and runCase (transition hybrids live here)
  symmetry.js   symmetries of the skeleton, families of equivalent structures
  pi.js         where pi bonds and lone pairs point
  smiles.js     molecules from SMILES (openchemlib reads the text; the engine redoes the chemistry)
  compounds.js  named compounds for the picker, as SMILES with counter-ions
test/      node test/run.js (known answers), families.js, pi.js, compounds.js, geometry.js, extra.js, rules-lab.js
web/       index.html, app.js (three.js viewer), vendor/ (three.js and openchemlib shipped with the page)
build.sh   assembles dist/ (web + core) for publishing

Run tests from this folder: node test/run.js
