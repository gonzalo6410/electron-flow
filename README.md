Electron Flow

A rules-based 3D tool for building and viewing molecules, resonance and electron movement. Design decisions live in the project doc 00 DESIGN; this file only says how the code is laid out.

core/      chemistry, no dependence on any renderer
  elements.js   valence electrons, octet ceiling, radii, colours
  lewis.js      Lewis structure enumeration, filters, equal weights, averages, transition hybrids
  geometry.js   electron groups, template shapes, placement, spring relaxation, lone pair directions
  cases.js      known-answer cases and runCase
test/      node test/run.js (known answers), node test/geometry.js, node test/extra.js
web/       index.html, app.js (three.js viewer), vendor/ (three.js shipped with the page)
build.sh   assembles dist/ (web + core) for publishing

Run tests from this folder: node test/run.js
