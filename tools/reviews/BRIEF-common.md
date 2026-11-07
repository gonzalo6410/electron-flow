(Paths below are the workspace layout used on 6 Nov 2026: a snapshot of the tree as ORIGINAL, the reviewer's own copy beside it, everything under /home/claude/agents. Recreate that layout or edit the paths.)

ELECTRON FLOW: common brief for reviewers

What this is. A rules-based 3D teaching tool for Lewis structures and resonance, built for a chemistry PhD student (Parkin group, Columbia). He designs, we code. The code is plain JavaScript ES modules, no framework, no build step beyond copying files. Chemistry core in core/, tests in test/ (run every file in test/ with node from the folder root, except rules-lab.js and rdkit-dump.js), viewer in web/ (three.js and openchemlib shipped in web/vendor). npm install has already been done in your copy.

The rules of the chemistry (his decisions; do not change them, flag anything that seems wrong):
- Lewis structures: the octet is a ceiling (never exceeded, second and third row alike; no expanded octets; hypervalent species are drawn with charge separation, SF6 as S2+ with four bonds and two F-) and a target (contributing structures are the ones with the fewest octet holes; holes only when nothing better exists, so carbocations and boranes are legal).
- Sigma bonds do not break inside a node: a pair of adjacent atoms loses its bond only where the octet ceiling leaves no room for it (SF6's two spare F-). So POCl3 is one structure, P+ with four single bonds and O-.
- No formal charge beyond +1 or -1 unless every structure has one.
- Ring rule (7 Nov 2026): in a ring of seven or fewer, two neighbouring ring atoms may share one bond in the ring's plane (benzyne's third bond); an atom with two pi parts and no such partner is refused (pyridine's C=N(+)=C); electrons in the ring's plane stay put across structures (pyridazine keeps its N lone pairs).
- All surviving structures weigh the same. No charge filter, no symmetry filter: charge-separated forms stay (formamide has O(-)-C=N(+) at 1/2; CO2 has O(-)-C#O(+) at 1/3 each).
- Group 1 and 2 metals and single-atom anions written as separate pieces are spectator ions: no bonds, fixed charge.
- Geometry is computed from the Lewis structure by rules (electron groups give the shape, springs relax lengths and angles), nothing stored per molecule.
- Structures equivalent under a relabelling of atoms that preserves elements and adjacency form a family.

Working rules for you:
- Work only in your own copy (the folder you are given). Do not touch /home/claude/ef-repo or /home/claude/electron-flow.
- Keep the code style: plain prose comments in full sentences, no jargon, short functions, no new dependencies unless they are shipped as a file in web/vendor and justified. Do not add TypeScript, bundlers or frameworks.
- Every change must keep all existing tests passing, and new behaviour needs a test with a known answer.
- Network: npm and pip work (pip needs --break-system-packages), GitHub clones work, other sites usually do not. RDKit is available with pip install rdkit; 3dmol and ngl are on npm.
- When you are done, write two things in your folder: PATCH.diff (a unified diff of your copy against the original, made with: diff -ruN /home/claude/agents/ORIGINAL . -x node_modules -x dist -x PATCH.diff -x REPORT.md > PATCH.diff) and REPORT.md (facts first, then what you changed, then what you did not change and why; tag each finding confirmed, candidate or provisional; keep it short and in plain prose).
- Never invent chemistry or numbers. If unsure, say so in the report.

The person and the screen (this applies to every reviewer, not only the viewer one). Everything ends up as what a person sees and does on the page, so judge your area by its effect on the page as well as by the code. Open the page yourself: ./build.sh then (cd dist && python3 -m http.server 8765) and a headless Chrome screenshot through playwright (python3 with playwright is installed; the browser is /usr/bin/google-chrome; launch with args --use-gl=swiftshader --enable-unsafe-swiftshader --no-sandbox; tools/shot.py is an example: python3 -I tools/shot.py "NAME OR SMILES" on|off WAIT_MS OUT.png). Look at several molecules, with the resonance switch on and off, after rotating. His taste: plain and simple, black text, one font, pictures over words, the important things big and the rest left out, nothing extravagant; alternation between structures slow and never snappy; nothing on screen that does not earn its place; no jargon anywhere in the interface. In REPORT.md include a section "What a person sees" with what you noticed and what you changed about it, and mention anything clunky you saw outside your area so it can be passed on.
Note on shot.py: it waits for the page, types the name or SMILES into the search box, switches resonance on when asked, waits WAIT_MS and saves a screenshot at double resolution; it prints console errors and the status line. Pick a free port if 8765 is taken (edit the script).
