// node test/names.js   the name search (core/names.js) on the shipped index (web/names.json.gz):
// what a few typed names find first, the counter-ions a bare ion gets, the skeleton key (the same for
// every way of writing a molecule and for every structure of it), and naming the products of a move.
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { nodeFromSmiles } from '../core/smiles.js';
import { buildNode, enumerate } from '../core/lewis.js';
import { makeIndex, search, namesFor, balanceSmiles, skeletonKey, normName } from '../core/names.js';
import { commit } from '../core/moves.js';
import { REACTIONS, resolveArrows } from '../core/reactions.js';

let bad = 0;
const say = (ok, text) => { if (!ok) bad++; console.log(`${ok ? 'ok  ' : 'FAIL'} ${text}`); };

const t0 = performance.now();
const data = JSON.parse(gunzipSync(readFileSync(new URL('../web/names.json.gz', import.meta.url))).toString());
const index = makeIndex(data);
const shown = index.entries.filter(e => e.smiles).length;
say(shown > 20000, `index: ${index.entries.length} entries, ${shown} that can be shown, ${index.flat.length} names, read in ${(performance.now() - t0).toFixed(0)} ms`);

// what comes first for a typed name: [typed, the SMILES the first hit must have (same skeleton)]
const first = [
  ['benzene', 'c1ccccc1'],
  ['Benzene', 'c1ccccc1'],
  ['acetic acid', 'CC(=O)O'],
  ['ethanoic acid', 'CC(=O)O'],
  ['acetone', 'CC(C)=O'],
  ['propan-2-one', 'CC(C)=O'],
  ['tert-butanol', 'CC(C)(C)O'],
  ['dmso', 'CS(C)=O'],
  ['caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'],
  ['benzyne', 'C1=CC#CC=C1'],
  ['sulfur hexafluoride', 'FS(F)(F)(F)(F)F'],
  ['acetate', 'CC(=O)[O-].[Na+]'],
  ['sodium acetate', 'CC(=O)[O-].[Na+]'],
  ['water', 'O'],
  ['THF', 'C1CCOC1'],
  ['LDA', 'CC(C)[N-]C(C)C.[Li+]'],
  ['sodium hydroxide', '[OH-].[Na+]'],
  ['butanol', 'CCCCO'],
];
const keyOf = smi => skeletonKey(buildNode(nodeFromSmiles(balanceSmiles(smi).smiles)));
for (const [typed, smi] of first) {
  const t = performance.now();
  const hits = search(index, typed, 8);
  const ms = performance.now() - t;
  const top = hits[0];
  const ok = top && top.entry.smiles && keyOf(top.entry.smiles) === keyOf(smi);
  say(ok, `"${typed}" finds ${top ? top.entry.name : 'nothing'} first (${ms.toFixed(1)} ms)`);
}
{
  const t = performance.now();
  search(index, 'chloro', 8);
  const ms = performance.now() - t;
  say(ms < 200, `a short common fragment ("chloro") searches the whole list in ${ms.toFixed(0)} ms`);
}
{
  // the name shown: a salt by its own name, water as water (not as heavy water, which the engine,
  // ignoring isotopes, cannot tell apart from it)
  for (const [typed, want] of [['sodium hydroxide', 'sodium hydroxide'], ['water', 'water'], ['sodium chloride', 'sodium chloride']]) {
    const top = search(index, typed, 1)[0];
    say(top && top.entry.name === want, `"${typed}" is shown as ${top ? top.entry.name : 'nothing'}`);
  }
  const moth = search(index, 'thf', 8).some(h => /moth/.test(h.matched));
  say(!moth, '"thf" does not find moth flakes (inside a name, only from the start of a word)');
}
{
  const hits = search(index, 'ferrocene', 8);
  const h = hits.find(x => x.entry.name === 'ferrocene');
  say(!!h && !h.entry.smiles && /^el:/.test(h.entry.key), `ferrocene is listed with the reason it cannot be shown: ${h ? h.entry.key : 'not found'}`);
}

// counter-ions
{
  const a = balanceSmiles('CC(=O)[O-]');
  say(a.smiles === 'CC(=O)[O-].[Na+]' && a.added.ion === 'Na+', `a bare acetate gets a sodium: ${a.smiles}`);
  const b = balanceSmiles('[O-]S(=O)(=O)[O-]');
  say(b.added.count === 2, `a bare sulfate gets two: ${b.smiles}`);
  const c = balanceSmiles('C[N+](C)(C)C');
  say(c.added.ion === 'Cl-', `a bare tetramethylammonium gets a chloride: ${c.smiles}`);
  const d = balanceSmiles('CC(=O)[O-].[Na+]');
  say(d.added === null, 'a balanced salt is left alone');
}

// the skeleton key
{
  const same = [['CCO', 'OCC', 'C(O)C'], ['c1ccccc1', 'C1=CC=CC=C1'], ['C1=CC#CC=C1', 'C1=CC=C=C=C1'], ['CC(=O)[O-].[Na+]', '[Na+].[O-]C(C)=O']];
  for (const group of same) {
    const keys = group.map(s => skeletonKey(buildNode(nodeFromSmiles(s))));
    say(keys.every(k => k === keys[0]), `one key for ${group.join(', ')}`);
  }
  const k1 = skeletonKey(buildNode(nodeFromSmiles('[CH3+].[Cl-]'))), k2 = skeletonKey(buildNode(nodeFromSmiles('[CH3].[Cl]')));
  say(k1 !== k2, 'a methyl cation beside chloride is not a methyl radical beside a chlorine atom');
  say(normName('2-Methylpropan-2-ol') === '2methylpropan2ol', 'names compare without case, spaces or hyphens');
}

// naming what is on screen after a move: SN2 gives methanol
{
  const rx = REACTIONS.find(x => x.name.startsWith('SN2'));
  const spec = nodeFromSmiles(rx.smiles), node = buildNode(spec), r = enumerate(node);
  const out = commit(node, r.structures[0], resolveArrows(spec, rx.arrows), r);
  // the products as one scene: methanol, and the sodium tied to the chloride
  const got = out.kind === 'move' ? namesFor(index, out.spec) : null;
  say(!!got && got.join(' + ') === 'methanol + sodium chloride', `the SN2 products are named: ${got ? got.join(' + ') : 'nothing'}`);
}

if (bad) { console.log(`${bad} failed`); process.exit(1); }
console.log('name checks pass');
