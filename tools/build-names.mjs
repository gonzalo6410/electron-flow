// Build the name search index (web/names.json.gz) from the compound lists of the `chemicals` Python
// package (MIT; Caleb Bell), whose names, synonyms and SMILES come from PubChem (public domain).
//
//   pip download chemicals --no-deps -d /tmp/chem && (cd /tmp/chem && python3 -m zipfile -e chemicals-*.whl x)
//   node tools/build-names.mjs /tmp/chem/x/chemicals/Identifiers
//   node tools/build-names.mjs --curated     (only redo our own entries, on the file already built)
//   node tools/build-names.mjs --add <folder> <list file>...   (run and merge in only what is new)
// A full build takes about fifty minutes on two cores.
//
// Our own entries come first: the ions, cations and radicals of core/compounds.js and a few more, under
// the plain names a student types ("acetate", "methyl cation"). PubChem's lists have few bare ions,
// and its synonyms point "acetate" at acetic acid.
// Every compound is run through the engine as the page would run it (read, counter-ions added to a
// bare ion, structures found, families, shape). Compounds the engine finishes within the time limit go
// in with their SMILES; the rest go in with their main names and the reason, so a search for them says
// why instead of saying nothing (an element not in the table yet, too slow for the engine yet).
// Compounds that are the same molecule to the engine (stereoisomers, which it ignores) are merged into
// one entry, and each entry carries a hash of the molecule's skeleton so the page can name what is on
// screen after an edit or a move (core/names.js).

import { isMainThread, Worker, parentPort } from 'node:worker_threads';
import { readFileSync, writeFileSync } from 'node:fs';
import { gzipSync, gunzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const LIMIT_MS = 1500;      // engine time allowed per compound (the page is slower than node)
const MAX_HEAVY = 70;       // heavier molecules are left out of the search altogether

// abbreviation, then the names it may be listed under (the first one found is used)
const ALIASES = [
  ['THF', 'tetrahydrofuran'], ['DCM', 'dichloromethane'], ['DMF', 'N,N-dimethylformamide', 'dimethylformamide'],
  ['DMA', 'N,N-dimethylacetamide', 'dimethylacetamide'], ['DMSO', 'dimethyl sulfoxide'], ['LDA', 'lithium diisopropylamide'],
  ['ATP', 'adenosine triphosphate', "adenosine 5'-triphosphate"], ['ADP', 'adenosine diphosphate', "adenosine 5'-diphosphate"],
  ['TFA', 'trifluoroacetic acid'], ['DME', '1,2-dimethoxyethane'], ['HMPA', 'hexamethylphosphoramide', 'hexamethylphosphoric triamide'],
  ['DMAP', '4-dimethylaminopyridine', '4-(dimethylamino)pyridine'], ['DBU', '1,8-diazabicyclo[5.4.0]undec-7-ene'],
  ['NBS', 'N-bromosuccinimide'], ['NCS', 'N-chlorosuccinimide'], ['mCPBA', '3-chloroperoxybenzoic acid', 'm-chloroperoxybenzoic acid', '3-chloroperbenzoic acid'],
  ['TsOH', 'p-toluenesulfonic acid', '4-methylbenzenesulfonic acid'], ['TsCl', 'p-toluenesulfonyl chloride', 'tosyl chloride'],
  ['MsCl', 'methanesulfonyl chloride'], ['Ac2O', 'acetic anhydride'], ['AcOH', 'acetic acid'], ['AcCl', 'acetyl chloride'],
  ['MeOH', 'methanol'], ['EtOH', 'ethanol'], ['iPrOH', 'isopropanol', 'isopropyl alcohol'], ['tBuOH', 'tert-butanol', 'tert-butyl alcohol'],
  ['EtOAc', 'ethyl acetate'], ['Et2O', 'diethyl ether'], ['Et3N', 'triethylamine'], ['TEA', 'triethylamine'],
  ['DIPEA', 'N,N-diisopropylethylamine', 'diisopropylethylamine'], ['MeCN', 'acetonitrile'], ['NaOEt', 'sodium ethoxide'],
  ['NaOMe', 'sodium methoxide'], ['KOtBu', 'potassium tert-butoxide'], ['NaBH4', 'sodium borohydride'],
  ['LiAlH4', 'lithium aluminum hydride', 'lithium aluminium hydride'], ['LAH', 'lithium aluminum hydride', 'lithium aluminium hydride'],
  ['NaH', 'sodium hydride'], ['BuLi', 'n-butyllithium', 'butyllithium'], ['nBuLi', 'n-butyllithium', 'butyllithium'],
  ['TMS', 'tetramethylsilane'], ['TMSCl', 'chlorotrimethylsilane', 'trimethylsilyl chloride'], ['TBAF', 'tetrabutylammonium fluoride'],
  ['DIBAL', 'diisobutylaluminum hydride', 'diisobutylaluminium hydride'], ['DCC', "N,N'-dicyclohexylcarbodiimide", 'dicyclohexylcarbodiimide'],
  ['EDTA', 'ethylenediaminetetraacetic acid', 'edetic acid'], ['Py', 'pyridine'],
];

if (!isMainThread) {
  const { nodeFromSmiles } = await import('../core/smiles.js');
  const { buildNode, enumerate } = await import('../core/lewis.js');
  const { structureFamilies } = await import('../core/symmetry.js');
  const { relax } = await import('../core/geometry.js');
  const { skeletonKey, keyHash, balanceSmiles } = await import('../core/names.js');
  parentPort.on('message', ({ id, smiles }) => {
    const t0 = performance.now();
    try {
      let spec;
      const balanced = balanceSmiles(smiles);
      const bare = !!balanced.added;
      try { spec = nodeFromSmiles(balanced.smiles); } catch (e) {
        const m = /element (\w+) is not in the table/.exec(e.message);
        parentPort.postMessage({ id, fail: m ? `el:${m[1]}` : 'read' });
        return;
      }
      const heavy = spec.atoms.filter(a => a.el !== 'H').length;
      if (heavy > MAX_HEAVY) { parentPort.postMessage({ id, fail: 'big' }); return; }
      const node = buildNode(spec);
      const r = enumerate(node);
      structureFamilies(node, r);
      relax(node, r);
      parentPort.postMessage({ id, ok: true, key: keyHash(skeletonKey(node)), bare, n: r.structures.length, ms: performance.now() - t0 });
    } catch (e) {
      parentPort.postMessage({ id, fail: /ring/.test(e.message) ? 'ring' : 'engine', msg: e.message });
    }
  });
} else {
  const here = fileURLToPath(import.meta.url);
  const out = path.join(path.dirname(here), '..', 'web', 'names.json.gz');
  const args = process.argv.slice(2);
  // the lists, most common first; a compound in two lists is read from the first
  const FILES = [
    ['chemical identifiers example user db.tsv', 0],   // the commonest few thousand (water, aspirin, chloride)
    ['chemical identifiers pubchem small.tsv', 0],
    ['Anion db.tsv', 0],                               // bare anions (acetate, nitrate, tetrafluoroborate)
    ['Cation db.tsv', 0],
    ['Inorganic db.tsv', 0],                           // water, ammonia, hydrogen chloride (and many with metals)
    ['chemical identifiers pubchem large.tsv', 1],
  ];
  if (args[0] === '--curated') {
    const data = JSON.parse(gunzipSync(readFileSync(out)).toString());
    data.entries = [...await curated(), ...data.entries.filter(e => e[4] !== 0)];
    write(out, data.entries);
  } else if (args[0] === '--add') {
    // only the compounds of the given lists that the file does not have yet are run, then merged in
    const dir = args[1];
    const data = JSON.parse(gunzipSync(readFileSync(out)).toString());
    const have = new Set(data.entries.map(e => e[4]));
    const rows = readRows(dir, FILES.filter(([f]) => args.slice(2).includes(f))).filter(r => !have.has(r.cid));
    console.log(`${rows.length} new compounds to run`);
    const results = await runEngine(rows, here);
    write(out, merge(data.entries.filter(e => e[4] !== 0), rows, results).concat(await curated()));
  } else {
    const dir = args[0];
    if (!dir) { console.error('usage: node tools/build-names.mjs <chemicals/Identifiers folder> | --add <folder> <file>... | --curated'); process.exit(1); }
    const rows = readRows(dir, FILES);
    console.log(`${rows.length} compounds read`);
    const results = await runEngine(rows, here);
    write(out, merge([], rows, results).concat(await curated()));
  }
}

function readRows(dir, files) {
  const rows = new Map(); // cid -> row
  for (const [file, tier] of files) {
    for (const line of readFileSync(path.join(dir, file), 'utf8').split('\n')) {
      const f = line.split('\t');
      if (f.length < 9 || !f[4]) continue;
      const cid = Number(f[0]);
      if (!(cid > 0) || rows.has(cid)) continue;
      rows.set(cid, { cid, tier, smiles: f[4], names: cleanNames([f[8], f[7], ...f.slice(9)]) });
    }
  }
  return [...rows.values()].filter(r => r.names.length);
}

// every row through the engine on two workers, each compound given LIMIT_MS (a worker that takes
// longer is stopped and a fresh one started)
// Results are kept in tools/.names-cache.json (not in git) under a hash of the core's code, so that
// a second build with the same engine runs only what it has not seen.
async function runEngine(all, here) {
  const { createHash } = await import('node:crypto');
  const { readdirSync, existsSync } = await import('node:fs');
  const coreDir = path.join(path.dirname(here), '..', 'core');
  const engine = createHash('sha1').update(readdirSync(coreDir).sort().map(f => readFileSync(path.join(coreDir, f), 'utf8')).join('\n')).digest('hex');
  const cacheFile = path.join(path.dirname(here), '.names-cache.json');
  let cache = { engine, results: {} };
  if (existsSync(cacheFile)) { try { const c = JSON.parse(readFileSync(cacheFile, 'utf8')); if (c.engine === engine) cache = c; } catch (e) { /* start afresh */ } }
  const results = all.map(row => cache.results[row.cid] || null);
  const todo = all.map((row, q) => q).filter(q => !results[q]);
  const list = todo.map(q => all[q]);
  console.log(`${all.length - list.length} from the cache, ${list.length} to run`);
  const fresh = new Array(list.length);
  let next = 0, done = 0;
  const t0 = Date.now();
  await new Promise(resolve => {
    if (!list.length) { resolve(); return; }
    let live = 0;
    const startLane = () => {
      live++;
      const w = new Worker(here);
      let current = -1, timer = null;
      const feed = () => {
        if (next >= list.length) { w.terminate(); live--; if (!live) resolve(); return; }
        current = next++;
        timer = setTimeout(() => {
          fresh[current] = { fail: 'slow' };
          done++;
          w.terminate();
          live--;
          startLane();
        }, LIMIT_MS + 2000); // reading and loading the worker take a moment of their own
        w.postMessage({ id: current, smiles: list[current].smiles });
      };
      w.on('message', msg => {
        clearTimeout(timer);
        fresh[msg.id] = msg.ok && msg.ms > LIMIT_MS ? { fail: 'slow' } : msg;
        if (++done % 2000 === 0) console.log(`${done} done, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
        feed();
      });
      w.on('error', e => { console.error('worker error', e.message); });
      feed();
    };
    for (let q = 0; q < 2; q++) startLane();
  });
  todo.forEach((q, t) => {
    const r = fresh[t] || { fail: 'lost' };
    results[q] = r.ok ? { ok: true, key: r.key, bare: r.bare } : { fail: r.fail };
    if (fresh[t]) cache.results[all[q].cid] = results[q];
  });
  writeFileSync(cacheFile, JSON.stringify(cache));
  return results;
}

// Entries from rows and their engine results, merged into `entries` already built. Compounds that
// are one molecule to the engine share an entry; the most common list, then the lowest PubChem
// number, gives it its SMILES and the name shown. Names are capped per compound (sixteen from a
// common list, five otherwise). Compounds the engine cannot do keep two names and the reason.
function merge(entries, rows, results) {
  const byKey = new Map();
  const out = [];
  for (const e of entries) {
    const entry = { smiles: e[0], names: [...e[1]], key: e[2], tier: e[3], cid: e[4] };
    out.push(entry);
    if (entry.smiles && !byKey.has(entry.key)) byKey.set(entry.key, entry);
  }
  const tally = {};
  // within a list, a compound written with its counter-ions comes before the bare ion that the page
  // would give the same counter-ions (sodium hydroxide before hydroxide), so the salt's name is shown
  const bare = q => (results[q] && results[q].bare ? 1 : 0);
  const order = rows.map((row, q) => q).sort((a, b) => rows[a].tier - rows[b].tier || bare(a) - bare(b) || rows[a].cid - rows[b].cid);
  for (const q of order) {
    const row = rows[q], res = results[q] || { fail: 'lost' };
    const what = res.ok ? 'ok' : res.fail.replace(/:.*/, '');
    tally[what] = (tally[what] || 0) + 1;
    if (!res.ok) {
      if (res.fail !== 'big' && res.fail !== 'read') out.push({ smiles: null, names: row.names.slice(0, 2), key: res.fail, tier: row.tier, cid: row.cid });
      if (process.env.SHOW && /slow|engine|ring/.test(res.fail)) console.log(`  ${res.fail}: ${row.names[0]} ${row.smiles}`);
      continue;
    }
    const cap = row.tier === 0 ? 16 : 5;
    const entry = byKey.get(res.key);
    if (!entry) {
      const made = { smiles: row.smiles, names: row.names.slice(0, cap), key: res.key, tier: row.tier, cid: row.cid };
      out.push(made); byKey.set(res.key, made);
      continue;
    }
    const seen = new Set(entry.names.map(norm));
    const fresh = row.names.filter(nm => !seen.has(norm(nm))).slice(0, cap);
    if (row.tier < entry.tier) {
      entry.smiles = row.smiles; entry.tier = row.tier; entry.cid = row.cid;
      entry.names = [...fresh, ...entry.names];
    } else entry.names.push(...fresh);
  }
  console.log('outcomes', tally);
  // the abbreviations a chemist types, added to the compound they stand for
  for (const [abbr, ...targets] of ALIASES) {
    const want = targets.map(norm);
    const hit = out.find(e => e.smiles && want.includes(norm(e.names[0]))) || out.find(e => e.smiles && e.names.some(nm => want.includes(norm(nm))));
    if (!hit) { console.log(`  no compound for ${abbr} (${targets.join(', ')})`); continue; }
    if (!hit.names.some(nm => norm(nm) === norm(abbr))) hit.names.splice(1, 0, abbr);
  }
  return out.map(e => [e.smiles, e.names, e.key, e.tier, e.cid]);
}


function write(file, entries) {
  entries.sort((a, b) => a[3] - b[3] || a[4] - b[4]);
  const data = {
    about: 'Names, synonyms and SMILES from PubChem, as packaged by the chemicals Python library (MIT), and a few of our own ' +
      '(PubChem number 0). Built by tools/build-names.mjs. Entry: [SMILES or null, names (first is shown), skeleton hash or ' +
      'the reason it is not shown, -1 ours, 0 a common chemical, 1 the rest, PubChem CID].',
    entries,
  };
  const json = JSON.stringify(data);
  const gz = gzipSync(json, { level: 9 });
  writeFileSync(file, gz);
  const shown = entries.filter(e => e[0]).length;
  console.log(`${entries.length} entries (${shown} shown, ${entries.length - shown} with a reason), ${(json.length / 1e6).toFixed(1)} MB, gzipped ${(gz.length / 1e6).toFixed(1)} MB`);
}

// Our own entries (tier -1, PubChem number 0): the ions, cations and radicals of the test list, named
// as a student says them, and the acylium and iminium ions of the ledger (item 41), each with its
// counter-ion, each run through the engine here.
async function curated() {
  const { COMPOUNDS } = await import('../core/compounds.js');
  const { nodeFromSmiles } = await import('../core/smiles.js');
  const { buildNode, enumerate } = await import('../core/lewis.js');
  const { skeletonKey, keyHash } = await import('../core/names.js');
  const extra = [
    { name: 'acylium (Cl-)', smiles: 'CC#[O+].[Cl-]', also: ['acetyl cation', 'acetylium'] },
    { name: 'iminium (Cl-)', smiles: 'C=[NH2+].[Cl-]', also: ['methaniminium', 'methyleneiminium'] },
  ];
  const ions = COMPOUNDS.filter(c => /\(|cation|anion|radical|tropylium|phenoxide|enolate/.test(c.name));
  const out = [];
  for (const c of [...ions, ...extra]) {
    const base = c.name.replace(/\s*\(.*\)\s*$/, '');
    const names = [base, ...(base !== c.name ? [c.name] : []), ...(c.also || [])];
    const node = buildNode(nodeFromSmiles(c.smiles));
    enumerate(node);
    out.push([c.smiles, names, keyHash(skeletonKey(node)), -1, 0]);
  }
  return out;
}

// lowercase, letters and digits only: how names are compared
function norm(s) { return s.toLowerCase().replace(/[^a-z0-9]/g, ''); }

// names worth searching: not a registry number, not PubChem's lambda notation or index-name tags,
// not too long, with some letters in them
function cleanNames(raw) {
  const out = [];
  const seen = new Set();
  for (let s of raw) {
    if (!s) continue;
    s = s.trim();
    if (!s || s.length > 60) continue;
    if (/^\d{2,7}-\d\d-\d$/.test(s)) continue;
    if (/[$;{}^]|\(\d+ci|ci\)|inverted|\bunii\b|\bnsc\b|\bbrn\b|\bec \d/.test(s)) continue;
    // a language or registry tag ("[polish]", "[inn]"), or words that name a mixture, not a compound
    if (/\[[a-z .,'-]+\]$/.test(s)) continue;
    if (/polymer|solution|mixture|technical|\bgrade\b|reagent|\(natural\)|monomer/.test(s)) continue;
    if (!/[a-z]{3}/.test(s)) continue;
    const k = norm(s);
    if (seen.has(k)) continue;
    seen.add(k);
    out.push(s);
  }
  return out;
}
