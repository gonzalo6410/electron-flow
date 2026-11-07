// Build the name search index (web/names.json.gz) from the compound lists of the `chemicals` Python
// package (MIT; Caleb Bell), whose names, synonyms and SMILES come from PubChem (public domain).
//
//   pip download chemicals --no-deps -d /tmp/chem && (cd /tmp/chem && python3 -m zipfile -e chemicals-*.whl x)
//   node tools/build-names.mjs /tmp/chem/x/chemicals/Identifiers
//
// Every compound is run through the engine as the page would run it (read, counter-ions added to a
// bare ion, structures found, families, shape). Compounds the engine finishes within the time limit go
// in with their SMILES; the rest go in with their main names and the reason, so a search for them says
// why instead of saying nothing (an element not in the table yet, too slow for the engine yet).
// Compounds that are the same molecule to the engine (stereoisomers, which it ignores) are merged into
// one entry, and each entry carries a hash of the molecule's skeleton so the page can name what is on
// screen after an edit or a move (core/names.js).

import { isMainThread, Worker, parentPort } from 'node:worker_threads';
import { readFileSync, writeFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const LIMIT_MS = 1500;      // engine time allowed per compound (the page is slower than node)
const MAX_HEAVY = 70;       // heavier molecules are left out of the search altogether

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
      try { spec = nodeFromSmiles(balanceSmiles(smiles).smiles); } catch (e) {
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
      parentPort.postMessage({ id, ok: true, key: keyHash(skeletonKey(node)), n: r.structures.length, ms: performance.now() - t0 });
    } catch (e) {
      parentPort.postMessage({ id, fail: /ring/.test(e.message) ? 'ring' : 'engine', msg: e.message });
    }
  });
} else {
  const dir = process.argv[2];
  if (!dir) { console.error('usage: node tools/build-names.mjs <chemicals/Identifiers folder>'); process.exit(1); }
  const files = [['chemical identifiers pubchem small.tsv', 0], ['chemical identifiers pubchem large.tsv', 1]];
  const rows = new Map(); // cid -> row
  for (const [file, tier] of files) {
    for (const line of readFileSync(path.join(dir, file), 'utf8').split('\n')) {
      const f = line.split('\t');
      if (f.length < 9 || !f[4]) continue;
      const cid = Number(f[0]);
      if (rows.has(cid)) continue;
      rows.set(cid, { cid, tier, smiles: f[4], names: cleanNames([f[8], f[7], ...f.slice(9)]) });
    }
  }
  const list = [...rows.values()].filter(r => r.names.length);
  console.log(`${list.length} compounds read`);

  const results = new Array(list.length);
  const here = fileURLToPath(import.meta.url);
  let next = 0, done = 0;
  const t0 = Date.now();
  await new Promise(resolve => {
    const lanes = 2;
    let live = 0;
    const startLane = () => {
      live++;
      const w = new Worker(here);
      let current = -1, timer = null;
      const feed = () => {
        if (next >= list.length) { w.terminate(); live--; if (!live) resolve(); return; }
        current = next++;
        timer = setTimeout(() => {
          results[current] = { fail: 'slow' };
          done++;
          w.terminate();
          live--;
          startLane();
        }, LIMIT_MS + 2000); // reading and loading the worker take a moment of their own
        w.postMessage({ id: current, smiles: list[current].smiles });
      };
      w.on('message', msg => {
        clearTimeout(timer);
        results[msg.id] = msg.ok && msg.ms > LIMIT_MS ? { fail: 'slow' } : msg;
        if (++done % 2000 === 0) console.log(`${done} done, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
        feed();
      });
      w.on('error', e => { console.error('worker error', e.message); });
      feed();
    };
    for (let q = 0; q < lanes; q++) startLane();
  });

  // merge compounds that are one molecule to the engine; the common-chemicals list first, then the
  // lowest PubChem number, gives the entry its SMILES and its display name
  const byKey = new Map();
  const failed = [];
  const tally = {};
  list.forEach((row, q) => {
    const res = results[q] || { fail: 'lost' };
    tally[res.ok ? 'ok' : res.fail.replace(/:.*/, '')] = (tally[res.ok ? 'ok' : res.fail.replace(/:.*/, '')] || 0) + 1;
    if (!res.ok) { if (res.fail !== 'big' && res.fail !== 'read') failed.push({ row, why: res.fail }); return; }
    const k = res.key;
    if (!byKey.has(k)) byKey.set(k, []);
    byKey.get(k).push(row);
  });
  console.log('outcomes', tally);
  if (process.env.SHOW) for (const { row, why } of failed) if (why === 'slow' || why === 'engine' || why === 'ring') console.log(`  ${why}: ${row.names[0]} ${row.smiles}`);
  const entries = [];
  for (const [k, group] of byKey) {
    group.sort((a, b) => a.tier - b.tier || a.cid - b.cid);
    const names = [];
    const seen = new Set();
    for (const g of group) {
      const cap = g.tier === 0 ? 16 : 5;
      let taken = 0;
      for (const nm of g.names) {
        const flat = norm(nm);
        if (!flat || seen.has(flat)) continue;
        if (taken >= cap) break;
        seen.add(flat); names.push(nm); taken++;
      }
    }
    entries.push([group[0].smiles, names, k, group[0].tier, group[0].cid]);
  }
  // failures keep two names each and the reason
  for (const { row, why } of failed) entries.push([null, row.names.slice(0, 2), why, row.tier, row.cid]);
  entries.sort((a, b) => a[3] - b[3] || a[4] - b[4]);
  const out = {
    about: 'Names, synonyms and SMILES from PubChem, as packaged by the chemicals Python library (MIT). ' +
      'Built by tools/build-names.mjs. Entry: [SMILES or null, names (first is shown), skeleton hash or the reason ' +
      'it is not shown, 0 if a common chemical, PubChem CID].',
    entries: entries.map(e => [e[0], e[1], e[2], e[3], e[4]]),
  };
  const json = JSON.stringify(out);
  const gz = gzipSync(json, { level: 9 });
  writeFileSync(path.join(path.dirname(here), '..', 'web', 'names.json.gz'), gz);
  console.log(`${entries.length} entries (${byKey.size} shown, ${failed.length} with a reason), ${(json.length / 1e6).toFixed(1)} MB, gzipped ${(gz.length / 1e6).toFixed(1)} MB`);
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
