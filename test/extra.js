import { buildNode, enumerate } from '../core/lewis.js';
import { relax, geometryReport } from '../core/geometry.js';
const cases = {
  ethane: { atoms: 'C C H H H H H H'.split(' ').map(el=>({el})), adj: [[0,1],[0,2],[0,3],[0,4],[1,5],[1,6],[1,7]], charge: 0 },
  butadiene: { atoms: 'C C C C H H H H H H'.split(' ').map(el=>({el})), adj: [[0,1],[1,2],[2,3],[0,4],[0,5],[1,6],[2,7],[3,8],[3,9]], charge: 0 },
  cyclohexane: { atoms: 'C C C C C C H H H H H H H H H H H H'.split(' ').map(el=>({el})), adj: [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0],[0,6],[0,7],[1,8],[1,9],[2,10],[2,11],[3,12],[3,13],[4,14],[4,15],[5,16],[5,17]], charge: 0 },
  ammonia: { atoms: 'N H H H'.split(' ').map(el=>({el})), adj: [[0,1],[0,2],[0,3]], charge: 0 },
  PCl5: { atoms: 'P Cl Cl Cl Cl Cl'.split(' ').map(el=>({el})), adj: [[0,1],[0,2],[0,3],[0,4],[0,5]], charge: 0 },
  XeF4: { atoms: 'Xe F F F F'.split(' ').map(el=>({el})), adj: [[0,1],[0,2],[0,3],[0,4]], charge: 0 },
};
for (const [name, c] of Object.entries(cases)) {
  try {
    const node = buildNode(c); const r = enumerate(node);
    const t0 = performance.now(); const geo = relax(node, r);
    console.log(`\n== ${name} (${(performance.now()-t0).toFixed(0)} ms) structures ${r.structures.length} holes ${r.holes}`);
    console.log(geometryReport(node, r, geo));
    // dihedral check for butadiene C0-C1-C2-C3
    if (name === 'butadiene' || name === 'ethane') {
      const p = geo.pos; const sub=(a,b)=>a.map((x,i)=>x-b[i]); const cr=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]; const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0); const ln=a=>Math.hypot(...a);
      const dih=(a,b,c,d)=>{const b1=sub(p[b],p[a]),b2=sub(p[c],p[b]),b3=sub(p[d],p[c]);const n1=cr(b1,b2),n2=cr(b2,b3);return Math.acos(dot(n1,n2)/(ln(n1)*ln(n2)))*180/Math.PI;};
      if (name==='butadiene') console.log('dihedral C0-C1-C2-C3', dih(0,1,2,3).toFixed(1));
      if (name==='ethane') console.log('dihedral H2-C0-C1-H5', dih(2,0,1,5).toFixed(1), 'H3-C0-C1-H5', dih(3,0,1,5).toFixed(1));
    }
  } catch (e) { console.log(name, 'ERROR', e.message); }
}
