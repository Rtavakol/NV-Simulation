const fs=require('fs'),vm=require('vm'),assert=require('assert');
const html=fs.readFileSync(__dirname+'/../electromagnet-studio.html','utf8'),els={};
for(const m of html.matchAll(/<(input|select)\b([^>]*)>([\s\S]*?)(?:<\/select>|(?=<))/g)){}
for(const m of html.matchAll(/<(input|select)\b([^>]*)>/g)){const id=m[2].match(/id="([^"]+)/)?.[1];if(!id)continue;let value=m[2].match(/value="([^"]*)/)?.[1]||'';if(m[1]==='select'){const rest=html.slice(m.index+m[0].length).split('</select>')[0];value=(rest.match(/<option value="([^"]*)" selected/)||rest.match(/<option value="([^"]*)"/))?.[1]||''}els[id]={value,disabled:false,hidden:false,checkValidity(){return true},parentElement:{textContent:id},textContent:'',innerHTML:'',clientWidth:1000};}
const element=id=>els[id]||(els[id]={value:'',textContent:'',innerHTML:'',hidden:false,disabled:false,checkValidity(){return true}});
const ctx={console,Math,JSON,Number,Array,Set,Map,Date,Error,Infinity,sessionStorage:{getItem(){return null}},document:{getElementById:element,querySelectorAll(){return[]}},window:{addEventListener(){},devicePixelRatio:1},setTimeout,clearTimeout,location:{},performance};vm.createContext(ctx);vm.runInContext(fs.readFileSync(__dirname+'/../em-physics.js','utf8')+'\n'+fs.readFileSync(__dirname+'/../em-studio.js','utf8'),ctx);
const run=s=>vm.runInContext(s,ctx);
run('emShow=()=>{};emCoil=()=>{}');
(async()=>{
console.time('tests');
run('var c=emConfig(); emModel(c); var ns=emAxes();');
const awg=run('emWinding(emCoilConfig(),8)');assert(Math.abs(awg.d-.25464)<.001);assert(awg.N>0);assert(awg.R>awg.R20);assert(Math.abs(awg.K_mT_A*awg.I-8)<1e-9);assert(Math.abs(awg.V-awg.I*awg.R)<1e-9);
const ideal=run('emWinding({...emCoilConfig(),wireOD:.285,packing:1,coilLength:.285,coilThickness:.285,leads:0},1)');assert.equal(ideal.N,1);const radius=(70/2+.285/2)/1000;assert(Math.abs(ideal.K_mT_A-4*Math.PI*1e-7/(2*radius)*1000)<1e-10);
const pp=run('emPseudo([[1,0,0],[0,2,0],[0,0,3],[1,1,1]])');assert.equal(pp.rank,3);const inp=[[1,0,0],[0,2,0],[0,0,3],[1,1,1]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)assert(Math.abs(pp.P[i].reduce((s,x,k)=>s+x*inp[k][j],0)-(i===j?1:0))<1e-9);
const pat=run('emPattern([{f:0,w:1}],{...c,fmdev:.01,linewidth:1.3,lintol:10})');assert(pat.half>.01&&pat.half<.3);assert(Math.abs(pat.left-pat.right)<1e-8);console.log('PASS AWG, discrete-loop field, Ohm law, temperature resistance, pseudoinverse identity, slope-defined linear interval',pat.half);
const wc=run('({...emCoilConfig(),vmax:100,turnMode:"minimum"})');ctx.wc=wc;
const minCoil=run('emWinding(wc,3)');const filled=run('emWinding({...wc,turnMode:"fill"},3)');assert(minCoil.N<filled.N);assert(minCoil.I<=minCoil.wireCurrentCap+1e-12);assert(Math.abs(minCoil.R-minCoil.wire*minCoil.ohmsPerM)<1e-10);ctx.nBefore=minCoil.N-1;
const iBefore=run('3/1000/EM.windingCache.K[nBefore]');assert(iBefore>minCoil.wireCurrentCap);
const cold=run('emWinding({...wc,copperTemp:20},3)');assert(Math.abs(minCoil.R/cold.R-(1+.00393*40))<1e-10);
const fine=run('emWire(wc)'),thick=run('emWire({...wc,awg:28})');assert(thick.reference>fine.reference);assert(thick.ohmsPerM20<fine.ohmsPerM20);assert(Math.abs(fine.reference-Math.pow(fine.d/.0254,2)/700)<1e-12);
assert.throws(()=>run('emWinding({...wc,wireOD:.01},3)'),/smaller/);
console.log('PASS minimum turns, AWG-based current cap, previous-turn rejection, gauge and temperature resistance updates');

console.time('search');await run('emSearch()');console.timeEnd('search');console.log(element('searchStatus').textContent);const r=run('EM.result');assert(r&&r.separation);console.log({mag:r.mag,capture:r.capture,linearIsolation:r.linearIsolation,drift:Math.max(...r.reserved),linear:Math.min(...r.patterns.map(p=>p.used)),condition:r.pseudo.condition,samples:r.samples,gap:Math.min(...r.gaps.map(g=>g.gap)),required:Math.max(...r.gaps.map(g=>g.required))});
assert(r.gaps.every(g=>g.gap>=g.required));assert(r.pseudo.rank===3);assert(r.A.length===8);assert(r.samples>=81);
// Zero dynamic/thermal range must pass linear capture.
const A=r.selected.map(i=>r.A[i]);for(let i=0;i<3;i++)for(let j=0;j<3;j++)assert(Math.abs(r.pseudo.P[i].reduce((sum,x,k)=>sum+x*A[k][j],0)-(i===j?1:0))<1e-8);
const overlap=await run('emEvaluate({mag:EM.result.mag,dir:EM.result.dir},{...c,manualLinear:30},ns)');assert(!overlap.linearIsolation);assert(!overlap.feasible);
const rz=await run('emEvaluate({mag:EM.result.mag,dir:EM.result.dir},{...c,ex:0,ey:0,ez:0,tmin:20,tmax:20,tref:20},ns)');assert(rz.capture);assert(rz.reserved.every(x=>x<1e-7));assert(rz.separation);assert(rz.feasible);assert.throws(()=>run('emWinding({...emCoilConfig(),wireOD:.01},8)'),/smaller/);console.log('PASS default full-model search and zero-range capture');console.timeEnd('tests');
})().catch(e=>{console.error(e);process.exitCode=1});
