const vm=require('vm'),fs=require('fs'),assert=require('assert');
const core=fs.readFileSync(__dirname+'/../em-physics.js','utf8'),app=fs.readFileSync(__dirname+'/../light-polarization.js','utf8');
const ctx={console,document:{getElementById(){return null}}};vm.createContext(ctx);vm.runInContext(core+'\n'+app.slice(0,app.indexOf('function lpModel()')),ctx);const run=s=>vm.runInContext(s,ctx);
run('var cfg={lightTheta:153,lightPhi:22,source:"linear",optic:"hwp",sourceAngle:30,excAngle:20,lightPower:100,detectTheta:47,detectPhi:18,analyzerAngle:35,analyzer:true,filterTransmission:100}; var ns=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(unit)');
assert(Math.abs(run('lpExcitation(cfg).angle')-10)<1e-10);assert.equal(run('lpExcitation(cfg).transmission'),1);
for(let t=0;t<=180;t+=15)for(let p=0;p<360;p+=45){ctx.t=t;ctx.p=p;assert(Math.abs(run('dot(lpAxis(lpFrame(t,p),37),lpFrame(t,p).k)'))<1e-12);assert(Math.abs(run('Math.hypot(...lpAxis(lpFrame(t,p),37))')-1)<1e-12)}
const a=run('lpOptics(cfg,ns)'),b=run('lpOptics({...cfg,excAngle:cfg.excAngle+90},ns)');assert(Math.abs(a.total-b.total)<1e-12);assert(a.total<=a.raw+1e-12);
const off=run('lpOptics({...cfg,analyzer:false},ns)');assert(Math.abs(off.total-off.raw)<1e-12);assert(a.total<off.total);
assert(Math.abs(run('lpOptics({...cfg,lightPower:50},ns).total')-a.total/2)<1e-12);assert(Math.abs(run('lpOptics({...cfg,filterTransmission:20},ns).total')-a.total*.2)<1e-12);
assert(run('lpOptics({...cfg,optic:"polarizer",sourceAngle:0,excAngle:90},ns).total')<1e-25);
assert(Math.abs(run('lpOptics({...cfg,optic:"polarizer",source:"unpolarized"},ns).exc.transmission')-.5)<1e-12);
const circ=run('lpOptics({...cfg,source:"circular",optic:"none"},ns)'),unpol=run('lpOptics({...cfg,source:"unpolarized",optic:"hwp"},ns)');assert(Math.abs(circ.total-unpol.total)<1e-12);
assert(Math.abs(run('lpOptics({...cfg,optic:"hwp"},ns).contributions.reduce((s,x)=>s+x.excitation,0)')-8/3)<1e-12);
run('var field=[3.54,1.73,6.95]; var before=ns.map(n=>lineAtField(field,n,[1,0,0])).map(l=>l.energies);lpOptics({...cfg,excAngle:125,lightPhi:62},ns);var after=ns.map(n=>lineAtField(field,n,[1,0,0])).map(l=>l.energies)');assert.equal(JSON.stringify(run('before')),JSON.stringify(run('after')));
assert(run('lpOptics({...cfg,source:"linear",optic:"none",sourceAngle:0,lightTheta:0,lightPhi:0},[[1,0,0],[0,1,0],[0,0,1],[0,1,0]]).contributions[0].excitation')<1e-12);
console.log('PASS HWP 2α−β rotation and power conservation; transverse fields at poles and arbitrary directions; 90° HWP periodicity; analyzer and filter bounds; Malus extinction; circular/unpolarized weights; tetrahedral absorption sum; unchanged spin energies.');
