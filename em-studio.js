'use strict';
const EM={running:false,cancel:false,result:null,context:null,coils:null};
const emNum=id=>{const el=$(id),v=+el.value;if(!el.checkValidity()||el.value===''||!Number.isFinite(v))throw Error('Check '+el.parentElement.textContent.trim());return v};
const emYield=()=>new Promise(r=>setTimeout(r,0));
function emConfig(){
 const c={};for(const id of ['ex','ey','ez','tmin','tmax','tref','linewidth','fmdev','lintol','manualLinear','guard','reserve','yaw','pitch','roll','ed','eapar','eaperp','eq','eth','eph','directions','bmin','bmax','bstep','egn'])c[id]=emNum(id);
 for(const id of ['domain','tones','strategy','subset','ecut','direction'])c[id]=$(id).value;
 if(c.tmin>c.tmax||c.tref<c.tmin||c.tref>c.tmax)throw Error('Reference temperature must be inside the ordered temperature range.');
 if(c.bmin>=c.bmax)throw Error('Search maximum must exceed minimum.');return c;
}
function emModel(c){cut=c.ecut;Object.assign(state,{rotation:c.yaw,tiltY:c.pitch,tiltX:c.roll,D:c.ed,apar:c.eapar,aperp:c.eaperp,quad:c.eq,gammaN:c.egn,temp:c.tref});mwMode='isotropic'}
function emAxes(){return [[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(v=>crystalToLab(unit(v)))}
function emDirection(theta,phi){return [Math.sin(theta*RAD)*Math.cos(phi*RAD),Math.sin(theta*RAD)*Math.sin(phi*RAD),Math.cos(theta*RAD)]}
function emElectronic(field,ns,temp,c){const D=c.ed+thermalPolynomial(temp)-thermalPolynomial(25);return ns.flatMap((n,axis)=>{const bp=dot(n,field),bt=Math.sqrt(Math.max(0,dot(field,field)-bp*bp)),z=28.025*bp,x=28.025*bt*Math.SQRT1_2,e=eigensystem([[D+z,x,0],[x,0,x],[0,x,D-z]]);let b=0;for(let i=1;i<3;i++)if(e.vectors[i][1]**2>e.vectors[b][1]**2)b=i;return e.values.map((v,i)=>({f:Math.abs(v-e.values[b]),axis,branch:i})).filter(v=>v.branch!==b)}).sort((a,b)=>a.f-b.f)}
function emDomainBound(row,c){const v=row.map((x,i)=>x*[c.ex,c.ey,c.ez][i]);return c.domain==='box'?v.reduce((s,x)=>s+Math.abs(x),0):Math.hypot(...v)}
function emScreen(mag,dir,c,ns){
 const field=dir.map(x=>mag*x),r=emElectronic(field,ns,c.tref,c),rows=r.map(()=>[]),h=.0001;
 for(let k=0;k<3;k++){const p=field.slice(),m=field.slice();p[k]+=h;m[k]-=h;const a=emElectronic(p,ns,c.tref,c),b=emElectronic(m,ns,c.tref,c);for(let i=0;i<8;i++){if(a[i].axis!==r[i].axis||b[i].axis!==r[i].axis)return null;rows[i][k]=(a[i].f-b[i].f)/(2*h*1000)}}
 const dt=Math.max(Math.abs(thermalPolynomial(c.tmin)-thermalPolynomial(c.tref)),Math.abs(thermalPolynomial(c.tmax)-thermalPolynomial(c.tref))),extra=Math.abs(c.eapar)*(c.tones==='triple'?2:1)+c.fmdev;
 const widths=rows.map(row=>c.reserve*(emDomainBound(row,c)+(c.strategy==='fixed'?dt:0))+extra);
 const surplus=Math.min(...r.slice(1).map((x,i)=>x.f-r[i].f-widths[i]-widths[i+1]-c.guard*c.linewidth));
 return {mag,dir,surplus};
}
function emNominal(field,temp,c,ns){state.temp=temp;const ls=ns.map(n=>lineAtField(field,n,[1,0,0]));const tracks=ls.flatMap((l,axis)=>{const refs=resonanceCentres(l);return refs.map((tr,b)=>{const other=refs[1-b];const group=l.transitions.filter(t=>Math.abs(t.f-tr.f)<8&&Math.abs(t.f-tr.f)<=Math.abs(t.f-other.f)+1e-7);const peak=Math.max(...group.map(t=>t.weight),1e-12);return {axis,i:tr.i,j:tr.j,f:tr.f,weight:tr.weight,group:group.filter(t=>t.weight>=peak*.05).map(t=>({f:t.f-tr.f,w:t.weight}))}})}).sort((a,b)=>a.f-b.f);return {ls,tracks,field,temp};}
function emTracked(field,temp,nom,c,ns){state.temp=temp;const ls=ns.map(n=>lineAtField(field,n,[1,0,0])),maps=ls.map((l,i)=>overlapAssignment(nom.ls[i],l));return nom.tracks.map(tr=>{const i=maps[tr.axis][tr.i],j=maps[tr.axis][tr.j],l=ls[tr.axis];return {f:Math.abs(l.energies[i]-l.energies[j]),weight:l.transitions.find(t=>t.i===Math.min(i,j)&&t.j===Math.max(i,j))?.weight||0}})}
function emJacobian(nom,c,ns){const A=Array.from({length:8},()=>[]),h=.0001;for(let k=0;k<3;k++){const p=nom.field.slice(),m=p.slice();p[k]+=h;m[k]-=h;const a=emTracked(p,c.tref,nom,c,ns),b=emTracked(m,c.tref,nom,c,ns);for(let i=0;i<8;i++)A[i][k]=(a[i].f-b[i].f)/(2*h*1000)}const p=emTracked(nom.field,c.tref+.01,nom,c,ns),m=emTracked(nom.field,c.tref-.01,nom,c,ns);return {A,kT:p.map((x,i)=>(x.f-m[i].f)/.02)}}
function emPseudo(A){const ata=Array.from({length:3},(_,i)=>Array.from({length:3},(_,j)=>A.reduce((s,row)=>s+row[i]*row[j],0))),e=eigensystem(ata),max=Math.max(...e.values),rank=e.values.filter(x=>x>max*1e-10).length;const inv=Array.from({length:3},(_,i)=>Array.from({length:3},(_,j)=>e.values.reduce((s,val,k)=>s+(val>max*1e-10?e.vectors[k][i]*e.vectors[k][j]/val:0),0)));return {rank,condition:rank===3?Math.sqrt(max/Math.min(...e.values)):Infinity,P:inv.map(row=>A.map(a=>dot(row,a)))}}
function emResponse(x,group,c){const N=Math.min(8192,Math.max(128,Math.ceil(24*c.fmdev/c.linewidth))),w=c.linewidth/2,offsets=c.tones==='triple'?[-Math.abs(c.eapar),0,Math.abs(c.eapar)]:[0];let sum=0;for(let j=0;j<N;j++){const u=Math.cos(2*Math.PI*j/N);for(const tone of offsets)for(const tr of group){const d=x+tone+c.fmdev*u-tr.f;sum-=tr.w*w*w/(w*w+d*d)*u*2/N}}return sum}
function emPattern(group,c){
 const dh=.002,slope=x=>(emResponse(x+dh,group,c)-emResponse(x-dh,group,c))/(2*dh),s0=slope(0),tol=c.lintol/100,step=.02;
 function edge(sign){let prev=0;if(Math.abs(s0)<1e-12)return 0;for(let x=step;x<=30;x+=step){if(Math.abs(slope(sign*x)/s0-1)>tol){let lo=prev,hi=x;for(let k=0;k<12;k++){const mid=(lo+hi)/2;if(Math.abs(slope(sign*mid)/s0-1)<=tol)lo=mid;else hi=mid}return lo}prev=x}return 30}
 const left=edge(-1),right=edge(1),half=Math.min(left,right);return {left,right,half,s0,group,measured:c.manualLinear>0,used:c.manualLinear||half};
}
function emSamples(c){const pts=[];if(c.domain==='box'){for(const a of [-1,0,1])for(const b of [-1,0,1])for(const d of [-1,0,1])pts.push([a*c.ex,b*c.ey,d*c.ez])}else{pts.push([0,0,0]);for(let k=0;k<3;k++)for(const sign of [-1,1]){const v=[0,0,0];v[k]=[c.ex,c.ey,c.ez][k]*sign;pts.push(v)}for(let k=0;k<42;k++){const z=1-2*(k+.5)/42,a=k*2.3999632297,r=Math.sqrt(1-z*z);pts.push([c.ex*r*Math.cos(a),c.ey*r*Math.sin(a),c.ez*z])}}return pts}
async function emEvaluate(candidate,c,ns){
 const field=candidate.dir.map(x=>candidate.mag*x),nom=emNominal(field,c.tref,c,ns),{A,kT}=emJacobian(nom,c,ns),selected=[...c.subset].map(x=>+x-1),pseudo=emPseudo(selected.map(i=>A[i]));
 const patterns=nom.tracks.map(tr=>emPattern(tr.group,c));const fieldBound=A.map(row=>emDomainBound(row,c)),temps=[...new Set([c.tmin,(c.tmin+c.tmax)/2,c.tmax,c.tref])],thermal=temps.map(t=>emTracked(field,t,nom,c,ns).map((x,i)=>x.f-nom.tracks[i].f)),thermalBound=A.map((_,i)=>Math.max(...thermal.map(r=>Math.abs(r[i]))));
 let drift=A.map((_,i)=>fieldBound[i]+(c.strategy==='fixed'?thermalBound[i]:0)),minimumGap=Infinity,observable=true,maxReconstruction=0,samples=0;const baseWeight=Math.max(...nom.tracks.map(t=>t.weight)),pts=emSamples(c);
 for(let ti=0;ti<temps.length;ti++)for(let pi=0;pi<pts.length;pi++){
  if(EM.cancel)throw Error('Search cancelled.');const p=pts[pi],values=emTracked(field.map((x,i)=>x+p[i]/1000),temps[ti],nom,c,ns),shifts=values.map((v,i)=>v.f-nom.tracks[i].f-(c.strategy==='temperature'?thermal[ti][i]:0));
  drift=drift.map((v,i)=>Math.max(v,Math.abs(shifts[i])));for(let i=0;i<7;i++)minimumGap=Math.min(minimumGap,values[i+1].f-values[i].f);if(values.some(v=>v.weight<baseWeight*.01))observable=false;
  const recon=pseudo.P.map(row=>dot(row,selected.map(i=>shifts[i])));maxReconstruction=Math.max(maxReconstruction,Math.hypot(...recon.map((x,i)=>x-p[i])));samples++;
  if(samples%20===0)await emYield();
 }
 const reserved=drift.map(x=>x*c.reserve),toneExtra=c.tones==='triple'?Math.abs(c.eapar):0;
 const hf=nom.tracks.map(tr=>Math.max(Math.abs(c.eapar),...tr.group.map(t=>Math.abs(t.f)))),alloc=reserved.map((x,i)=>x+hf[i]+c.fmdev+toneExtra),gaps=nom.tracks.slice(1).map((t,i)=>({gap:t.f-nom.tracks[i].f,required:alloc[i]+alloc[i+1]+c.guard*c.linewidth,linear:patterns[i].used+patterns[i+1].used+reserved[i]+reserved[i+1]}));
 const separation=gaps.every(g=>g.gap>=g.required)&&observable&&pseudo.rank===3;
 const capture=reserved.every((d,i)=>d<=patterns[i].used),linearIsolation=gaps.every(g=>g.gap>=g.linear);
 state.temp=c.tref;return {...candidate,nom,A,kT,selected,pseudo,patterns,fieldBound,thermalBound,drift,reserved,hf,alloc,gaps,separation,capture,linearIsolation,feasible:separation&&capture&&linearIsolation,minimumGap,observable,samples,maxReconstruction,c};
}
async function emSearch(){
 if(EM.running)return;let c;try{c=emConfig()}catch(e){$('searchStatus').textContent=e.message;return}EM.running=true;EM.cancel=false;$('findField').disabled=true;$('cancelSearch').disabled=false;EM.result=null;EM.context=null;$('designResults').hidden=$('matrixSection').hidden=true;
 document.querySelectorAll('#requirements input,#requirements select').forEach(el=>el.disabled=true);
 try{
 emModel(c);const ns=emAxes(),dirs=[emDirection(c.eth,c.eph)];if(c.direction==='search')for(let k=0;k<c.directions;k++){const z=(k+.5)/c.directions,a=k*2.3999632297,r=Math.sqrt(1-z*z);dirs.push([r*Math.cos(a),r*Math.sin(a),z])}
 const candidates=[],coarse=Math.max(.2,c.bstep);for(let d=0;d<dirs.length;d++){
  if(EM.cancel)throw Error('Search cancelled.');for(let mag=c.bmin;mag<=c.bmax+1e-9;mag+=coarse){const a=emScreen(mag,dirs[d],c,ns);if(a&&a.surplus>=0){candidates.push(a);break}}
  if(d%4===0){$('searchStatus').textContent='Screening direction '+(d+1)+' / '+dirs.length+'…';await emYield()}
 }
 candidates.sort((a,b)=>a.mag-b.mag||b.surplus-a.surplus);if(!candidates.length)throw Error('No separated electronic candidate in this search range. Increase the maximum bias or reduce the frequency allocation.');
 let best=null,split=null;const shortlist=candidates.slice(0,c.direction==='fixed'?1:6);
 for(let j=0;j<shortlist.length;j++){
  const seed=shortlist[j];let verified=null;for(let k=0;k<7;k++){const mag=Math.min(c.bmax,seed.mag+k*.25);$('searchStatus').textContent='Full ¹⁴N validation '+(j+1)+' / '+shortlist.length+' at '+mag.toFixed(3)+' mT…';await emYield();const r=await emEvaluate({...seed,mag},c,ns);if(r.separation){verified=r;break}if(mag===c.bmax)break}
  if(!verified)continue;
  // Descend on the chosen magnitude grid until the first failed full-model test.
  let mag=Math.floor((verified.mag-c.bmin)/c.bstep-1e-7)*c.bstep+c.bmin;
  for(let k=0;k<30&&mag>=Math.max(c.bmin,seed.mag-.5)-1e-8;k++,mag-=c.bstep){$('searchStatus').textContent='Refining direction '+(j+1)+' at '+mag.toFixed(3)+' mT…';await emYield();const r=await emEvaluate({...seed,mag},c,ns);if(!r.separation)break;verified=r}
  if(!split||verified.mag<split.mag)split=verified;
  let usable=verified.feasible?verified:null;
  if(!usable){let lower=verified.mag;const probes=[.5,1,2,4,8,c.bmax-verified.mag].map(delta=>Math.min(c.bmax,verified.mag+delta)).filter((m,i,a)=>m>verified.mag&&a.indexOf(m)===i).sort((a,b)=>a-b);
   for(const m of probes){if(best&&m>=best.mag)break;$('searchStatus').textContent='Checking linear feasibility for direction '+(j+1)+' at '+m.toFixed(3)+' mT…';await emYield();const r=await emEvaluate({...seed,mag:m},c,ns);if(r.feasible){usable=r;break}lower=m}
   if(usable){let upper=usable.mag;while(upper-lower>c.bstep){const mid=(lower+upper)/2;const r=await emEvaluate({...seed,mag:mid},c,ns);if(r.feasible){usable=r;upper=mid}else lower=mid}}
  }
  if(usable&&(!best||usable.mag<best.mag))best=usable;
 }
 const result=best||split;if(!result)throw Error('Electronic screen found candidates, but none passed full-model separation/observability. Try a finer direction grid or wider field range.');
 result.searchedDirections=dirs.length;result.validatedDirections=shortlist.length;EM.result=result;EM.context=JSON.stringify(c);emShow(result);$('searchStatus').textContent=result.feasible?'Search complete: sampled separation, linear capture and linear-interval isolation passed.':'Search complete: separation candidate found, but the requested linear operating range is not feasible. See the failure details.';$('coilTargetSource').value='search';emCoil();
 }catch(e){$('searchStatus').textContent=e.message}finally{state.temp=c.tref;EM.running=false;$('findField').disabled=false;$('cancelSearch').disabled=true;document.querySelectorAll('#requirements input,#requirements select').forEach(el=>el.disabled=false)}
}
function emFresh(){if(!EM.result)throw Error('Run the field search first.');if(EM.context!==JSON.stringify(emConfig()))throw Error('Inputs changed. Run the field search again.');return EM.result}
const emFmt=(v,n=3)=>Number.isFinite(v)?v.toFixed(n):'∞';
function emShow(r){
 $('designResults').hidden=$('matrixSection').hidden=false;const bad=[];if(!r.capture)bad.push('Field/temperature drift exceeds at least one linear half-range. Increasing bias does not guarantee a wider linear interval: change FM deviation/linewidth, reduce the required range, or use carrier retuning.');if(!r.linearIsolation)bad.push('Adjacent linear intervals overlap. Linear half-ranges plus reserved drift must fit between adjacent resonance centres.');
 $('fieldSummary').innerHTML='<p class="'+(r.feasible?'ok':'warn')+'"><b>'+(r.feasible?'Candidate meets the modeled checks.':'Separation-only candidate: no fully feasible design found.')+'</b></p>'+bad.map(s=>'<p class="warn">'+s+'</p>').join('');
 $('fieldLabel').textContent=r.feasible?'Smallest verified usable candidate':'Smallest verified separation candidate';$('bestField').textContent=emFmt(r.mag)+' mT';$('driftReadout').textContent=emFmt(Math.max(...r.reserved))+' MHz';$('linearReadout').textContent=emFmt(Math.min(...r.patterns.map(p=>p.used)))+' MHz';
 const theta=Math.acos(r.dir[2])/RAD,phi=(Math.atan2(r.dir[1],r.dir[0])/RAD+360)%360;$('directionResult').textContent='Lab θ = '+emFmt(theta,2)+'°, φ = '+emFmt(phi,2)+'°. Bias [Bx, By, Bz] = ['+r.nom.field.map(x=>emFmt(x)).join(', ')+'] mT. '+r.searchedDirections+' directions screened; '+r.validatedDirections+' shortlisted for full-model checks. This is the smallest verified candidate found, not a certified global optimum.';
 $('budgetRows').innerHTML=r.nom.tracks.map((tr,i)=>'<tr><td>R'+(i+1)+' · NV'+(tr.axis+1)+'</td><td>'+emFmt(tr.f)+'</td><td>'+emFmt(r.fieldBound[i])+'</td><td>'+emFmt(r.thermalBound[i])+'</td><td>'+emFmt(r.reserved[i])+'</td><td>'+emFmt(r.patterns[i].used)+'</td><td class="'+(r.reserved[i]<=r.patterns[i].used?'ok':'warn')+'">'+(r.reserved[i]<=r.patterns[i].used?'PASS':'FAIL')+'</td></tr>').join('');
 $('gapRows').innerHTML=r.gaps.map((g,i)=>'<tr><td>R'+(i+1)+' / R'+(i+2)+'</td><td>'+emFmt(g.gap)+'</td><td>'+emFmt(g.required)+'</td><td>'+emFmt(g.gap-g.required)+'</td><td>'+emFmt(g.linear)+'</td><td class="'+(g.gap>=g.linear?'ok':'warn')+'">'+(g.gap>=g.linear?'PASS':'FAIL')+'</td></tr>').join('');
 $('validationResult').textContent=r.samples+' field/temperature combinations sampled. Minimum instantaneous centre gap = '+emFmt(r.minimumGap)+' MHz. All eight representative transitions observable: '+(r.observable?'yes':'no')+'. The reserve factor covers a design allowance, not a mathematical bound on all unsampled nonlinear cases.';
 $('matrixRows').innerHTML=r.A.map((row,i)=>'<tr><td>R'+(i+1)+'</td>'+row.map(x=>'<td>'+emFmt(x,6)+'</td>').join('')+'<td>'+emFmt(r.kT[i],6)+'</td></tr>').join('');
 $('pinvResult').textContent='Columns: '+r.selected.map(i=>'R'+(i+1)).join('   ')+'\n'+r.pseudo.P.map((row,i)=>['Bx','By','Bz'][i]+'  '+row.map(x=>emFmt(x,5).padStart(12)).join(' ')).join('\n');
 $('conditionResult').textContent='Selected A rank: '+r.pseudo.rank+' / 3. Singular-value condition number: '+emFmt(r.pseudo.condition,2)+'. Units: A in MHz/µT, A⁺ in µT/MHz.';
 const leak=r.pseudo.P.map(row=>dot(row,r.selected.map(i=>r.kT[i])));$('leakageResult').textContent='Uncompensated local temperature leakage A⁺kT = ['+leak.map(x=>emFmt(x,3)).join(', ')+'] µT/°C for [Bx, By, Bz].';
 $('reconstructionResult').textContent='Largest sampled vector reconstruction error using this fixed A⁺ and '+(r.c.strategy==='temperature'?'ideal thermal feed-forward':'uncompensated frequency shifts')+' = '+emFmt(r.maxReconstruction,3)+' µT. This includes nonlinear model error and thermal leakage, but no measurement noise.';
 $('showBranch').innerHTML=r.nom.tracks.map((tr,i)=>'<option value="'+i+'">R'+(i+1)+' · NV'+(tr.axis+1)+'</option>').join('');emDrawPattern();emDrawBudget();
}
function emCanvas(id){const cv=$(id),w=Math.max(600,cv.clientWidth||1000),h=350,d=window.devicePixelRatio||1;cv.width=w*d;cv.height=h*d;const c=cv.getContext('2d');c.scale(d,d);c.fillStyle='#0c1521';c.fillRect(0,0,w,h);c.font='13px system-ui';return {c,w,h}}
function emDrawPattern(){if(!EM.result)return;const r=EM.result,i=+$('showBranch').value||0,p=r.patterns[i],{c,w}=emCanvas('dispersivePlot'),span=Math.max(2,r.reserved[i]*1.3,p.half*2.5),N=260,points=Array.from({length:N+1},(_,j)=>{const x=-span+2*span*j/N,y=emResponse(x,p.group,r.c),d=(emResponse(x+.002,p.group,r.c)-emResponse(x-.002,p.group,r.c))/.004;return {x,y,d:p.s0?d/p.s0:0}}),ym=Math.max(...points.map(v=>Math.abs(v.y)),1e-5),L=65,R=w-25,X=x=>L+(x+span)/(2*span)*(R-L);
 for(const [top,bottom,quantity,scale] of [[30,170,'y',ym*1.15],[210,310,'d',2]]){const Y=y=>bottom-(y+scale)/(2*scale)*(bottom-top);c.fillStyle='#173e36';c.fillRect(X(-p.half),top,X(p.half)-X(-p.half),bottom-top);c.strokeStyle='#50667a';c.beginPath();c.moveTo(L,Y(0));c.lineTo(R,Y(0));c.stroke();if(quantity==='d')for(const val of [1-r.c.lintol/100,1+r.c.lintol/100]){c.strokeStyle='#ffb86b';c.setLineDash([4,4]);c.beginPath();c.moveTo(L,Y(val));c.lineTo(R,Y(val));c.stroke();c.setLineDash([])}c.strokeStyle=quantity==='y'?'#6de4c0':'#70b6ff';c.beginPath();points.forEach((v,j)=>{const y=Math.max(-scale,Math.min(scale,v[quantity]));j?c.lineTo(X(v.x),Y(y)):c.moveTo(X(v.x),Y(y))});c.stroke();c.fillStyle='#adc0d1';c.fillText(quantity==='y'?'FM first harmonic (relative PL)':'Local slope / slope at carrier',L,top-8);c.setLineDash([3,3]);c.strokeStyle='#d49bff';for(const edge of [-r.reserved[i],r.reserved[i]]){c.beginPath();c.moveTo(X(edge),top);c.lineTo(X(edge),bottom);c.stroke()}c.setLineDash([])}c.fillStyle='#adc0d1';c.textAlign='center';for(let j=0;j<=6;j++){const x=-span+2*span*j/6;c.fillText(x.toFixed(2),X(x),331)}c.fillText('Carrier offset / MHz · mint: calculated linear interval · purple: required drift',w/2,348);c.textAlign='left';
 $('patternInfo').textContent='R'+(i+1)+': slope-tolerance interval −'+emFmt(p.left)+' to +'+emFmt(p.right)+' MHz; symmetric half-range '+emFmt(p.half)+' MHz; full symmetric width '+emFmt(2*p.half)+' MHz. '+(p.measured?'Capture checks use your measured half-range of '+emFmt(p.used)+' MHz. ':'')+'FM deviation '+r.c.fmdev+' MHz, FWHM '+r.c.linewidth+' MHz, ±'+r.c.lintol+'% local slope tolerance. Slope units are relative PL/MHz.';
}
function emDrawBudget(){if(!EM.result)return;const r=EM.result,{c,w}=emCanvas('budgetPlot'),lo=Math.min(...r.nom.tracks.map((t,i)=>t.f-r.alloc[i]))-3,hi=Math.max(...r.nom.tracks.map((t,i)=>t.f+r.alloc[i]))+3,L=65,R=w-20,X=f=>L+(f-lo)/(hi-lo)*(R-L);c.fillStyle='#adc0d1';c.fillText('Swept frequency allocations (including reserved drift) · MHz',L,22);r.nom.tracks.forEach((tr,i)=>{const y=48+i*33;c.fillStyle=colors[tr.axis];c.fillRect(X(tr.f-r.alloc[i]),y,Math.max(1,X(tr.f+r.alloc[i])-X(tr.f-r.alloc[i])),18);c.fillStyle='#edf4fa';c.fillText('R'+(i+1),12,y+14);c.fillRect(X(tr.f)-1,y-3,2,24)});for(let j=0;j<=5;j++){const f=lo+(hi-lo)*j/5;c.fillStyle='#adc0d1';c.fillText(f.toFixed(1),X(f)-20,337)}}
function emCoilConfig(){const c={};for(const id of ['awg','coilID','coilLength','coilThickness','wireOD','packing','leads','copperTemp','imax','vmax','headroom','jmax','manualB','sensorHalf','wireDerating','wireManufacturer'])c[id]=emNum(id);if(!Number.isInteger(c.awg))throw Error('AWG must be an integer.');c.turnMode=$('turnMode').value;return c}
function emWire(c){
 const d=.127*Math.pow(92,(36-c.awg)/39),area=Math.PI*d*d/4,od=c.wireOD||d+.03;
 if(od<d)throw Error('Insulated OD cannot be smaller than the bare copper diameter.');
 const reference=Math.pow(d/.0254,2)/700,raw=c.wireManufacturer||reference,rated=raw*c.wireDerating,density=c.jmax*area,allowed=Math.min(c.imax,rated,density),rho=1.7241e-8,factor=1+.00393*(c.copperTemp-20);
 return {d,area,od,reference,raw,rated,density,allowed,ohmsPerM20:rho/(area*1e-6),ohmsPerM:rho*factor/(area*1e-6)};
}
function emWireInfo(c){const a=emWire(c);$('wireInfo').textContent='AWG '+c.awg+': bare diameter '+emFmt(a.d,4)+' mm; copper area '+emFmt(a.area,5)+' mm². AWG reference current '+emFmt(a.reference,4)+' A'+(c.wireManufacturer?' (replaced by your manufacturer coil rating '+emFmt(c.wireManufacturer,4)+' A)':' (700 circular mils/A wiring guideline)')+'. After '+emFmt(c.wireDerating*100,0)+'% derating: '+emFmt(a.rated,4)+' A. Coil density cap: '+emFmt(a.density,4)+' A. Design current cap including supply: '+emFmt(a.allowed,4)+' A. Resistance per metre: '+emFmt(a.ohmsPerM20,5)+' Ω/m at 20 °C; '+emFmt(a.ohmsPerM,5)+' Ω/m at '+c.copperTemp+' °C.';return a;}
function emWinding(c,target){
 const a=emWire(c),{d,od,area}=a,pitch=od*c.packing,perLayer=Math.floor(c.coilLength/pitch),capacityLayers=Math.floor(c.coilThickness/pitch),capacityTurns=perLayer*capacityLayers;
 if(!capacityTurns)throw Error('No turns fit in this winding envelope.');if(capacityTurns>1000000)throw Error('Winding exceeds one million turns; increase wire diameter or reduce the envelope.');
 const key=JSON.stringify([pitch,perLayer,capacityLayers,c.coilID,c.sensorHalf,c.leads]);let sums=EM.windingCache;
 if(!sums||sums.key!==key){
  const mu=4*Math.PI*1e-7,K=new Float64Array(capacityTurns+1),plus=new Float64Array(capacityTurns+1),minus=new Float64Array(capacityTurns+1),wire=new Float64Array(capacityTurns+1);wire[0]=c.leads;
  const zs=Array.from({length:perLayer},(_,k)=>(k-(perLayer-1)/2)*pitch/1000).sort((x,y)=>Math.abs(x)-Math.abs(y)||x-y);let n=0;
  for(let j=0;j<capacityLayers;j++){const r=(c.coilID/2+pitch*(j+.5))/1000;for(const z of zs){n++;const ring=u=>mu*r*r/(2*(r*r+u*u)**1.5);K[n]=K[n-1]+ring(z);plus[n]=plus[n-1]+ring(z-c.sensorHalf/1000);minus[n]=minus[n-1]+ring(z+c.sensorHalf/1000);wire[n]=wire[n-1]+2*Math.PI*r}}
  sums=EM.windingCache={key,K,plus,minus,wire};
 }
 const fieldT=target/1000,at=n=>{const wire=sums.wire[n],R20=wire*a.ohmsPerM20,R=wire*a.ohmsPerM,I=fieldT/sums.K[n],V=I*R;return {wire,R20,R,I,V,supply:V+c.headroom}};
 let low=1,high=capacityTurns;while(low<high){const mid=Math.floor((low+high)/2);if(fieldT/sums.K[mid]<=a.allowed)high=mid;else low=mid+1}const minimumTurnsCurrent=low;
 let N=capacityTurns;
 if(c.turnMode==='minimum'){
  N=minimumTurnsCurrent;let bestV=Infinity,bestN=N;
  for(let n=N;n<=capacityTurns;n++){const w=at(n);if(w.supply<bestV){bestV=w.supply;bestN=n}if(w.I<=a.allowed+1e-12&&w.supply<=c.vmax+1e-10){N=n;break}if(n===capacityTurns)N=bestN}
 }
 const w=at(N),K=sums.K[N],K_mT_A=K*1000,P=w.I*w.I*w.R,J=w.I/area,limit=Math.min(a.allowed,Math.max(0,c.vmax-c.headroom)/w.R),layers=Math.ceil(N/perLayer),lastLayerTurns=N-(layers-1)*perLayer;
 return {target,d,od,pitch,perLayer,layers,capacityLayers,capacityTurns,minimumTurnsCurrent,lastLayerTurns,N,...w,K_mT_A,area,P,J,limit,maxB:K_mT_A*limit,wireReference:a.reference,wireRated:a.rated,wireCurrentCap:a.allowed,currentDensityCap:a.density,ohmsPerM20:a.ohmsPerM20,ohmsPerM:a.ohmsPerM,axialVariation:Math.max(Math.abs(sums.plus[N]/K-1),Math.abs(sums.minus[N]/K-1))*100,outerDiameter:c.coilID+2*layers*pitch,c};
}
function emCoil(){try{
 let r=null;if(EM.result&&$('coilTargetSource').value==='search')r=emFresh();
 if(r){$('manualB').value=r.mag.toFixed(6);$('manualB').disabled=true;$('coilTargetInfo').textContent='Automatically transferred from field search: |B| = '+emFmt(r.mag,6)+' mT'+(r.feasible?' · all modeled operating-range checks passed.':' · separation-only candidate; linear operating-range checks failed.');}
 else{$('manualB').disabled=false;$('coilTargetInfo').textContent=EM.result?'Using your manual field target; it is independent of the resonance search.':'No field-search result yet. The editable manual target is used until a search completes.'}
 const c=emCoilConfig();emWireInfo(c);const b=r?r.mag:c.manualB,w=emWinding(c,b);EM.coils=w;const pass=w.I<=w.limit+1e-12;
 const winding=w.layers===1?w.lastLayerTurns+' turns in one partial/full layer':(w.layers-1)+' full layers × '+w.perLayer+' turns + '+w.lastLayerTurns+' turns in the outer layer';
 $('coilResult').innerHTML='<p class="'+(pass?'ok':'warn')+'"><b>'+(pass?'Coil can reach the target within the entered electrical limits.':'Coil cannot reach the target within the entered limits.')+'</b></p>'+(r&&!r.feasible?'<p class="warn">This coil targets a separation-only candidate. It does not resolve the failed dispersive operating-range check.</p>':'')+'<div class="readouts"><div><small>Calculated winding turns</small><strong>'+w.N+' turns</strong></div><div><small>Required regulated current</small><strong>'+emFmt(w.I,4)+' A</strong></div><div><small>Supply voltage with headroom</small><strong>'+emFmt(w.supply,2)+' V</strong></div></div><p>'+winding+'. '+(c.turnMode==='minimum'?'Smallest turn count found that meets the current and voltage constraints when feasible.':'Full winding envelope used to reduce current.')+' Winding capacity: '+w.capacityTurns+' turns ('+w.capacityLayers+' layers). Actual outer diameter '+emFmt(w.outerDiameter,2)+' mm. Wire including '+c.leads+' m of leads: '+emFmt(w.wire,2)+' m.</p><p>Field/current coefficient '+emFmt(w.K_mT_A,4)+' mT/A; target '+emFmt(b,6)+' mT. R20 = '+emFmt(w.R20,2)+' Ω; R('+c.copperTemp+' °C) = '+emFmt(w.R,2)+' Ω = '+emFmt(w.ohmsPerM,5)+' Ω/m × '+emFmt(w.wire,2)+' m. Coil voltage '+emFmt(w.V,2)+' V; copper dissipation '+emFmt(w.P,2)+' W; current density '+emFmt(w.J,2)+' A/mm².</p><p>Wire/driver design current cap '+emFmt(w.wireCurrentCap,4)+' A; actual allowed current after voltage compliance '+emFmt(w.limit,4)+' A, giving at most '+emFmt(w.maxB,3)+' mT. Centre-to-axial-endpoint field variation over ±'+c.sensorHalf+' mm: '+emFmt(w.axialVariation,3)+'%. Transverse uniformity has not been evaluated.</p><p>Use a regulated constant-current driver. Voltage and power use the entered copper temperature; temperature rise and coil ampacity still need thermal validation.</p>';
 }catch(e){EM.coils=null;$('coilResult').textContent=e.message}}
function emDownload(data,type,name){const url=URL.createObjectURL(new Blob([data],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}
$('findField').onclick=emSearch;$('cancelSearch').onclick=()=>{EM.cancel=true};$('calculateCoil').onclick=emCoil;$('showBranch').onchange=emDrawPattern;
$('applyField').onclick=()=>{try{const r=emFresh();sessionStorage.setItem('NV_EM_APPLY',JSON.stringify({mag:r.mag,theta:Math.acos(r.dir[2])/RAD,phi:(Math.atan2(r.dir[1],r.dir[0])/RAD+360)%360,rotation:r.c.yaw,tiltX:r.c.roll,tiltY:r.c.pitch,cut:r.c.ecut,D:r.c.ed,apar:r.c.eapar,aperp:r.c.eaperp,quad:r.c.eq,gammaN:r.c.egn,temp:r.c.tref,linewidth:r.c.linewidth}));location.href='index.html'}catch(e){$('searchStatus').textContent=e.message}};
$('exportDesign').onclick=()=>{try{const r=emFresh();emCoil();emDownload(JSON.stringify({schema:'MagLab-electromagnet-v1',created:new Date().toISOString(),modelSource:'e354868701d753ce320d2720d3cd58bafafded24',requirements:r.c,field_mT:r.nom.field,magnitude_mT:r.mag,feasible:r.feasible,separation:r.separation,linearCapture:r.capture,linearIsolation:r.linearIsolation,resonances_MHz:r.nom.tracks.map(t=>t.f),A_MHz_per_uT:r.A,pseudoinverse_uT_per_MHz:r.pseudo.P,selected:r.selected.map(i=>i+1),thermal_MHz_per_C:r.kT,reservedDrift_MHz:r.reserved,linearHalfRange_MHz:r.patterns.map(p=>p.used),coil:EM.coils,validationSamples:r.samples,limitations:'Sampled numerical candidate, not a certified global minimum or continuous-domain guarantee; no thermal-rise model.'},null,2),'application/json','NV_Electromagnet_Design.json')}catch(e){$('searchStatus').textContent=e.message}};
$('exportMatrix').onclick=()=>{try{const r=emFresh(),rows=['# A: MHz/uT; pseudoinverse: uT/MHz; frequencies: MHz','# selected: '+r.selected.map(i=>i+1).join(' '),'resonance,frequency_MHz,dBx,dBy,dBz,dT_MHz_C,fieldBound_MHz,thermalBound_MHz,reservedDrift_MHz,linearHalf_MHz'];r.A.forEach((a,i)=>rows.push([i+1,r.nom.tracks[i].f,...a,r.kT[i],r.fieldBound[i],r.thermalBound[i],r.reserved[i],r.patterns[i].used].join(',')));rows.push('# pseudoinverse rows Bx By Bz');r.pseudo.P.forEach(row=>rows.push(row.join(',')));emDownload(rows.join('\n'),'text/csv','NV_Electromagnet_Matrix.csv')}catch(e){$('searchStatus').textContent=e.message}};
try{const saved=JSON.parse(sessionStorage.getItem('NV_EM_SETTINGS')||'null');if(saved){const map={ecut:'cut',yaw:'rotation',pitch:'tiltY',roll:'tiltX',ed:'D',eapar:'apar',eaperp:'aperp',eq:'quad',egn:'gammaN',eth:'theta',eph:'phi',linewidth:'width'};for(const [id,key]of Object.entries(map))if(saved[key]!==undefined)$(id).value=saved[key];$('importStatus').textContent='Imported the diamond orientation, model constants, linewidth and bias direction from the NV simulator. Test-field offsets are centred on zero for this design.'}else $('importStatus').textContent='Using the simulator default diamond and Hamiltonian settings. Open this studio from the simulator to import its current settings.'}catch(e){$('importStatus').textContent='Using default settings.'}
window.addEventListener('resize',()=>{emDrawPattern();emDrawBudget()});

document.querySelectorAll('#requirements input,#requirements select').forEach(el=>el.addEventListener('input',()=>{if(EM.result&&!EM.running)$('searchStatus').textContent='Inputs changed. Plots show the previous design; run the search again before using or exporting it.'}));

let emCoilTimer=0;document.querySelectorAll('#coilSection input,#coilSection select').forEach(el=>{el.addEventListener('input',()=>{clearTimeout(emCoilTimer);emCoilTimer=setTimeout(emCoil,180)});el.addEventListener('change',()=>{clearTimeout(emCoilTimer);emCoil()})});
emCoil();
