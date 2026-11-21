// Shared model snapshot from NV-Simulation index.html, blob e354868701d753ce320d2720d3cd58bafafded24. Frequencies MHz, field mT.
const $=id=>document.getElementById(id), RAD=Math.PI/180, colors=['#6de4c0','#70b6ff','#ffb86b','#d49bff'];
const state={mag:Math.hypot(3.54,1.73,6.95),theta:Math.acos(6.95/Math.hypot(3.54,1.73,6.95))/RAD,phi:Math.atan2(1.73,3.54)/RAD,width:1.3,contrast:2,D:2870,apar:-2.16,aperp:-2.7,quad:-4.95,gammaN:3.077,temp:25,testX:0,testY:0,testZ:0,mwTheta:90,mwPhi:0,rotation:45,tiltY:0,tiltX:0};
const enabled=[true,true,true,true];let lines=[],axes=[],B=[],range=null,view=[],pending=false;
const camera={yaw:0,pitch:0,zoom:1};let interaction='view',b1Handle=null;
let cut='001',mwMode='isotropic';
function thermalPolynomial(celsius){const t=celsius+273.15;return .097*t-.00037*t*t+.00000017*t*t*t}
function effectiveD(){return state.D+thermalPolynomial(state.temp)-thermalPolynomial(25)}
function temperatureSlope(){const t=state.temp+273.15;return 1000*(.097-.00074*t+.00000051*t*t)}
let biasB=[0,0,0];
const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const unit=v=>v.map(x=>x/Math.hypot(...v));
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
function cutBasis(){let n=unit(cut==='111'?[1,1,1]:cut==='110'?[1,1,0]:[0,0,1]),u=unit(cut==='001'?[1,0,0]:[1,-1,0]);return [u,cross(n,u),n]}
function labRotate([x,y,z]){let a=state.tiltX*RAD,b=state.tiltY*RAD,c=state.rotation*RAD;[y,z]=[y*Math.cos(a)-z*Math.sin(a),y*Math.sin(a)+z*Math.cos(a)];[x,z]=[x*Math.cos(b)+z*Math.sin(b),-x*Math.sin(b)+z*Math.cos(b)];return [x*Math.cos(c)-y*Math.sin(c),x*Math.sin(c)+y*Math.cos(c),z]}
function crystalToLab(v){return labRotate(cutBasis().map(row=>dot(row,v)))}
// Real symmetric Jacobi eigensystem, eigenvectors stored as columns.
// Rotating each local transverse axis onto B_perp makes the full S=1, I=1 H real.
function eigensystem(matrix){
 const n=matrix.length,a=matrix.map(r=>r.slice()),v=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>+(i===j)));
 for(let it=0;it<4000;it++){
  let p=0,q=1,max=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(a[i][j])>max){max=Math.abs(a[i][j]);p=i;q=j}
  if(max<1e-10)break;
  if(it===3999)throw Error('Eigensolver did not converge');
  const angle=.5*Math.atan2(2*a[p][q],a[q][q]-a[p][p]),c=Math.cos(angle),s=Math.sin(angle),ap=a[p][p],aq=a[q][q],off=a[p][q];
  for(let k=0;k<n;k++){
   if(k!==p&&k!==q){let x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y}
   let x=v[k][p],y=v[k][q];v[k][p]=c*x-s*y;v[k][q]=s*x+c*y;
  }
  a[p][p]=c*c*ap-2*s*c*off+s*s*aq;a[q][q]=s*s*ap+2*s*c*off+c*c*aq;a[p][q]=a[q][p]=0;
 }
 const order=Array.from({length:n},(_,i)=>i).sort((i,j)=>a[i][i]-a[j][j]);
 return {values:order.map(i=>a[i][i]),vectors:order.map(i=>v.map(row=>row[i]))};
}
const ident=[[1,0,0],[0,1,0],[0,0,1]],spinX=[[0,Math.SQRT1_2,0],[Math.SQRT1_2,0,Math.SQRT1_2],[0,Math.SQRT1_2,0]],spinY=[[0,-Math.SQRT1_2,0],[Math.SQRT1_2,0,-Math.SQRT1_2],[0,Math.SQRT1_2,0]],spinZ=[[1,0,0],[0,0,0],[0,0,-1]];
function kron(a,b){return Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>a[Math.floor(i/3)][Math.floor(j/3)]*b[i%3][j%3]))}
const SX=kron(spinX,ident),SY=kron(spinY,ident),SZ=kron(spinZ,ident),IX=kron(ident,spinX),IZ=kron(ident,spinZ),ZZ=kron(spinZ,spinZ),XX=kron(spinX,spinX),YY=kron(spinY,spinY);
function hamiltonian(bp,bt){return Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>{
 const ms=1-Math.floor(i/3),mi=1-i%3;
 return (i===j?effectiveD()*ms*ms+state.quad*mi*mi:0)+28.025*(bp*SZ[i][j]+bt*SX[i][j])-state.gammaN/1000*(bp*IZ[i][j]+bt*IX[i][j])+state.apar*ZZ[i][j]+state.aperp*(XX[i][j]-YY[i][j]);
}))}
function matrixElement(a,op,b){let val=0;for(let i=0;i<9;i++)for(let j=0;j<9;j++)val+=a[i]*op[i][j]*b[j];return val}
function mwDirection(){const t=state.mwTheta*RAD,p=state.mwPhi*RAD;return [Math.sin(t)*Math.cos(p),Math.sin(t)*Math.sin(p),Math.cos(t)]}
let physicsKey='',drive=[];
function physics(){
 const key=JSON.stringify([state.mag,state.theta,state.phi,state.rotation,state.tiltX,state.tiltY,state.D,state.apar,state.aperp,state.quad,state.gammaN,state.temp,state.testX,state.testY,state.testZ,state.mwTheta,state.mwPhi,cut,mwMode]);
 if(key===physicsKey)return;physicsKey=key;
 const t=state.theta*RAD,p=state.phi*RAD;biasB=[state.mag*Math.sin(t)*Math.cos(p),state.mag*Math.sin(t)*Math.sin(p),state.mag*Math.cos(t)];B=biasB.map((v,i)=>v+[state.testX,state.testY,state.testZ][i]/1000);
 axes=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(v=>crystalToLab(unit(v)));drive=mwDirection();
 lines=axes.map(n=>lineAtField(B,n,drive));
}
function lineAtField(field,n,drive){
  const bp=dot(n,field),transverse=field.map((v,i)=>v-bp*n[i]);let bt=Math.hypot(...transverse),x;
  if(bt>1e-9)x=unit(transverse);else{bt=0;const ref=Math.abs(n[2])<.9?[0,0,1]:[1,0,0];x=unit(cross(ref,n))}
  const y=cross(n,x),mw=[dot(drive,x),dot(drive,y),dot(drive,n)],H=hamiltonian(bp,bt),e=eigensystem(H);
  const brightness=e.vectors.map(v=>v.slice(3,6).reduce((sum,a)=>sum+a*a,0)),transitions=[];
  for(let i=0;i<9;i++)for(let j=i+1;j<9;j++){
   const f=e.values[j]-e.values[i];if(f<1e-6)continue;
   const sx=matrixElement(e.vectors[i],SX,e.vectors[j]),sy=matrixElement(e.vectors[i],SY,e.vectors[j]),sz=matrixElement(e.vectors[i],SZ,e.vectors[j]);
   const strength=mwMode==='isotropic'?(sx*sx+sy*sy+sz*sz)/3:(mw[0]*sx+mw[2]*sz)**2+(mw[1]*sy)**2;
   const brightnessContrast=(brightness[i]-brightness[j])**2,weight=strength*brightnessContrast;
   if(weight>1e-8)transitions.push({f,i,j,strength,brightnessContrast,weight});
  }
  transitions.sort((a,b)=>a.f-b.f);return {bp,bt,transitions,energies:e.values,brightness,eigenvectors:e.vectors,frameX:x,frameY:y};
}


function overlapAssignment(oldLine,newLine){
 const delta=Math.atan2(dot(newLine.frameX,oldLine.frameY),dot(newLine.frameX,oldLine.frameX));
 const score=oldLine.eigenvectors.map(a=>newLine.eigenvectors.map(b=>{let re=0,im=0;for(let k=0;k<9;k++){const phase=(2-Math.floor(k/3)-k%3)*delta,p=a[k]*b[k];re+=p*Math.cos(phase);im+=p*Math.sin(phase)}return re*re+im*im}));
 const n=9,u=Array(10).fill(0),v=Array(10).fill(0),p=Array(10).fill(0),way=Array(10).fill(0);
 for(let i=1;i<=n;i++){p[0]=i;let j0=0;const minv=Array(10).fill(Infinity),used=Array(10).fill(false);
  do{used[j0]=true;const i0=p[j0];let delta=Infinity,j1=0;for(let j=1;j<=n;j++)if(!used[j]){const cur=1-score[i0-1][j-1]-u[i0]-v[j];if(cur<minv[j]){minv[j]=cur;way[j]=j0}if(minv[j]<delta){delta=minv[j];j1=j}}
   for(let j=0;j<=n;j++)if(used[j]){u[p[j]]+=delta;v[j]-=delta}else minv[j]-=delta;j0=j1;
  }while(p[j0]!==0);
  do{const j1=way[j0];p[j0]=p[j1];j0=j1}while(j0!==0);
 }
 const assignment=Array(n);for(let j=1;j<=n;j++)assignment[p[j]-1]=j-1;return assignment;
}

function resonanceCentres(line){
 const z=28.025*line.bp,x=28.025*line.bt/Math.sqrt(2),D=effectiveD();
 const e=eigensystem([[D+z,x,0],[x,0,x],[0,x,D-z]]);
 let anchor=0;for(let i=1;i<3;i++)if(e.vectors[i][1]**2>e.vectors[anchor][1]**2)anchor=i;
 const refs=e.values.map((v,i)=>({i,f:Math.abs(v-e.values[anchor])})).filter(v=>v.i!==anchor).sort((a,b)=>a.f-b.f);
 return refs.map((ref,k)=>{
  const sep=Math.abs(refs[0].f-refs[1].f),window=Math.max(12,4*Math.abs(state.apar)+4*Math.abs(state.aperp));
  let group=line.transitions.filter(tr=>Math.abs(tr.f-ref.f)<=window&&(sep<1e-5||Math.abs(tr.f-ref.f)<=Math.abs(tr.f-refs[1-k].f)+1e-8));
  if(!group.length)group=line.transitions.slice().sort((a,b)=>Math.abs(a.f-ref.f)-Math.abs(b.f-ref.f)).slice(0,3);
  if(!group.length){for(let i=0;i<9;i++)for(let j=i+1;j<9;j++)group.push({i,j,f:Math.abs(line.energies[j]-line.energies[i]),weight:0});return group.sort((a,b)=>Math.abs(a.f-ref.f)-Math.abs(b.f-ref.f))[0]}
  const peak=Math.max(...group.map(tr=>tr.weight));group=group.filter(tr=>tr.weight>=.2*peak);
  return group.sort((a,b)=>Math.abs(a.f-ref.f)-Math.abs(b.f-ref.f)||b.weight-a.weight)[0];
 });
}

function solveInverse3(matrix,rhs){
 const a=matrix.map((row,i)=>[...row,rhs[i]]);for(let p=0;p<3;p++){let pivot=p;for(let i=p+1;i<3;i++)if(Math.abs(a[i][p])>Math.abs(a[pivot][p]))pivot=i;if(Math.abs(a[pivot][p])<1e-18)return null;[a[p],a[pivot]]=[a[pivot],a[p]];const norm=a[p][p];for(let j=p;j<4;j++)a[p][j]/=norm;for(let i=0;i<3;i++)if(i!==p){const c=a[i][p];for(let j=p;j<4;j++)a[i][j]-=c*a[p][j]}}
 return a.map(row=>row[3]);
}
