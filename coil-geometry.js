/* MagLab coil former geometry. Coordinates in mm; one connected watertight mesh per segment. */
(function(root){'use strict';
function validate(p){const errors=[];for(const k of ['id','width','build','wall','flange','edge','clearance'])if(!Number.isFinite(p[k]))errors.push(k+' must be a finite number.');if(errors.length)return errors;if(p.id<10||p.id>500)errors.push('Winding inner diameter must be 10–500 mm.');if(p.width<3||p.width>200)errors.push('Winding width must be 3–200 mm.');if(p.wall<1||p.wall>10||2*p.wall>=p.id-4)errors.push('Choose a 1–10 mm wall and leave at least a 4 mm bore.');if(p.build<1||p.build>80)errors.push('Radial winding thickness must be 1–80 mm.');if(p.flange<1||p.flange>10)errors.push('Flange thickness must be 1–10 mm.');if(p.edge<1||p.edge>20)errors.push('Flange height above winding must be 1–20 mm.');if(p.clearance<0||p.clearance>1)errors.push('Total seam clearance must be 0–1 mm.');if(![1,2,4].includes(p.parts))errors.push('Choose one, two or four parts.');return errors}
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],sub=(a,b)=>a.map((v,i)=>v-b[i]);
function dimensions(p){return{bore:p.id-2*p.wall,windingOD:p.id+2*p.build,flangeOD:p.id+2*(p.build+p.edge),length:p.width+2*p.flange}}
function mesh(p,part=0){const e=validate(p);if(e.length)throw Error(e.join(' '));const rb=p.id/2,ri=rb-p.wall,rf=rb+p.build+p.edge,span=rf-rb,keyed=p.keyed&&p.parts>1;
const rr=keyed?[ri,rb,rb+.2*span,rb+.38*span,rb+.68*span,rb+.86*span,rf]:[ri,rb,rf],zz=[0,p.flange,p.flange+p.width,2*p.flange+p.width],N=Math.max(24,Math.round(288/p.parts)),turn=2*Math.PI/p.parts,shift=part*turn,tris=[];
function key(r){if(!keyed)return 0;const u=(r-rb)/span;let f;if(u<=.2||u>=.86)f=0;else if(u<=.38)f=(u-.2)/.18;else if(u<=.68)f=1;else f=(.86-u)/.18;return Math.min(3,span*.18)*f/r}
function point(i,j,k){const r=rr[i],g=p.parts>1?p.clearance/(2*r):0,th=shift+key(r)+g+(turn-2*g)*j/N;return[r*Math.cos(th),r*Math.sin(th),zz[k]]}
function active(i,k){return i>=0&&i<rr.length-1&&k>=0&&k<3&&(k!==1||i===0)}
function quad(q){tris.push([q[0],q[1],q[2]],[q[0],q[2],q[3]])}
for(let i=0;i<rr.length-1;i++)for(let j=0;j<N;j++)for(let k=0;k<3;k++){if(!active(i,k))continue;const v=(a,b,c)=>point(i+a,j+b,k+c),v000=v(0,0,0),v100=v(1,0,0),v010=v(0,1,0),v110=v(1,1,0),v001=v(0,0,1),v101=v(1,0,1),v011=v(0,1,1),v111=v(1,1,1);
if(!active(i-1,k))quad([v000,v001,v011,v010]);if(!active(i+1,k))quad([v100,v110,v111,v101]);if(j===0&&p.parts>1)quad([v000,v100,v101,v001]);if(j===N-1&&p.parts>1)quad([v010,v011,v111,v110]);if(!active(i,k-1))quad([v000,v010,v110,v100]);if(!active(i,k+1))quad([v001,v101,v111,v011]);}
return tris}
function stl(tris){const buffer=new ArrayBuffer(84+50*tris.length),view=new DataView(buffer),label='MagLab coil former | millimetres';for(let i=0;i<label.length;i++)view.setUint8(i,label.charCodeAt(i));view.setUint32(80,tris.length,true);let o=84;for(const tri of tris){let n=cross(sub(tri[1],tri[0]),sub(tri[2],tri[0])),l=Math.hypot(...n);n=n.map(x=>x/l);for(const v of [n,...tri])for(const x of v){view.setFloat32(o,x,true);o+=4}view.setUint16(o,0,true);o+=2}return buffer}
function printMesh(tris){let min=[Infinity,Infinity,Infinity];for(const t of tris)for(const v of t)for(let j=0;j<3;j++)min[j]=Math.min(min[j],v[j]);return tris.map(t=>t.map(v=>v.map((x,j)=>x-min[j])))}
const api={validate,dimensions,mesh,stl,printMesh};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.CoilGeometry=api;
})(typeof globalThis!=='undefined'?globalThis:this);
