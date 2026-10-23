'use strict';
// Independent actual-ideal check. All ring products are obtained by rectangular
// Schur Pieri rules; no companion recurrence, resultant colon or target Hilbert
// formula is used. BigInt fraction-free elimination is exact over Q.
const fs = require('node:fs');
const path = require('node:path');
function abs(a) { return a < 0n ? -a : a; }
function gcd(a,b) { a=abs(a);b=abs(b);while(b){let c=a%b;a=b;b=c;}return a; }
function primitive(v) { let g=0n;for(const a of v)g=gcd(g,a);if(g>1n)for(let i=0;i<v.length;i++)v[i]/=g;let p=v.findIndex(a=>a!==0n);if(p>=0&&v[p]<0n)for(let i=0;i<v.length;i++)v[i]=-v[i];return v; }
class Basis { constructor(n){this.n=n;this.rows=[];} reduce(v){v=v.slice();for(const {p,r} of this.rows){if(v[p]){const a=r[p],b=v[p];for(let i=0;i<v.length;i++)v[i]=a*v[i]-b*r[i];primitive(v);}}return v;} add(v){v=primitive(this.reduce(v));const p=v.findIndex(a=>a!==0n);if(p<0)return null;this.rows.push({p,r:v});this.rows.sort((a,b)=>a.p-b.p);return v;} has(v){return this.reduce(v).every(a=>a===0n);} get rank(){return this.rows.length;} }
function permutations(n){if(!n)return [{p:[],sign:1n}];const out=[];function rec(p,left,sign){if(!left.length){out.push({p,sign});return;}left.forEach((v,i)=>rec(p.concat(v),left.filter((_,j)=>j!==i),i%2?-sign:sign));}rec([],Array.from({length:n},(_,i)=>i),1n);return out;}
function combinations(n,r){const out=[];function rec(a,start){if(a.length===r){out.push(a);return;}for(let i=start;i<=n-r+a.length;i++)rec(a.concat(i),i+1);}rec([],0);return out;}
function ring(k,ell){
 const ps=[];function rec(a,max){if(a.length===k){ps.push(a);return;}for(let t=max;t>=0;t--)rec(a.concat(t),t);}rec([],ell);
 const weights=ps.map(a=>a.reduce((x,y)=>x+y,0));const N=ps.length;
 const unit=Array(N).fill(0n);unit[ps.findIndex(a=>a.every(x=>!x))]=1n;
 const hTargets=[],eTargets=[];
 for(let r=0;r<=ell;r++)hTargets[r]=ps.map((lam,i)=>ps.map((mu,j)=>weights[j]-weights[i]===r&&mu.every((x,t)=>x>=lam[t]&&(t===k-1||mu[t+1]<=lam[t]))?j:-1).filter(j=>j>=0));
 for(let r=1;r<=k;r++)eTargets[r]=ps.map((lam,i)=>ps.map((mu,j)=>weights[j]-weights[i]===r&&mu.every((x,t)=>x>=lam[t]&&x<=lam[t]+1)?j:-1).filter(j=>j>=0));
 const zero=()=>Array(N).fill(0n);
 function applyTargets(v,ts){const out=zero();for(let i=0;i<N;i++)if(v[i])for(const j of ts[i])out[j]+=v[i];return out;}
 function hop(r,v){if(r<0||r>ell)return zero();return applyTargets(v,hTargets[r]);}
 const aCache=new Map();function aop(r,v){if(r<0||r>2*ell)return zero();if(!aCache.has(r)){const columns=ps.map((_,i)=>{const b=zero();b[i]=1n;const out=zero();for(let t=0;t<=r;t++){const w=hop(t,hop(r-t,b));for(let j=0;j<N;j++)out[j]+=w[j];}return out;});aCache.set(r,columns);}const out=zero(),cols=aCache.get(r);for(let i=0;i<N;i++)if(v[i])for(let j=0;j<N;j++)if(cols[i][j])out[j]+=v[i]*cols[i][j];return out;}
 const permCache=new Map();function determinant(indices){const n=indices.length;if(!permCache.has(n))permCache.set(n,permutations(n));const out=zero();for(const {p,sign} of permCache.get(n)){let v=unit;for(let i=0;i<n;i++)v=aop(indices[i][p[i]],v);for(let j=0;j<N;j++)out[j]+=sign*v[j];}return out;}
 function weight(v){let w=null;for(let i=0;i<N;i++)if(v[i]){if(w!==null&&w!==weights[i])throw Error('nonhomogeneous vector');w=weights[i];}return w;}
 const allTransforms=ps.map(lam=>determinant(Array.from({length:k},(_,i)=>Array.from({length:k},(_,j)=>lam[i]-i+j))));
 const fullBasis=new Map();for(const v of allTransforms){const w=weight(v);if(w===null)throw Error('zero full transform');if(!fullBasis.has(w))fullBasis.set(w,new Basis(N));if(!fullBasis.get(w).add(v))throw Error('dependent full transformed rectangular basis');}
 const tails=[];const L=ell-k+1;
 for(let s=1;s<=k;s++){
  const d=s-1,e=ell-k+d,ideal=new Map(),queue=[];
  function addIdeal(v){const w=weight(v);if(w===null)return;if(!ideal.has(w))ideal.set(w,new Basis(N));const b=ideal.get(w).add(v);if(b)queue.push({v:b,w});}
  const minors=combinations(k,s).map(rows=>determinant(rows.map(i=>Array.from({length:s},(_,j)=>L+i+j))));for(const v of minors)addIdeal(v);
  for(let z=0;z<queue.length;z++){const {v,w}=queue[z];for(let r=1;r<=k;r++)if(w+r<=k*ell)addIdeal(applyTargets(v,eTargets[r]));}
  const outside=new Map();let outsideCount=0;for(let i=0;i<N;i++)if(ps[i][d]>e){outsideCount++;const v=allTransforms[i],w=weight(v);if(!ideal.has(w)||!ideal.get(w).has(v))throw Error('outside-hook transform is not in actual ideal: '+JSON.stringify({k,ell,s,lambda:ps[i]}));if(!outside.has(w))outside.set(w,new Basis(N));if(!outside.get(w).add(v))throw Error('dependent outside-hook transform');}
  const idealDim=Array.from(ideal.values()).reduce((a,b)=>a+b.rank,0);
  if(idealDim!==outsideCount)throw Error('actual ideal / outside span mismatch '+JSON.stringify({k,ell,s,idealDim,outsideCount}));
  const byWeight=[];for(let w=0;w<=k*ell;w++){const a=ideal.get(w)?.rank||0,b=outside.get(w)?.rank||0;if(a!==b)throw Error('graded mismatch');if(a)byWeight.push({weight:w,ideal_rank:a,outside_rank:b});}
  if(s===k&&(idealDim!==1||!ideal.has(k*ell)))throw Error('terminal ideal is not the rectangular top line');
  tails.push({s,d,e,maximal_minor_count:minors.length,ideal_dimension:idealDim,outside_count:outsideCount,quotient_dimension:N-idealDim,by_weight:byWeight});
 }
 return {k,ell,ring_dimension:N,full_transformed_basis_rank:Array.from(fullBasis.values()).reduce((a,b)=>a+b.rank,0),tails};
}
const start=Date.now();const cases=[];
for(let k=1;k<=4;k++)for(let ell=k;ell<=Math.min(5,k===1?5:5);ell++)cases.push(ring(k,ell));
const result={method:'rectangular Schur Pieri multiplication, actual Hankel minors, elementary-generator ideal closure, BigInt exact rank; no shifted recurrence or resultant formula',case_count:cases.length,tail_count:cases.reduce((n,c)=>n+c.tails.length,0),runtime_ms:Date.now()-start,passed:true,cases};
const out=path.join(__dirname,'fresh_ordinary_pieri_round4.json');fs.writeFileSync(out,JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify({passed:result.passed,case_count:result.case_count,tail_count:result.tail_count,runtime_ms:result.runtime_ms,output:out}));
