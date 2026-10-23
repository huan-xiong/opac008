// Independent symbolic verification for the revised two-row appendix.
// No external packages, previous checker code, or old test outputs are used.
const ratTR = (n,d=1n)=>{n=BigInt(n);d=BigInt(d);if(d<0n){n=-n;d=-d;}var a=n<0n?-n:n,b=d;while(b){var z=a%b;a=b;b=z;}return[n/a,d/a];};
const raddTR = (a,b)=>ratTR(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const rmulTR = (a,b)=>ratTR(a[0]*b[0],a[1]*b[1]);
const ptrimTR = a=>{while(a.length>1&&a.at(-1)[0]===0n)a.pop();return a;};
const paddTR = (a,b)=>ptrimTR(Array.from({length:Math.max(a.length,b.length)},(_,i)=>raddTR(a[i]||[0n,1n],b[i]||[0n,1n])));
const pmulTR = (a,b)=>{var c=Array.from({length:a.length+b.length-1},()=>[0n,1n]);for(var i=0;i<a.length;i++)for(var j=0;j<b.length;j++)c[i+j]=raddTR(c[i+j],rmulTR(a[i],b[j]));return ptrimTR(c);};
const pscaleTR = (a,c)=>a.map(v=>rmulTR(v,c));
const factTR = b=>{var f=1n;for(var j=2;j<=b;j++)f*=BigInt(j);return f;};
const binTR = (n,b)=>factTR(n)/(factTR(b)*factTR(n-b));
const etaTR = (s,b)=>{if(b<0)return[[0n,1n]];var out=[[0n,1n]];for(var j=0;j<=b+1;j++){var w=s-2*b-j;if(w%3)continue;var v=w/3,p=[[1n,1n]];for(var h=0;h<b;h++)p=pmulTR(p,[ratTR(b+v-h),[1n,1n]]);var a=binTR(b+1,j)*(Math.abs(v)%2?-1n:1n);out=paddTR(out,pscaleTR(p,ratTR(a,factTR(b))));}return out;};
const detPolyTR = A=>{var n=A.length,out=[[0n,1n]];function rec(p,left){if(!left.length){var t=[[1n,1n]],inv=0;for(var i=0;i<n;i++){t=pmulTR(t,A[i][p[i]]);for(var j=i+1;j<n;j++)if(p[i]>p[j])inv++;}out=paddTR(out,pscaleTR(t,ratTR(inv%2?-1n:1n)));return;}for(var j=0;j<left.length;j++)rec([...p,left[j]],left.filter((_,i)=>i!==j));}rec([],Array.from({length:n},(_,i)=>i));return out;};
const matPolyTR = (c,d)=>{var cols=[];for(var j=0;j<=Math.floor((d-1)/2);j++)cols.push([c+1,j]);for(var j=0;j<=Math.floor((d-2)/2);j++)cols.push([c+2,j]);return Array.from({length:d},(_,i)=>cols.map(([s,j])=>etaTR(s,i+1-j)));};
const prodTR = (f,den=1n,sgn=1n)=>{var a=[[sgn,den]];for(var [c,b,n] of f)for(var j=0;j<n;j++)a=pmulTR(a,[ratTR(c),ratTR(b)]);return a;};
const eqPolyTR = (a,b)=>paddTR(a,pscaleTR(b,[-1n,1n])).every(v=>v[0]===0n);

// Sparse exact Z[m,t] polynomials. The near-top comparison clears
// only the positive denominators D1^2 D2^2 D3 from the Pearson recurrence.
const bC = n => new Map(BigInt(n) === 0n ? [] : [['0,0', BigInt(n)]]);
const bM = new Map([['1,0',1n]]), bT = new Map([['0,1',1n]]);
const bAdd = (...polys) => { const r = new Map(); for (const p of polys) for (const [k,v] of p) { const a=(r.get(k)||0n)+v; if(a) r.set(k,a); else r.delete(k); } return r; };
const bScale = (p,n) => new Map([...p].map(([k,v])=>[k,v*BigInt(n)]).filter(([,v])=>v!==0n));
const bMul = (...polys) => polys.reduce((a,b)=>{ const r=new Map();for(const [ak,av] of a)for(const [bk,bv] of b){const [am,at]=ak.split(',').map(Number),[bm,bt]=bk.split(',').map(Number),key=(am+bm)+','+(at+bt),v=(r.get(key)||0n)+av*bv;if(v)r.set(key,v);else r.delete(key);}return r;},bC(1));
const bLin = (a,b) => bAdd(bScale(bM,a),bC(b));
const bNeg = p => bScale(p,-1);
const bCoeffT = (p,j) => new Map([...p].filter(([k])=>Number(k.split(',')[1])===j).map(([k,v])=>[k.split(',')[0]+',0',v]));
function nearTopAudit(sigma) {
 const D1=bLin(1,2+sigma),D2=bLin(1,3+sigma),D3=bLin(1,4+sigma);
 const A=bAdd(bLin(2,3*sigma),bMul(bScale(bM,3),bT));
 const B=bAdd(bMul(bLin(2,2+3*sigma),A),bMul(bScale(bLin(1,1),3),D1));
 const C=bAdd(bMul(bLin(2,4+3*sigma),B),bMul(bScale(bLin(1,2),3),A,D2));
 const P=bScale(bLin(1,1+sigma),3);
 const R=bAdd(bScale(bMul(D2,D3),9),bScale(D3,-12),bScale(D2,-24),bC(36));
 const determinantNumerator=bAdd(bMul(A,C,D2),bNeg(bMul(B,B,D3)),bNeg(bMul(P,C,D1)),bMul(P,A,B,D3),bMul(R,B,D1),bNeg(bMul(R,A,A,D2)));
 let Q, discExpected;
 if(sigma===0){
  Q=bAdd(bMul(bScale(bM,9),bLin(3,7),bT,bT),bMul(bScale(bAdd(bScale(bMul(bM,bM),6),bScale(bM,13),bC(-6)),3),bT),bNeg(bAdd(bScale(bMul(bM,bM),9),bScale(bM,38),bC(48))));
  discExpected=bMul(bC(81),bLin(1,2),bLin(1,2),bAdd(bScale(bMul(bM,bM),16),bScale(bM,32),bC(1)));
 } else {
  Q=bAdd(bMul(bScale(bM,9),bLin(2,7),bT,bT),bMul(bScale(bAdd(bScale(bMul(bM,bM),2),bScale(bM,9),bC(6)),6),bT),bNeg(bAdd(bScale(bMul(bM,bM),6),bScale(bM,37),bC(60))));
  discExpected=bMul(bC(144),bLin(1,3),bLin(1,3),bAdd(bScale(bMul(bM,bM),4),bScale(bM,14),bC(1)));
 }
 const expected=bMul(bC(sigma===0?2:3),bM,bLin(1,1),Q);
 const a=bCoeffT(Q,2),b=bCoeffT(Q,1),c=bCoeffT(Q,0);
 const disc=bAdd(bMul(b,b),bScale(bMul(a,c),-4));
 return {sigma, determinantResidualTerms:bAdd(determinantNumerator,bNeg(expected)).size, discriminantResidualTerms:bAdd(disc,bNeg(discExpected)).size};
}
export function runAudit(){
 const three=[prodTR([[0,1,2],[1,1,2],[2,3,1]],4n),prodTR([[0,1,1],[1,1,3],[4,3,1]],4n),[[0n,1n]]];
 const five=[prodTR([[0,1,3],[-1,1,1],[1,1,3],[2,1,1],[2,3,1],[4,3,1],[5,3,1]],5760n,-1n),prodTR([[0,1,2],[-1,1,1],[1,1,4],[2,1,2],[4,3,1],[5,3,1]],1920n),prodTR([[0,1,2],[1,1,3],[2,1,3],[4,3,1],[5,3,1],[7,3,1]],2880n)];
 const excessThree=[0,1,2].map(c=>({residue:c,identity:eqPolyTR(detPolyTR(matPolyTR(c,3)),three[c])}));
 const excessFive=[0,1,2].map(c=>{const d=detPolyTR(matPolyTR(c,5));return{residue:c,identity:eqPolyTR(d,five[c]),degree:d.length-1};});
 const nearTop=[0,1].map(nearTopAudit);
 return {arithmetic:'BigInt exact rational Q[r] and sparse Z[m,t]',excessThree,excessFive,nearTop,pass:excessThree.every(v=>v.identity)&&excessFive.every(v=>v.identity)&&nearTop.every(v=>v.determinantResidualTerms===0&&v.discriminantResidualTerms===0)};
}
