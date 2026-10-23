// Independent exact interpolation proof for excess five.
// No import or read of the parent's symbolic program or result.
const gcd = (a,b) => {a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b];}return a;};
const Q = (n,d=1n) => {n=BigInt(n);d=BigInt(d);if(d<0n){n=-n;d=-d;}const g=gcd(n,d);return {n:n/g,d:d/g};};
const add=(x,y)=>Q(x.n*y.d+y.n*x.d,x.d*y.d);
const neg=x=>Q(-x.n,x.d);
const mul=(x,y)=>Q(x.n*y.n,x.d*y.d);
const div=(x,y)=>Q(x.n*y.d,x.d*y.n);
const trim=p=>{while(p.length>1&&p.at(-1).n===0n)p.pop();return p;};
const pa=(p,q)=>trim(Array.from({length:Math.max(p.length,q.length)},(_,i)=>add(p[i]||Q(0n),q[i]||Q(0n))));
const pm=(p,q)=>{let r=Array(p.length+q.length-1).fill(null).map(()=>Q(0n));for(let i=0;i<p.length;i++)for(let j=0;j<q.length;j++)r[i+j]=add(r[i+j],mul(p[i],q[j]));return trim(r);};
const pe=(p,x)=>p.reduceRight((a,b)=>add(mul(a,Q(x)),b),Q(0n));
const scale=(p,q)=>p.map(x=>mul(x,q));
function choose(n,k){n=BigInt(n);if(k<0||n<BigInt(k))return 0n;let v=1n;for(let i=1;i<=k;i++)v=v*(n-BigInt(i-1))/BigInt(i);return v;}
function coeffGF(n,a){if(a<0||2*a>n)return 0n;let v=0n;for(let j=0;j<=a+1;j++){const rem=n-2*a-j;if(rem<0||rem%3)continue;const q=rem/3;v+=choose(a+1,j)*(q%2?-1n:1n)*choose(q+a,a);}return v;}
function coeffRecurrence(n,a){if(a<0)return 0n;let f=[1n],g=[];for(let i=1;i<=n;i++){let h=Array(Math.floor(i/2)+1).fill(0n);for(let j=0;j<h.length;j++)h[j]=(f[j]||0n)-(g[j]||0n)+(j?(g[j-1]||0n):0n);g=f;f=h;}return f[a]||0n;}
function det(a){a=a.map(r=>r.slice());let old=1n,sign=1n;for(let k=0;k<a.length-1;k++){let pivot=k;while(pivot<a.length&&a[pivot][k]===0n)pivot++;if(pivot===a.length)return 0n;if(pivot!==k){[a[k],a[pivot]]=[a[pivot],a[k]];sign=-sign;}const p=a[k][k];for(let i=k+1;i<a.length;i++){for(let j=k+1;j<a.length;j++){const v=p*a[i][j]-a[i][k]*a[k][j];if(v%old)throw Error('nonexact Bareiss division');a[i][j]=v/old;}a[i][k]=0n;}old=p;}return sign*a.at(-1).at(-1);}
function matrix(r,c,d=5){const k=3*r+c,sign=r%2?-1n:1n;return Array.from({length:d},(_,i)=>{const a=i+1,out=[];for(let j=0;j<=Math.floor((d-1)/2);j++)out.push(sign*coeffGF(k+1,a-j));for(let j=0;j<=Math.floor((d-2)/2);j++)out.push(sign*coeffGF(k+2,a-j));return out;});}
function interpolate(xs,ys){let p=[Q(0n)];for(let i=0;i<xs.length;i++){let b=[Q(1n)],d=1n;for(let j=0;j<xs.length;j++)if(i!==j){b=pm(b,[Q(-xs[j]),Q(1n)]);d*=xs[i]-xs[j];}p=pa(p,scale(b,Q(ys[i],d)));}return trim(p);}
function divideLinear(p,root){if(pe(p,root).n!==0n)return null;const q=Array(p.length-1);q[q.length-1]=p.at(-1);for(let i=q.length-2;i>=0;i--)q[i]=add(p[i+1],mul(Q(root),q[i+1]));return trim(q);}
const serQ=q=>q.d===1n?q.n.toString():q.n.toString()+'/'+q.d.toString();
function reconstruct(){const samples=Array.from({length:12},(_,i)=>BigInt(i+2));const branches=[];for(let c=0;c<3;c++){const ys=samples.map(r=>det(matrix(Number(r),c)));let p=interpolate(samples,ys);const original=p;let factors=[];for(let rr=-12;rr<=12;rr++){let multiplicity=0;while(p.length>1){let q=divideLinear(p,BigInt(rr));if(!q)break;p=q;multiplicity++;}if(multiplicity)factors.push({root:rr,multiplicity});}const checks=[];for(let r=0;r<=16;r++){const a=matrix(r,c);for(let i=0;i<5;i++)for(let j=0;j<5;j++){const blocks=[[1,0],[1,1],[1,2],[2,0],[2,1]];const [off,shift]=blocks[j];const direct=(r%2?-1n:1n)*coeffRecurrence(3*r+c+off,i+1-shift);if(a[i][j]!==direct)throw Error('GF/recurrence disagreement');}const exact=det(a),predicted=pe(original,BigInt(r));if(predicted.d!==1n||predicted.n!==exact)throw Error('determinant disagreement');if(r<2||r>13)checks.push({r,k:3*r+c,det:exact.toString()});}branches.push({c,samples:samples.map((r,i)=>({r:Number(r),det:ys[i].toString()})),degree:original.length-1,coefficients:original.map(serQ),integerLinearFactors:factors,residual:p.map(serQ),checks});}return {method:'exact degree-11 Lagrange interpolation on r=2,...,13; independent GF coefficients and Bareiss determinants; additional recurrence checks',degreeBound:11,branches};}

const formulas=['-r^3(r-1)(r+1)^3(r+2)(3r+2)(3r+4)(3r+5)/5760','r^2(r-1)(r+1)^4(r+2)^2(3r+4)(3r+5)/1920','r^2(r+1)^3(r+2)^3(3r+4)(3r+5)(3r+7)/2880'];
function factored(r,c){r=BigInt(r);let num,den;if(c===0){num=-(r**3n)*(r-1n)*(r+1n)**3n*(r+2n)*(3n*r+2n)*(3n*r+4n)*(3n*r+5n);den=5760n;}else if(c===1){num=r**2n*(r-1n)*(r+1n)**4n*(r+2n)**2n*(3n*r+4n)*(3n*r+5n);den=1920n;}else{num=r**2n*(r+1n)**3n*(r+2n)**3n*(3n*r+4n)*(3n*r+5n)*(3n*r+7n);den=2880n;}if(num%den)throw Error('factored value is nonintegral');return num/den;}
function run(){const out=reconstruct();for(const b of out.branches){for(const s of b.samples)if(factored(s.r,b.c)!==BigInt(s.det))throw Error('factor identity failed');b.factoredFormula=formulas[b.c];b.factorDegree=11;b.factorVerification='degree at most 11 and equality at all twelve exact interpolation nodes';}out.boundaries=[5,6,7].map(k=>{const r=Math.floor(k/3),c=k%3,A=matrix(r,c),D=det(A);return{k,r,c,A:A.map(row=>row.map(String)),normalizedDet:D.toString(),detM:((r%2?-1n:1n)*D).toString()};});return out;}
export {run,coeffGF,coeffRecurrence,det,matrix,interpolate};
import {readFile,writeFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
const result=run();
result.sourceSHA256=createHash('sha256').update(await readFile(new URL(import.meta.url))).digest('hex');
await writeFile(new URL('./odd_excess_five_independent.json',import.meta.url),JSON.stringify(result,null,2)+'\n','utf8');

