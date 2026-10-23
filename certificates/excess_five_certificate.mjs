// Exact certificate for the three excess-five determinant identities.
// Node.js, no external modules. All arithmetic below is integral.
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
const dir=path.dirname(fileURLToPath(import.meta.url));
function choose(n,k){if(k<0||n<0||k>n)return 0n;let a=1n;for(let i=1;i<=k;i++)a=a*BigInt(n-k+i)/BigInt(i);return a;}
function signedPow(n){return (n%2)?-1n:1n;}
function eta(r,s,b){if(b<0)return 0n;let z=0n;for(let j=0;j<=b+1;j++){let t=s-2*b-j;if(t%3)continue;t/=3;z+=choose(b+1,j)*signedPow(t)*choose(r+b+t,b);}return z;}
function coefficients(n){let a=[1n],b=[1n];if(n===0)return a;if(n===1)return b;for(let j=2;j<=n;j++){let z=Array(Math.floor(j/2)+1).fill(0n);for(let i=0;i<b.length;i++)z[i]+=b[i];for(let i=0;i<a.length;i++){z[i]-=a[i];z[i+1]+=a[i];}a=b;b=z;}return b;}
function matrixFromEta(r,c){return Array.from({length:5},(_,i)=>{let b=i+1;return [eta(r,c+1,b),eta(r,c+1,b-1),eta(r,c+1,b-2),eta(r,c+2,b),eta(r,c+2,b-1)];});}
function matrixFromRecurrence(r,c){let f=coefficients(3*r+c+1),g=coefficients(3*r+c+2),sg=signedPow(r);let coef=(h,b)=>b<0?0n:(h[b]??0n)*sg;return Array.from({length:5},(_,i)=>{let b=i+1;return [coef(f,b),coef(f,b-1),coef(f,b-2),coef(g,b),coef(g,b-1)];});}
function det(a){a=a.map(row=>row.slice());const n=a.length;let prev=1n,sg=1n;for(let k=0;k<n-1;k++){let p=k;while(p<n&&a[p][k]===0n)p++;if(p===n)return 0n;if(p!==k){[a[p],a[k]]=[a[k],a[p]];sg=-sg;}let pivot=a[k][k];for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++){let num=a[i][j]*pivot-a[i][k]*a[k][j];if(num%prev!==0n)throw new Error('Inexact Bareiss division');a[i][j]=num/prev;}for(let i=k+1;i<n;i++)a[i][k]=0n;prev=pivot;}return sg*a[n-1][n-1];}
function product(r,c){let R=BigInt(r),num,den;if(c===0){num=-(R**3n)*(R-1n)*(R+1n)**3n*(R+2n)*(3n*R+2n)*(3n*R+4n)*(3n*R+5n);den=5760n;}else if(c===1){num=R**2n*(R-1n)*(R+1n)**4n*(R+2n)**2n*(3n*R+4n)*(3n*R+5n);den=1920n;}else{num=R**2n*(R+1n)**3n*(R+2n)**3n*(3n*R+4n)*(3n*R+5n)*(3n*R+7n);den=2880n;}if(num%den!==0n)throw new Error('Inexact product division');return num/den;}
const records=[];for(let r=2;r<=13;r++)for(let c=0;c<3;c++){let a=matrixFromEta(r,c),b=matrixFromRecurrence(r,c);if(a.some((row,i)=>row.some((x,j)=>x!==b[i][j])))throw new Error('Eta/recurrence mismatch');let actual=det(a),expected=product(r,c),residual=actual-expected;if(residual!==0n)throw new Error('Nonzero residual');records.push({r,c,actual:String(actual),expected:String(expected),residual:String(residual)});}
let boundary=det(matrixFromRecurrence(1,2));if(boundary!==42n||boundary!==product(1,2))throw new Error('Boundary failed');
// Newton coefficients at r=2. All twelve residual values are zero,
// hence all finite differences through order eleven are zero.
const residualNewton=Array.from({length:3},(_,c)=>{let q=records.filter(x=>x.c===c).map(x=>BigInt(x.residual)),z=[];while(q.length){z.push(String(q[0]));q=q.slice(1).map((x,i)=>x-q[i]);}return z;});
fs.writeFileSync(path.join(dir,'excess_five_certificate.csv'),'r,c,determinant,product,residual\n'+records.map(x=>[x.r,x.c,x.actual,x.expected,x.residual].join(',')).join('\n')+'\n');
fs.writeFileSync(path.join(dir,'excess_five_certificate.json'),JSON.stringify({method:'Exact recurrence coefficients, independent eta formula, Bareiss determinant; twelve nodes and Newton residuals',nodes:[2,3,4,5,6,7,8,9,10,11,12,13],degree_bound:11,records,residualNewton,boundary_c2_r1:String(boundary),all_pass:true},null,2)+'\n');
console.log(JSON.stringify({records:records.length,residualNewton,boundary_c2_r1:String(boundary),all_pass:true}));
