// Independent rational checks of the two Lagrangian transition windows.
// Source SHA256: 4b4af4f8b7ed48700074a42625c384b7d6bf383f9503006216f8ce4018a812fc
function qnorm(a,b=1n){if(b<0n){a=-a;b=-b;} var x=a<0n?-a:a,y=b;while(y){var z=x%y;x=y;y=z;}return [a/x,b/x];}

function qa(a,b){return qnorm(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);}

function qm(a,b){return qnorm(a[0]*b[0],a[1]*b[1]);}

function qn(a){return [-a[0],a[1]];}

function qi(a){return qnorm(a[1],a[0]);}

function qd(a,b){return qm(a,qi(b));}

function qstr(a){return a[1]===1n?String(a[0]):`${a[0]}/${a[1]}`;}

function qdet(A){A=A.map(r=>r.map(x=>[...x]));var n=A.length,s=[1n,1n];for(var i=0;i<n;i++){var j=i;while(j<n&&A[j][i][0]===0n)j++;if(j===n)return [0n,1n];if(j!==i){[A[i],A[j]]=[A[j],A[i]];s=qn(s);}var v=A[i][i];s=qm(s,v);for(j=i+1;j<n;j++){var c=qd(A[j][i],v);for(var k=i+1;k<n;k++)A[j][k]=qa(A[j][k],qn(qm(c,A[i][k])));}}return s;}

function transitionMatrix(M,d,g){var cols=[];for(var seed=1;seed<=d;seed++){var u=Array(2*d+1).fill(null).map(()=>[0n,1n]),v=Array(2*d+1).fill(null).map(()=>[0n,1n]);u[seed]=[1n,1n];for(var t=1;t<=(M%2?2*d-1:2*d-2);t++){if(t>d){for(var i=1;i<=d;i++){u[t]=qa(u[t],qn(qm([g[i],1n],u[t-i])));u[t]=qa(u[t],qm([BigInt(i)*g[i],1n],v[t-i]));}}if((M+t)%2===0)v[t]=qd(qn(u[t]),[BigInt(M+t),1n]);}var col=[];if(M%2){for(var t=1;t<=2*d-1;t+=2)col.push(v[t]);}else{for(var t=2;t<=2*d-2;t+=2)col.push(v[t]);var b=u[d];for(var i=1;i<d;i++){b=qa(b,qm([g[i],1n],u[d-i]));b=qa(b,qn(qm([BigInt(i)*g[i],1n],v[d-i])));}col.push(qd(b,[BigInt(M+d),1n]));}cols.push(col);}return Array.from({length:d},(_,i)=>cols.map(c=>c[i]));}

let determinantCount = 0, windowCount = 0, rowCount = 0;
for (let d = 1; d <= 7; d++) {
  let g = [1n];
  for (let root = 1; root <= d; root++) {
    let next = Array(g.length + 1).fill(0n);
    g.forEach((a, i) => { next[i] += a; next[i+1] += BigInt(root)*a; });
    g = next;
  }
  let rho = 1n;
  for (let i = 1; i <= d; i++) for (let j = i+1; j <= d; j++) rho *= BigInt(i+j);
  let W = Array(2*d+1).fill(0n);
  for (let i = 0; i <= d; i++) for (let j = 0; j <= d; j++) W[i+j] += g[i]*g[j]*(j%2 ? -1n : 1n);
  for (let M = 1; M <= 8; M++) {
    let A = transitionMatrix(M,d,g), den = 1n, sign;
    if (M%2) {
      for (let j = 1; j <= d; j++) den *= BigInt(M+2*j-1);
      sign = (d*(d+1)/2)%2 ? -1n : 1n;
    } else {
      den = BigInt(M+d);
      for (let j = 1; j < d; j++) den *= BigInt(M+2*j);
      sign = ((d-1)*(d-2)/2)%2 ? -1n : 1n;
    }
    let actual = qdet(A), expected = qnorm(sign*rho,den);
    if (actual[0]*expected[1] !== expected[0]*actual[1]) throw new Error('Determinant mismatch M='+M+', d='+d);
    determinantCount++;
    let limit = M+2*d;
    let D = Array.from({length:limit+1},()=>[0n,1n]); D[0] = [1n,1n];
    for (let a = 1; 2*a <= limit; a++) {
      let val = [0n,1n];
      for (let i = 1; i <= Math.min(d,a); i++) val = qa(val,qm([BigInt(2*a-i)*W[2*i],2n],D[2*(a-i)]));
      D[2*a] = qd(qn(val),[BigInt(a),1n]);
    }
    let U = Array.from({length:limit+1},()=>[0n,1n]);
    for (let k = 1; k <= limit; k++) {
      let val = [0n,1n];
      for (let i = 1; i <= Math.min(d,k); i++) {
        val = qa(val,qm([BigInt(i)*g[i],1n],D[k-i]));
        val = qa(val,qn(qm([g[i],1n],U[k-i])));
      }
      U[k] = val;
    }
    let seeds = Array.from({length:d},(_,i)=>U[M+i+1]);
    let f = A.map(row=>row.reduce((s,a,i)=>qa(s,qm(a,seeds[i])),[0n,1n]));
    let wanted = M%2 ? Array.from({length:d},(_,j)=>D[M+2*j+1]) : [...Array.from({length:d-1},(_,j)=>D[M+2*j+2]),qm([g[d],1n],D[M])];
    for (let i = 0; i < d; i++) {
      if (f[i][0]*wanted[i][1] !== wanted[i][0]*f[i][1]) throw new Error('Window mismatch M='+M+', d='+d+', row='+i);
      rowCount++;
    }
    windowCount++;
  }
}
console.log(JSON.stringify({determinantCount,windowCount,rowCount,range:'M=1..8, d=1..7',allPass:true}));
