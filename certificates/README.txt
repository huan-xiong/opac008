Excess-five determinant certificate
==================================
Run from this directory:
  node excess_five_certificate.mjs

Requirements: Node.js supporting BigInt and standard ES modules; no third-party package.
Outputs: excess_five_certificate.csv and excess_five_certificate.json.

The three determinants are built twice:
(1) from the manuscript's explicit finite eta formula;
(2) independently from the recurrence h_n=u*h_(n-1)-(u^2-v)*h_(n-2), with normalization (-1)^r.
The matrices must agree entry by entry. Bareiss elimination uses BigInt arithmetic and verifies every division is exact. At each of the twelve nodes r=2,...,13 it compares each determinant to the displayed factored product. The CSV records all 36 exact comparisons. The JSON additionally records the zero Newton residual coefficient vectors and the boundary r=1,c=2, determinant 42.

Proof scope: the twelve nodes suffice only because the manuscript proves that each determinant and product has degree at most eleven. Together they establish polynomial identities for every r. The legal widths remain r>=2 for c=0,1 and r>=1 for c=2. Finite checks alone do not establish the other results in the paper.

This runner and data were regenerated in the October 15 polishing pass. No previous PASS, hand-copied table, or floating-point arithmetic is used.
