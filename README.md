# Computational supplement

This package reproduces the finite exact checks reported in the computational appendix of the revised manuscript. Python computations use rational coefficients or exact integers; Node computations use BigInt rational arithmetic. The unrestricted theorems have mathematical proofs in the manuscript. These programs are not proof-assistant certificates.

## Requirements

- Python 3.11 with SymPy 1.14.0 (the recorded Python runs used Python 3.11.3).
- Node.js 18 or later for the separate independent Node suite. No Node packages are required.

Install the Python requirement in your chosen environment with `python -m pip install -r requirements.txt`. Assertions must remain enabled: do not use `python -O`.

## Reproduce

Run from this directory:

```text
python reproduce.py --list
python reproduce.py --suite core --output run_core
python reproduce.py --suite extensions --output run_extensions
python reproduce.py --suite reviewed --output run_reviewed
python reproduce.py --suite near_top --output run_near_top
node reproduce_node.mjs --list
node reproduce_node.mjs --output run_node
```

Every output directory must be new. The runners copy the source files into an isolated work directory, record their SHA256 hashes and the runtime environment, save each task log, and stop on a failed assertion or nonzero exit. The Python `--suite all` selection combines the four Python suites (23 tasks); the Node suite is separate. `--program NAME` selects one Python program in a chosen suite.

## What the suites check

| Suite | Tasks | Scope |
|---|---:|---|
| Python core | 8 | Ordinary Pieri subalgebras and actual Jacobian-minor ideal closure; finite resultant and colon maps; Lagrangian subalgebras and all potentially nonzero transformed strict-partition generators through n=10; counterexample and negative controls. |
| Python extensions | 8 | Finite conductors with endpoint controls, intrinsic coefficient-algebra obstructions, the infinite two-row failure family, two-row candidates through width 30, and Catalan moment identities. |
| Python reviewed | 5 | Coefficient criterion at 527 parameter pairs; symbolic Lagrangian transition identities; 661 actual top intersection values; the rational versus Q(sqrt(3)) algebra examples; an independent complete Q-flag and derivative implementation through n=5. |
| Python near_top | 2 | Actual defining-ideal ranks in 113 cases in the six highest weights for widths 2 through 20; an independent symbolic derivation of the two near-top 3-by-3 determinants, discriminants and square gaps. |
| Independent Node | 3 | Actual ordinary Pieri ideal closure in 14 rings and 30 tails; 56 Lagrangian determinants and 224 coefficient rows; three degree-bounded interpolation certificates for excess five. |

The two older ordinary Pieri programs share their multiplication implementation. The two older Lagrangian ring programs share their squarefree reducer. The new reviewed and Node programs supply separate implementations where stated; programs in one family are not advertised as fully independent of each other.

## Recorded executions

The records directory contains the actual October 14, 2026 isolated executions, including manifests, source hashes, logs and outputs. The 16 core-plus-extension Python tasks, five reviewed Python tasks, two near-top Python tasks, and three Node tasks were executed as separate batches. Their wall-clock intervals may overlap; they are not added to the manuscript revision duration. 

The excess-five interpolation check proves the stated finite polynomial identities using a degree bound of 11 and twelve exact nodes. It does not classify all odd excesses. The uniform defect bound in the manuscript is proved analytically by positivity of a Gram submatrix. The classification of the six highest weights is also proved for every width: its two exceptional pairs come from a moment recurrence and a rational-square obstruction, rather than from the finite rank table. The higher odd-excess exploration is not part of this supplement or of the paper's claims.

## Version 15 October 2026

This version accompanies the reorganized manuscript. The two-row identity
certificate is in Appendix B, and the computational description is in
Appendix C. The mathematical results tested by the older suites are retained;
their original October 14 records and source hashes are preserved as historical
execution records, not relabelled as new executions.

The referee_20261015 directory adds two independent checks executed against
the revised manuscript:

- two_row_symbolic_referee.mjs uses only BigInt rational and sparse polynomial
  arithmetic to prove six excess-three/five determinant identities and four
  near-top determinant/discriminant identities exactly.
- ordinary_macaulay_referee.py uses direct weighted Macaulay ideal multiples,
  actual Jacobian minors, rational kernels and multiplication constraints in
  the cases (k,ell,m)=(2,2,1),(2,3,2),(2,2,2). It includes derivative source
  weights above the ambient top. Its optional first argument supplies only
  the manuscript hash and is not used to calculate the ideals.

To reproduce these checks, run from the referee_20261015 directory:

    node run_symbolic_referee.mjs /path/to/grassmannian_frontiers_revised.tex

    python ordinary_macaulay_referee.py /path/to/grassmannian_frontiers_revised.tex

The source argument may be omitted; the resulting record then has a null
manuscript hash. Node.js requires no packages. The Python check requires
SymPy 1.14.0, as specified above. Keep assertions enabled.

The runners write JSON records beside their source files. The distributed
records are the actual successful executions; source hashes identify their
versions. These finite and exact symbolic checks supplement the mathematical
proofs and do not certify a general classification of all odd excesses.

## Minor revision on 15 October 2026

The abstract and contribution discussion have been revised, the core
maximum-hook bounds have been stated explicitly, and the ordinary proof
completion has been shortened by applying the common rigidity theorem.
The coefficient-algebra calculation formerly in Example 10.8 is now in
Appendix D; its statement remains in Example 10.8. The mathematical
statements tested by this package are unchanged.

The file `minor_revision_provenance.json` links the revised manuscript to
the previously reviewed source. All distributed execution records keep
their original dates and hashes; none has been relabelled as a new run.

## Polish after the strict review on 15 October 2026

Appendix B.4 retains the explicit coefficient formulas, degree bound of
eleven, and the twelve-node zero-residual certificate for the three
excess-five identities. The large integer table has moved to
`certificates/excess_five_certificate.csv`. The same directory contains
the executable certificate, its exact JSON results and a README.

To reproduce this certificate from the package root, run:

    node certificates/excess_five_certificate.mjs

No external modules are needed. It independently compares the coefficient
recurrence with the eta formula, computes the five-by-five determinants
by exact BigInt Bareiss elimination, verifies all 36 determinant/product
values, reconstructs the zero differences in Newton form, and checks
the boundary value 42. All required divisions are exact. The runner
writes the CSV and JSON beside itself. The polynomial degree argument
and legal parameter ranges remain in the manuscript.

`polish_provenance.json` records the current manuscript hash and the new
execution. `minor_revision_provenance.json` remains as a historical record
of the earlier revision. Older program execution dates and source hashes
have not been changed or advertised as new runs.
