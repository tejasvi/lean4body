# lean4body

Programs, certificate files, logs and a Lean 4 formalization for the paper [*Uniqueness of convex
central configurations in the planar four-body problem*](https://doi.org/10.5281/zenodo.23004307)
by Tejasvi S. Tomar. This release is archived at
[doi:10.5281/zenodo.23004297](https://doi.org/10.5281/zenodo.23004297).

The paper proves:

> **Theorem A.** Let `m1, ..., m4 > 0`. For each cyclic order of the four bodies there is exactly
> one convex central configuration with that cyclic order, up to similarity. It is a nondegenerate
> local minimum of `U I^(1/2)` on the shape space, and it depends real-analytically on the masses.

This answers Problem 10 in the list of Albouy, Cabral and Santos (2012). The main step is
Theorem B, the bound `Q ≥ K/4` for the Hessian at every convex central configuration. Its proof
comes down to Proposition 5.1: at every normalized convex central configuration there is a vector
`y` with `tr S(y) < 3/4`. Proposition 5.1 is proved by interval arithmetic on the
three-dimensional set of these configurations. The numbers of sections, propositions, tables and
checks below refer to the paper.

## Contents

| Path | Contents |
|---|---|
| [`work/rig`](work/rig) | the program that proves Proposition 5.1, its certificate files and the logs of its runs (Table 1); the scripts of checks (2), (3) and (4) |
| [`work/cc4.py`](work/cc4.py) | an independent implementation of the geometry in `mpmath` (positions, areas, masses, and `G` by the cotangent formula), used by checks (2) and (3) |
| [`paper`](paper) | computer-algebra checks of identities used in the paper, and the script of Figure 1 |
| [`paper/checks`](paper/checks) | check (1); the comparison of the chart formulas with the direct formulas in check (3); the scripts of Proposition 7.2 and of the remarks after Conjecture 7.3 |
| [`work/indep`](work/indep) | the second implementation (Section 5.6, Table 2), the specification it was written from (`SPEC.md`), the record of its development (`README.md`) and its logs |
| [`C4`](C4) | the Lean formalization (Section 6.2); see [`C4/README.md`](C4/README.md) |

The files of the proof of Proposition 5.1 and of the checks of Section 5.5:

| File | Role |
|---|---|
| `work/rig/bb.py` | the procedure of Section 5.3 |
| `work/rig/cert.py` | the functions of Section 5.1, in the direct coordinates (`direct`, `direct_P`) and in the blow-up chart with the formulas of the Appendix (`chart`, `chart_P`) |
| `work/rig/jet.py` | the jet evaluation and the mean-value enclosures of Section 5.2 |
| `work/rig/cert_direct_075.txt`, `cert_chart_075.txt` | the certificate files |
| `paper/checks/verify_cert.py` | check (1): rebuilds every leaf of a certificate file and verifies it again |
| `work/rig/check_cert.py` | check (2): the structure of a certificate file, and sampling |
| `work/rig/check_identity.py` | check (3): the identity `F ∏ s_e = -P` at random points, and the eigenvalues at random convex central configurations |
| `paper/checks/chart_identity.py`, `chart_identity_psi.py` | check (3): the chart formulas against the direct formulas, in computer algebra |
| `work/rig/sample_sup.py`, `work/indep/maxsearch.py` | check (4): `tr G` at random points, and its local maximization (floating point) |
| `paper/sym_checks.py` | computer-algebra checks of identities used in the proofs: `F ∏ s_e = -P`, the formulas for `w_e`, the factorizations and signed areas in the coordinates of Corbera, Cors and Roberts, the derivative of Heron's formula, and the chart formulas `h_i` |
| `paper/fig_lmax.py` | Figure 1 |
| `paper/checks/corner_G0.py` | the computer algebra in the proof of Proposition 7.2 |
| `paper/checks/corner_lmax.py`, `lmax_locate.py`, `lmax_boundary.py` | the floating-point observations after Conjecture 7.3 |

## Software

The computations use Python 3.14.4, `python-flint` 0.9.0 (the ball arithmetic library Arb),
`mpmath` 1.3.0, `sympy` 1.14.0 and `numpy` 2.5.3. Figure 1 also needs `matplotlib` 3.11.2, and
`lmax_boundary.py` needs `scipy` 1.18.1.

```bash
pip install python-flint==0.9.0 mpmath==1.3.0 sympy==1.14.0 numpy==2.5.3 matplotlib==3.11.2 scipy==1.18.1
```

The formalization uses Lean 4.35.0-rc3 (`C4/lean-toolchain`) and Mathlib at the tag `v4.35.0-rc3`
(`C4/lake-manifest.json`).

## Reproducing the results

Each block of commands below starts in the root of the repository. The comments give the running
times on a machine with 32 hardware threads. The logs and the certificate files in the repository
were produced by running these commands on the files of this commit, with the versions above. A
rerun reproduces them except for running times, dates and the order of progress lines. The
certificate files should be compared after sorting.

### The computation (Proposition 5.1)

```bash
cd work/rig
THRESH=0.75 python bb.py direct 40 28 56 8 1e-5 32 cert_direct_075.txt > run_direct_075.log   # 4 min
THRESH=0.75 python bb.py chart 8 96 96 8 1e-6 32 cert_chart_075.txt > run_chart_075.log       # 10 s
```

The arguments are the region (`direct`, or `chart` for the blow-up chart), the numbers of initial
boxes in the three coordinate directions, the reciprocal of the value `b = 1/8` that separates the
two regions, the width below which a box that would have to be bisected counts as a failure, the
number of processes and the output file. `THRESH` is the threshold `3/4`.

In a certificate file, a header line with the parameters is followed, for each initial box, by a
line `B` with its bounds and one line for each of its leaves: `D path` for a discarded leaf, and
`A path run ...` for an accepted leaf. The path is the sequence of bisections from the initial box
(`-` if there are none), and each run is written as `xlo,xhi,y1,y2,t`: the ends of the run of the
dependent coordinate, the vector `y` (`Y` in the chart) and the bound `t_{B,X}`. The initial boxes
appear in the order in which the processes finish them, so a rerun writes the same lines in another
order. Sorted, the two files have these hashes:

```bash
LC_ALL=C sort cert_direct_075.txt | sha256sum   # c4c119fbf15945cefd4c93c15e73c118d01e15f173caef45e657e92d2550830e
LC_ALL=C sort cert_chart_075.txt | sha256sum    # 7d4768f24f851abb2f561f3cfa0150ce2d623285b608aaceebe2e28ce0feea1d
```

### The checks (Section 5.5)

The checks are not part of the proof.

```bash
cd paper/checks
python verify_cert.py ../../work/rig/cert_direct_075.txt 32 > verify_direct.log   # check (1); 2 min
python verify_cert.py ../../work/rig/cert_chart_075.txt 32 > verify_chart.log     # 7 s
cd ../../work/rig
python check_cert.py cert_direct_075.txt 1500 4 > check_cert_direct.log   # check (2), 1500 points, seed 4; 5 min
python check_cert.py cert_chart_075.txt 1500 3 > check_cert_chart.log     # seed 3; 6 min
python check_identity.py > check_identity.log                             # check (3); 2 s
python sample_sup.py direct 1 > sample_sup_direct.log                     # check (4), 200 000 points, seed 1; 19 s
python sample_sup.py chart 2 > sample_sup_chart.log                       # seed 2; 34 s
cd ../../paper/checks
python chart_identity.py > chart_identity.log                             # check (3); under 1 s
python chart_identity_psi.py > chart_identity_psi.log                     # check (3); 8 min
cd ..
python sym_checks.py > sym_checks.log                                     # under 1 s
python fig_lmax.py > fig_lmax.log                                         # Figure 1; 20 s
```

`chart_identity.py` compares `T_e` and `β_e` for the six edges and the functions `g_i`.
`chart_identity_psi.py` compares `D_e` for the six edges and `P̂`, with four of the five values of
`ψ` that occur replaced by independent symbols; an identity that holds for all values of these
symbols holds in particular for the true ones. It runs the seven comparisons in parallel. The
local maximization of check (4) is `maxsearch.py` in the second implementation, below.
`fig_lmax.py` writes `fig_lmax.pdf` and `fig_lmax.png`.

### Proposition 7.2 and Conjecture 7.3

```bash
cd paper/checks
python corner_G0.py > corner_G0.log           # 1 s
python corner_lmax.py > corner_lmax.log       # 6 s
python lmax_locate.py > lmax_locate.log       # 30 s
python lmax_boundary.py > lmax_boundary.log   # 66 s
```

`corner_G0.py` does the computer algebra in the proof of Proposition 7.2 and compares its limits
with the chart formulas at a small `b`. The other three scripts compute the floating-point values
quoted after Conjecture 7.3, which are not part of any proof: `corner_lmax.py` maximizes
`λ_max(G)` over `γ` near the corner, `lmax_locate.py` samples `λ_max(G)` in the two regions, and
`lmax_boundary.py` maximizes it locally at fixed `b`.

### The second implementation (Section 5.6)

```bash
cd work/indep
python checks.py 300 > checks.log                    # the checks done before the branch and bound; 3 s
bash run_all.sh                                      # threshold 3/4, then the sampling test of the leaves; 3 min
bash run_tight.sh                                    # threshold 19/32, then the sampling test; 6 min
python negative_control.py > negative_control.log    # thresholds below sup tr G, which must fail; 24 s
python explore.py direct 60 > explore_direct.log     # under 1 s
python explore.py chart 60 > explore_chart.log       # under 1 s
python maxsearch.py direct > maxsearch_direct.log    # 20 s
python maxsearch.py chart > maxsearch_chart.log      # 16 s
```

The two shell scripts run the interpreter named in `PY` (default `python`). They write the logs
and `.json` files of the runs in Table 2 and of the sampling tests, and the accepted leaves in
`.npy` files, which are not included in the repository (the largest has 129 MB).

`README.md` in that directory is the record of the development: the formulas, the method, the
checks and the results. `SPEC.md` is the specification the program was written from. In its paths,
`~/convex4` stands for the root of this repository and `~/sundman/venv/bin/python` for a Python
interpreter with the versions listed above.

### The Lean formalization (Section 6.2)

```bash
cd C4
lake exe cache get                                                       # fetch the Mathlib build cache
scripts/build_certs.sh > logs/build.log 2>&1                             # 21 min
lake env lean Axioms.lean > logs/axioms.log                              # 6 s
scripts/nanoda_check.sh nanoda-out 16 > logs/nanoda.log 2>&1             # 29 min
python3 scripts/nanoda_control.py nanoda-out > logs/nanoda_control.log   # 84 s
```

`build_certs.sh` builds everything, the certificates in batches that keep the memory well below
128 GB. `Axioms.lean` prints the axioms of the main statements, `nanoda_check.sh` checks the whole
proof again with the type checker nanoda, and `nanoda_control.py` checks that nanoda rejects a
certificate in which one number is changed. The 3418 certificates of the computation are in
`C4/C4Cert`, and `C4/Gen.lean` is the unverified program that writes their hints.
[`C4/README.md`](C4/README.md) lists the formal statements, describes the trust base and the
check with nanoda, and gives the commands that regenerate the certificates.

## SHA-256

The hashes and sizes of the programs outside `C4`, the certificate files and the specification.
The commit of the release identifies all files.

| File | Bytes | SHA-256 |
|---|--:|---|
| `work/rig/bb.py` | 10093 | `b936e84984bb652ab8906d54a93a2bd8a667eee61bd7290a8868394d47cca8df` |
| `work/rig/cert.py` | 11260 | `a93f9174d9d70a544d2628d67c9a13094f48788603bfd08199d178173d160984` |
| `work/rig/cert_chart_075.txt` | 5613606 | `f0a3cc871749280df6b62801125a9554c4db092d63c2c261fc0fde98edb6f2d8` |
| `work/rig/cert_direct_075.txt` | 39500202 | `187a4eebba4895a51c065eb89be4b5c17c6c84be3f72dffd5a8bf6654cfeddbe` |
| `work/rig/check_cert.py` | 8419 | `626b1eaeb9655e37cefcf63050130c1993911315731d0ce25ed915f8a61a1858` |
| `work/rig/check_identity.py` | 3483 | `e4e131078792c0c3492adf6371010d2fd884d6c0afaa785268c666a6958b7cdc` |
| `work/rig/jet.py` | 2144 | `5017362da35f3b96514f131521e66c0e78f2bdce1244195d603d9d529838595b` |
| `work/rig/sample_sup.py` | 2397 | `75fdde3da680210327fae32ca38328b0992761797007227fee2b9748cb88c094` |
| `work/cc4.py` | 5911 | `aa10b4b662a3d54f73609f87000fed73994f76e2ea4846ca2498352e23cca11b` |
| `paper/fig_lmax.py` | 4023 | `9db3b57aeff163ed2d52ca98e8d25223ca10d51c29f899f04aae74b8db4ae10c` |
| `paper/sym_checks.py` | 4759 | `455881031584946b383ffd7f1811a427a283447a2fd4942d8f1af8e4e70219e6` |
| `paper/checks/chart_identity.py` | 2258 | `e6ef5dc813778b726848f11451a448e0a44fe9b2612e8f075a3ce8c17f88eb12` |
| `paper/checks/chart_identity_psi.py` | 2608 | `a78becb5c9e3ab82add7af1066f28a4dc4167c7bc8c56ad0dc23d3d558e1db51` |
| `paper/checks/corner_G0.py` | 4094 | `d62df029aa85424f2e20a652107fb6824e8ab0729c262ff5d0e0589cd34e1c58` |
| `paper/checks/corner_lmax.py` | 2416 | `37e3365b7b73026c59dac68d66fa227c1e8ccbdd48e1bc49a34d0d74720f357c` |
| `paper/checks/lmax_boundary.py` | 3243 | `84cf50b84a3a258e0647573570f7463c8a7fe9b7ac98f0303a49ca7e0ce55f34` |
| `paper/checks/lmax_locate.py` | 2787 | `36bf9a27099855bb54258b90ff6cd712cb7c2f6cb4d0da0c10e53e6bfcc271fe` |
| `paper/checks/verify_cert.py` | 4254 | `3cd2cac76d2ea47f1984c432070903f1f64e6294f84a6da00a34b8d72d092d71` |
| `work/indep/SPEC.md` | 6582 | `0e4095a1579c209f18b109cf188669e517161db02b0c3a1dd73042d6a0695fe5` |
| `work/indep/ad.py` | 2988 | `0f4a5f1c6fd9cf3b71b1bc13bac5ddb979e4ab54deab6e2f560d040298326a07` |
| `work/indep/check_leaves.py` | 7249 | `54c7e5e29c96c3dd228b4188005da6d9fc6a1c217bc2d1d3030d3e3299df5b2d` |
| `work/indep/checks.py` | 14371 | `cadad72f2b7dfdd265adb46522dff5548008490f2f514f57f3cbdee8e3c88b9b` |
| `work/indep/core.py` | 12918 | `77b2e4d2feee7018d8a4f938acf5dbeabd2c1c7722f9def3568de456dfca43fe` |
| `work/indep/explore.py` | 5444 | `efeef2ccc72bb56eb409f5579ad82130ee9e87cb1dd2b217fece3c2e13fbb968` |
| `work/indep/formulas.py` | 12216 | `ec2fbec7fb3fbeef9b2b84ea6c1b8cbb50983737d1a79dd809c67992868aed93` |
| `work/indep/maxsearch.py` | 3854 | `f0d08ceef70e8453debbfc5a3502e13332fbb6c4e30cc3286d43c172b9b0ff8c` |
| `work/indep/negative_control.py` | 2959 | `34b8d33a7a99cbefcfd32f948a068a5f0eda1345c10754764f7ccc1dce5e67ef` |
| `work/indep/run_all.sh` | 738 | `2c32e75552260e0a2b8a8080df7c473de6b44a0c728e18787bb70cde06929837` |
| `work/indep/run_tight.sh` | 874 | `0bf041daf2484624f76605b9b4c34c517dfeb243aa60ab81f3248868b21090fe` |
| `work/indep/verify.py` | 10190 | `25193cea9f654787c82d91cf7e4deeafb3c82801fb520632de197a6286cffdf5` |

## Use of AI tools

Generative AI tools were used to write the programs and the Lean formalization, as stated in the
paper.

## License

Apache License 2.0; see [`LICENSE`](LICENSE).

## Citation

If you use these files, please cite the paper:

```bibtex
@misc{tomar2026uniqueness,
  author    = {Tomar, Tejasvi Singh},
  title     = {Uniqueness of convex central configurations in the planar four-body problem},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23004307}
}
```

To refer to the files themselves, cite this release:

```bibtex
@software{lean4body,
  author    = {Tomar, Tejasvi Singh},
  title     = {lean4body: programs, certificates and {Lean} formalization for the uniqueness of
               convex four-body central configurations},
  year      = {2026},
  version   = {v1.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23004297},
  url       = {https://github.com/tejasvi/lean4body}
}
```

[`CITATION.cff`](CITATION.cff) gives the same information in the Citation File Format.
