# Integer multiplication: a conditional witness above 2^-17

**Conditional research draft by Swapnil Jain.**

This draft builds on OpenAI's *Integer multiplication below n log n* and on Douglas
Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). In
that model (a fixed finite alphabet and a fixed number of one-dimensional tapes),
it gives the conditional witness

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\boxed{\kappa=\frac{59861145819}{5\cdot10^{15}}\approx1.19722\times10^{-5}>2^{-17}}.
$$

That is about 2.86 times our round-three witness `4.1878e-6`, and about 1600 times our first
witness `7499/10^12`. These figures
compare asymptotic exponents, not practical runtimes.

## Round four

`scripts/certificate_round4.py` assembles the witness from:

- the stack below, which drives `eps` toward 1 (the bit side binds);
- PR #24's endpoint-gauge bit network, built on PR #18's partial-swap networks. We take its published
  counts and child-width multiplicities (`certificates/external/pr24-bit-network.json`), check the rank sum
  and deficit, and certify `a_b = 4788949/(4*10^11)` with our own moment bisection;
- PR #7's complex network at `h=28`, fully batched, with **complex source frames**
  (`notes/complex-source-frames.tex`): starting each stage-two auxiliary role in the phase frame of `D0`
  merges its entrance and exit into one child of rank `m-h`. Our histogram certifies
  `a_c = 4079603/(2.5*10^11)`, more than the bit side needs. PR #21 contains the same translation,
  independently and in a more complete form.

The same script also reports a witness using only our own bit side: the two-stage interchange with
**data-edge batching** (`notes/data-edge-batching.tex`, `independent/two-stage-bit`). Both stage-two data
entrances have idempotent difference `Q1 (x) Q2`, of rank `(h-1)^2 > m/2`, so the partial-swap batching
lemma applies to them too; at `h=47` this certifies `a_b = 10033/10^9` and `kappa = 1003289933971/10^17`.

## Round three (previous witness, `13086957581/(3.125*10^15)`)

`notes/round3-combination.tex` assembles the witness from four components:

- the stack below, which drives `eps` toward 1;
- a batched two-stage bit interchange, using the partial-swap batching lemma in
  `notes/partial-swap-batching.tex` (`a_b = 22157/(5*10^9)`);
- PR #7's complex network at `h=28`, fully batched as in PR #15 (every residual edge is one
  whole-residual child), rebuilt and certified independently in `independent/complex-network`
  (`a_c = 1048009/(2.5*10^11)`; this side binds);
- a coefficient-depth guard that needs no path-topology assumption.

`scripts/certificate_round3.py` checks every row and constraint exactly.

## The ceiling and how this gets past it

After the recent pull requests, every witness is capped at about half the
bit-network saving, `kappa < a_b/2`. Three costs each tie a saving to either the
number of axes `d` or the axis width `l`, and `d*l ~ log n`:

- the butterflies save on `d`;
- the chunk layout and the axis layout save on `l`.

On top of that, the Gaussian precision budget `~8d^2` must fit in `O(log n)`.

Five changes remove those caps.

| Change | What it does | Note |
| --- | --- | --- |
| D. Longer digits | Digits of `(log n)^(1+x)` bits separate the precision budget from the address length. The volume stays `Theta(n)`. | `notes/stack-notes.tex` §1 |
| B. Segmented inverse | Inside any wrap-free block, the Gaussian correction `N` is exactly `D_c T_c D_c^{-1}` for any centre `c`, and every `T_c` is a geometric rescaling of one Toeplitz matrix. Recentred sub-blocks keep the chirp precision bounded. Gohberg–Semencul inverts each block, and Woodbury couples the cuts and wraps. With oversampling `1/(4 d L^eps)`, this replaces PR #5's `O(d)` Neumann terms per output (its counts are at most `30d`; the burn-in alone is about `1/theta = 4d`) with `O(1)` work per output. | `notes/segmented-inverse.tex` |
| C. One chunk per axis | `K = l - 1` is possible once compact controls remove the `K^tau` factor. The transform layout then needs no exchanges. | `notes/stack-notes.tex` §2 |
| A. Fine-bit exposure | Only `O(log L)` low bits of an axis move innermost for resampling. The line passes stream over slabs. | `notes/stack-notes.tex` §3 |
| E'. Permuted selected-bit swaps | Colkitt's selected-bit addition with a permuted pairing reverses the CRT axis order at cost `l (d log p)^tau`. | `notes/stack-notes.tex` §4 |

Only the combination moves the headline. D alone, for example, leaves PR #5's
Gaussian row capping `eps < 1/2`. The binding constraint is now the butterfly row,
`eps(1-lambda')`, at `eps = 0.9999`.

## Evidence and scope

| Component | Evidence |
| --- | --- |
| Parameters, margins, side constraints | Exact rational certificate (`scripts/certificate.py`) with negative controls |
| Closed form `X = n(n+2 beta_j)`, segment and recentred factorisations, cross-wrap bound `X >= 2` | Exact rational checks on five `(s,t)` pairs (`scripts/check_identities.py`) |
| Segmented inverse | Written proof (`notes/segmented-inverse.tex`); 260-digit comparisons against a direct solve, with cyclic indices, for whole segments (`scripts/check_segmented_inverse.py`), recentred sub-blocks (`scripts/check_subblock_inverse.py`), and the reuse path, which uses one stored inverse rescaled per block, elimination without pivoting, and a fixed-point negative control that must fail (`scripts/check_reuse_inverse.py`) |
| A, C, D, E' | Written arguments, each attacked by a refute-by-default check, with the resulting fixes in the notes |
| Prior results from PR #3, #5 and #7 on integer-mult-bounds (networks, chirped Gaussian lemma) | Assumed, and unmerged there. We recomputed PR #7's W, m, N, L, s, eta and their complex counterparts exactly from its formulas. An independent re-implementation (`independent/pr7-role-counts`, `make roles`) reproduces its h=28 role counts, 11840940 and 93838, exactly. The finite-alphabet version of the interchange lemmas that its ternary payloads need is proved in `notes/stack-notes.tex`, Appendix A. |
| Partial-swap batching (bit side) | Written proof (`notes/partial-swap-batching.tex`); exact Bruhat profiles on random, sparse and adversarial idempotents, and two-stage frames at h=6,7,8 (`independent/partial-swap`) |
| Fully batched complex saving | Independent rebuild of PR #7's producer, exact label checks (0 bad edges at h=28), the residual rank of every edge with an exact sum check against s, and a certified moment bisection (`independent/complex-network`) |
| PR #24's bit network | Its published counts and child-width multiplicities, pinned by commit and SHA-256; we check that the widths sum to the total rank, the deficit, and that every child is narrower than m, and certify the moment ourselves. We also re-derived its frame identities and confirmed it is compatible with our stack (separate bit and complex arities, crude guard, payload and prime choice); its producers and common basis are assumed |
| Data-edge batching (bit side) | Written proof (`notes/data-edge-batching.tex`): closed form `Q1 (x) Q2` checked in exact rationals; corners invertible under one common basis for every auxiliary edge and every data pair at h=6,7,8 (all 3136 pairs at h=8, `corners.py`); exact Bruhat profiles (`bruhat8.py`) |
| Two-stage side circuit at h=47 | Our generator (`independent/two-stage-bit/sidegen.py`) reproduces the published h=32 count, 123157, and an exact checker (`checkside.py`) verifies supports, disjoint children, common points and every output at h=28 to 50; `moment.py` builds the child histogram role by role and checks it sums to s |
| Complex source frames | Written proof (`notes/complex-source-frames.tex`); exact label and phase checks of every endpoint case (`independent/complex-network/sourceframe_labels.py`); a three-stage scalar simulation of PR #7's network with source frames, arbitrary scratch and a negative control (`sourceframe_sim.py`); the histogram option `sf` keeps the exact rank sum s |
| Topology-free guard | Written argument (docstring of `scripts/certificate_round3.py`); the witness is stated with this guard |
| Prior results from PR #10, #13 and #15, and the two-stage motif | Assumed. We reproduced their rank moments, certified at least PR #15's fully batched saving from our own histogram, and found controlled-basis witnesses for PR #10 at h=8 to 14 |
| Full upstream multiplication theorem | Assumed |
| Independent review / formalisation | Not supplied |

## Reproduce

Requires Python 3.9 or newer, standard library only.

```sh
python3 independent/complex-network/fullbatch_hist.py 28 /tmp/fb28sf.json sf   # ~4 min
python3 scripts/certificate_round4.py /tmp/fb28sf.json
python3 independent/two-stage-bit/checkside.py 47
make verify      # exact checks, numerical inverse checks, certificate, tests
make roles       # independent recount of PR #7's role counts (clang++, ~1 GB)
make notes       # PDF notes (tectonic)
```

The PDFs are in `artifacts/`.

## Attribution

Author: **Swapnil Jain**. Research, implementation and drafting were done with
assistance from Claude (Anthropic). This builds on OpenAI's manuscript and Douglas
Colkitt's framework. It also uses published results from open pull requests on
Colkitt's repository, cited in `NOTICE`; those are citations of prior work, not
collaborators. Licensed under Apache-2.0.
