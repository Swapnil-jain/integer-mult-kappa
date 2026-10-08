# Past the a_b/2 ceiling for integer multiplication

**Conditional research draft by Swapnil Jain.**

This draft builds on OpenAI's *Integer multiplication below n log n* and on Douglas
Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). In
that model (a fixed finite alphabet and a fixed number of one-dimensional tapes),
it gives the conditional witness

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\boxed{\kappa=\frac{7499}{10^{12}}=7.499\times10^{-9}>2^{-27}}.
$$

That is about 2.01 times the strongest open pull request on that repository
(#7, `373/10^11`), and about 1.51 times our previous `249/(5*10^10)`. These figures
compare asymptotic exponents, not practical runtimes.

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
| B. Segmented inverse | Inside any wrap-free block, the Gaussian correction `N` is exactly `D_c T_c D_c^{-1}` for any centre `c`, and every `T_c` is a geometric rescaling of one Toeplitz matrix. Recentred sub-blocks keep the chirp precision bounded. Gohberg–Semencul inverts each block, and Woodbury couples the cuts and wraps. With oversampling `1/(4 d L^eps)`, this replaces about `30d` Neumann terms with `O(1)` work per output. | `notes/segmented-inverse.tex` |
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
| Full upstream multiplication theorem | Assumed |
| Independent review / formalisation | Not supplied |

## Reproduce

Requires Python 3.9 or newer, standard library only.

```sh
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
