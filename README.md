# Past the a_b/2 ceiling for integer multiplication

**Conditional research draft by Swapnil Jain.**

This draft builds on OpenAI's *Integer multiplication below n log n* and on Douglas
Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). In
that model (a fixed finite alphabet and a fixed number of one-dimensional tapes),
it gives the conditional witness

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\boxed{\kappa=\frac{249}{5\cdot10^{10}}=4.98\times10^{-9}>2^{-28}}.
$$

That is about 1.34 times the strongest open pull request on that repository
(#7, `373/10^11`). These figures compare asymptotic exponents, not practical
runtimes.

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
| B'. Segmented inverse | Between wraps, the Gaussian correction `N` is exactly diagonal × Toeplitz × diagonal-inverse. Gohberg–Semencul handles each segment, and Woodbury handles the cyclic coupling between them. This replaces about `30d` Neumann terms with `O(sqrt d)` work per output. | `notes/segmented-inverse.tex` |
| C. One chunk per axis | `K = l - 1` is possible once compact controls remove the `K^tau` factor. The transform layout then needs no exchanges. | `notes/stack-notes.tex` §2 |
| A. Fine-bit exposure | Only `O(log L)` low bits of an axis move innermost for resampling. The line passes stream over slabs. | `notes/stack-notes.tex` §3 |
| E'. Permuted selected-bit swaps | Colkitt's selected-bit addition with a permuted pairing reverses the CRT axis order at cost `l (d log p)^tau`. | `notes/stack-notes.tex` §4 |

Only the combination moves the headline. D alone, for example, leaves PR #5's
Gaussian row capping `eps < 1/2`. The binding constraint is now the butterfly row,
`eps(1-lambda')`, at `eps = 0.666`.

## Evidence and scope

| Component | Evidence |
| --- | --- |
| Parameters, margins, side constraints | Exact rational certificate (`scripts/certificate.py`) with negative controls |
| Closed form `X = n(n+2 beta_j)`, segment factorisation, cross-wrap bound `X >= 2` | Exact rational checks on five `(s,t)` pairs (`scripts/check_identities.py`) |
| Segmented inverse | Written proof (`notes/segmented-inverse.tex`); a 260-digit comparison against a direct solve, with cyclic indices (`scripts/check_segmented_inverse.py`) |
| A, C, D, E' | Written arguments, each attacked by a refute-by-default check, with the resulting fixes in the notes |
| Prior results from PR #3, #5 and #7 on integer-mult-bounds (networks, chirped Gaussian lemma) | Assumed, unmerged there |
| Full upstream multiplication theorem | Assumed |
| Independent review / formalisation | Not supplied |

## Reproduce

Requires Python 3.9 or newer, standard library only.

```sh
make verify
```

## Attribution

Author: **Swapnil Jain**. Research, implementation and drafting were done with
assistance from Claude (Anthropic). This builds on OpenAI's manuscript and Douglas
Colkitt's framework. It also uses published results from open pull requests on
Colkitt's repository, cited in `NOTICE`; those are citations of prior work, not
collaborators. Licensed under Apache-2.0.
