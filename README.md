# Integer multiplication: a conditional witness above 2^-13

**Conditional research draft by Swapnil Jain.**

This draft builds on OpenAI's *Integer multiplication below n log n* and on Douglas
Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). In
that model (a fixed finite alphabet and a fixed number of one-dimensional tapes),
it gives the conditional witness

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\boxed{\kappa={KAPPA_TEX}\approx{KAPPA_SCI}>2^{-13}}.
$$

That is about {FOLD_R7} times our round-seven witness `6.55177e-5`, {FOLD_R6} times our round-six witness `3.66657e-5`,
{FOLD_R5} times our round-five witness `1.54788e-5`, and about {FOLD_R1} times our first witness `7499/10^12`. These
figures compare asymptotic exponents, not practical runtimes.

## Round eight

`scripts/certificate_round8.py` certifies both savings from frozen child-width histograms
(`certificates/round8/`) and assembles them with `scripts/certificate_round3.py` (`beta = 1/1000`, crude guard).

**Bit side: opposite bank orders on our round-seven witness 2.** The schedule, frames and payload word of round
seven's headline witness (`h = 23`, `R = 27794`; Avi Eisenberg's producer, ikeboy, PR #62, compiled by eumemic's
joint frame compiler, PR #57) are unchanged; only the address geometry of each interchange is
recompiled with opposite bank orders (IceKylin, icekylinx, PR #104):

- In one common generic rational basis, every residual idempotent of rank `r` has nonzero leading principal minors
  `1..r` (a Zariski-density argument; all our residuals are rational idempotents because the frames are nested and
  `G`-nondegenerate). Storing the tail bank in reverse atom order turns the generic anti-diagonal Bruhat cell of each
  partial swap into **one contiguous reversed run**, so a residual of rank `r` is a single child of width `r`. The
  family closes under this reversal, and an outer wrapper of three lower updates restores the ordinary interchange.
- So the aux edge and its corner merge into one child `m - r_u`, every chain, `y_T`, `X_S` and centre step is one
  child, and the data entrance is one run of `m - 2h + 1`. `W` and `s` are exactly round seven's.
- The cross-bank adapters read a reversed atom index and are paid as atom loops, a linear toll. Stopping the
  recursion at atom width `e^theta`, `theta = 1/1000`, makes it subordinate, and the usable saving is
  `a_b = (1 - theta) a* + theta a_old` with `a_old = 12899/(2*10^8)`, our certified round-seven saving of the same
  witness (valid because `theta > a_b`).

This certifies the one-run saving `a* = {A_STAR}` and the stopped `a_b = {A_B}`.

**Complex side: completed-core sharing on our two-stage complex word.** Our `h = 24` complex word (retained
centres, the normal form for alternating residuals (IceKylin, icekylinx, PR #24), carrier matching (IceKylin,
icekylinx, PR #104), lifted frames, deferred
readouts, V leaves and the phase-1 closure) is built on the frozen `h = 24` addition DAG of eumemic (PR #117), read as
data, with every support recomputed. On top of it we reimplemented completed-core sharing (Andrey Mas, an664,
PR #128): the triples are split into 87
binary-orthonormal groups (83 of size 24, 4 of size 8; Andrey Mas's partition, after Zhang and Ge), the cores of a group run consecutively
on one bank of auxiliaries with dirty scratch, and the per-core wrap child is replaced by one fix-up child per
group and stream, merged into the core's last step where that is legal. This certifies `a_c = {A_C}`.

Together: `kappa = {KAPPA_FRAC}` ({KAPPA_SCI}), bound by the {BINDING} side.

Checks (`make round8`, standard library, minutes; `make round8-heavy`, a many-core machine with numpy and a C
compiler, hours):

- `scripts/certificate_round8.py`: both moment certificates with rigorous rational bounds
  (`ln(m/w)` from its series, `exp x <= 1 + x + x^2/2 + x^3/(6(1 - x/4))`), the grid point checked to be the largest,
  the stopped mix, and the assembly; `tests/test_round8.py`;
- `independent/round8-oppbank/`: `bit_hist.py` rebuilds the bit histogram from the frozen round-seven schedule and
  first reproduces round seven's certified saving; `run_gate.sh` checks, on **every** charged residual of the
  witness, that one 60-bit random integer basis gives nonzero leading minors (fast primes, every zero re-checked
  mod `2^61 - 1`), samples the reversed-run Bruhat profile, and runs negative controls; `wrap.py` checks the outer
  wrapper and remainder pairing exactly;
- `independent/round8-coreshare/`: `complex_hist.py` rebuilds the complex histogram from the word,
  `e2e_share.py` replays the shared word exactly over `Q(i)` with dirty scratch and controls, and
  `partition_check.py` checks the partition (exact cover, Gram `I`, even overlaps);
- `independent/round8-complex-gate/`: a walk of the actual op sequence and an exact `Q(i)` dirty-scratch replay of
  both stages of the unshared word, written without importing the word's builder.

The moment certificates, the stopped mix and the assembly are also checked in Lean 4's kernel (`lean/Round8.lean`,
`make lean`). Stage two as the complement time-reversal of stage one is assumed, as in rounds six and seven. The
opposite-bank compile needs only the existence of the common basis; the all-edges check exhibits one.

## Round seven (previous witness, `3275885357429/(5*10^16)`)

`scripts/certificate_round7.py` rebuilds the bit side's child-width histogram from frozen schedules
(`certificates/round7/`, about 8.2 MB compressed in all) and certifies it; the complex side is round six's. The bit
side at `h = 23` (`notes/deferred-readout.tex`) keeps the round-six interchange and adds:

- **Lifted frames and late copies.** Every addition acts at `M_n = span(n) + (U_n cap F_j)`, after Paureel's
  complement-frame construction, and a node with several free uses copies itself late, at the intersection of the
  frames of its remaining users, so the copy starts higher up its chain.
- **Deferred readouts (B-defer).** The second pass of stage one runs in two phases: first everything the retained
  totals need (their copies then read every target at `0`), then every slot that phase did not touch reads its
  garbage out at a nonzero frame `sigma_u` inside its start frame and every target's `t_T^perp`. Only the order of
  the word changes; it is still a mirror. Readouts go in increasing `dim sigma_u`, which keeps every target's frame
  sequence nested.
- **V leaves.** Each use of a leaf gets its own `V` gate on the data wire `X_S`, as a late copy there: the use slot
  starts at the intersection `s_i` of the frames of the remaining uses, and the deferred `V` gates follow in
  increasing dimension.

Two witnesses, with two side programs, both built on PR #62's producer (interval strips and core-aware pair
assembly) with PR #41's alternating order and links:

- **Witness 2 (the headline).** PR #69's balanced coarse sums in the side DAG, compiled by PR #57's joint frame
  compiler in PR #60's rank-first reclamation order (`R = 27794` roles). Here the first phase must also respect the
  order of frames along every role: it is the closure of read/write edges **and frame-order edges** (an op touching
  a role at a lower chain position runs first). It certifies `a_b = 32761/(5*10^8)` and
  `kappa = 3275885357429/(5*10^16)` with the staircase data-entrance corner (runs `h-2, h-5, 1^6`), and
  `a_b = 12899/(2*10^8)`, `kappa = 3224542033151/(5*10^16)` with the round-six entrance corner.
- **Witness 1 (an independent second construction).** PR #62's producer alone, with a moment-weighted choice of
  links after PR #44's weighted matching (`R = 28866` slots). It certifies `a_b = 31987/(5*10^8)`,
  `kappa = 1599247689723/(2.5*10^16)` (staircase corner) and `a_b = 12599/(2*10^8)`,
  `kappa = 6299103187973/10^17` (round-six corner).

In both, the rank sum is exactly round six's budget `s = Wm - N + L`. The checks are standard-library Python, one
process each (`make round7`):

- witness 2, `independent/joint-frame-stack/`: `check_frames.py` rebuilds every frame of the lifted program exactly
  over Q from the definitions and compares it with the frozen bases; `check_schedule.py` replays the compiled, lifted
  and final programs and the reordered stage-1 word symbolically over F2 with dirty scratch, recomputes the phase-1
  closure, and walks the actual reordered op sequence through every role's frame chain; `check_design.py` checks the
  deferred and V-leaf frames exactly, the late copies, the side lemma on every high-rank step the certificate
  charges, and rebuilds the histogram; all with negative controls;
- witness 1, `independent/deferred-readout/`: `check_word.py` (replays over F2 and Z), `check_frames.py` (exact
  frames, nesting in time order, side lemma on the new steps, histogram) and `check_lifted.py` (the base program);
- `independent/deferred-readout/check_stair.py` certifies the staircase entrance corner at the integer point used by
  the side-lemma checks.

Stage two is the complement time-reversal of stage one, as in round six; this is assumed, not machine-checked. The
moment certificates and all four assemblies are also checked in Lean 4's kernel (`lean/Round7.lean`, `make lean`).

## Round six (previous witness, `3666565558019/10^17`)

`scripts/certificate_round6.py` (bit side) and `independent/complex-twostage/` (complex side) assemble the witness:

- **Copied centres** (`notes/copied-centres.tex`), after PR #36's copied retained-centre schedule. A centre must be
  read at `D0` by the second scatter and then gathered at `D1`; instead, a temporary copy goes down to `D0` for the
  reads and is erased, and the original never leaves `D1`. Each centre loses one rank-`h` return per stage, the rank
  budget becomes `s = Wm - N + L`, and the deficit `v(v - 2h^2)` is positive at much smaller `h`. On our bit
  interchange (flag basis, gm side circuit, side roles batched one level down, data-entrance run) this certifies
  `a_b = 34919/10^9` at `h = 25`.
- **Retained point totals** replace the direct centre wires: each total is a side role built from existing side
  nodes, its copy pays rank `h-1` and the original rank `1`, so the loss per centre drops from `h` to `h-1` and the
  centre roles leave `W` (`independent/two-stage-bit/rtgm.py`). This certifies `a_b = 36667/10^9` at `h = 23`.
- **The two-stage complex interchange** with our pair-exclusion producer, copied retained centres, complex data-edge
  batching and every residual batched as one whole-residual child, certifies `a_c = 36926111/(5*10^11)` at `h = 24`
  (`cd independent/complex-twostage && python3 run.py 24`: rank sum exact, labels 0 bad, `a_c >= 73861113/10^12`).

The bit side binds: `kappa = 3666565558019/10^17`. Both moment certificates and the assembly are also checked in
Lean 4's kernel (`lean/Round6.lean`, `make lean`), as in round five.

## Round five (previous witness, `309575208081/(2*10^16)`)

`scripts/certificate_round5.py` assembles the witness from our own two-stage bit interchange (Paureel's
two-stage motif), with no outside bit network:

- **A common flag basis** (`notes/flag-basis.tex`). One rational basis makes the corners of every auxiliary and
  centre edge of the two-stage motif lower triangular, in both stages at once, so each role's `h` corner
  transpositions become one contiguous run instead of `h` singleton calls. The same basis gives the data edges
  their batched profiles one level down, including one run of width `h-2` inside the data-entrance corner (a
  triangular corner plus a rank-one term). `independent/two-stage-bit/flagbasis.py` checks every edge class
  exactly, and `flag_existence.py` checks the nonvanishing conditions at the working dimensions `h = 46, 47`.
- **A cheaper side circuit** (`independent/two-stage-bit/gmside.py`): a global matching of the points shares
  pair sums and exclusion chains between centres, cutting the side roles at `h = 47` from 418678 to 403248;
  `checkgm.py` re-derives every support and output exactly.
- PR #7's complex network, fully batched with complex source frames, as in round four.

- **Side roles batched one level down.** In the flag basis every side-role edge has the form `pi (x) P_b`, whose
  profile is that of the `h x h` idempotent `pi`: an edge of rank `r > h/2` becomes `h-r` singletons and one
  block of width `2r-h`. `independent/two-stage-bit/side_chains.py` derives every side role's frame chain from
  the compiled side program.

The two moment certificates are also checked in Lean 4's kernel (`lean/Round5.lean`, `make lean`, core Lean
only): both child-width histograms sum to their total rank, the rational upper bounds on `ln(m/w)` are computed from
their series, the moment bound holds at the claimed savings, and the assembly (every constraint row and cost margin of
`scripts/certificate_round3.py`, with the guard constant bounded through the same ln series) gives the headline
`kappa`. The analytic facts behind those bounds (the atanh
tail and `(m/w)^a <= 1/(1 - a ln(m/w))`) are stated premises, not formalised.

Together these certify `a_b = 15479/10^9` at `h = 47` and `kappa = 309575208081/(2*10^16) > 2^-16`.

## Round four (previous witness, `59861145819/(5*10^15)`)

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
| Round-seven bit side (deferred readouts, V leaves, lifted frames, late copies; two side programs) | Frozen schedules and frames (`certificates/round7/`); stdlib replays over F2 (and Z for witness 1) with arbitrary scratch, walks of the reordered op sequence, exact frame checks over Q, side lemma at one integer point, histogram rebuilt from the schedule, negative controls (`independent/deferred-readout/`); Lean kernel check of the moment certificates and assembly. Stage two as the complement time-reversal of stage one is assumed |
| Round-eight bit side (opposite bank orders, stopped interchange) | Same frozen schedule as round seven; histogram rebuilt from it; all-edges common-basis check with fast-prime passes re-checked mod `2^61 - 1`, Bruhat sample and negative controls (`independent/round8-oppbank/`); Lean kernel check of the moment, stopped mix and assembly. The existence of the common basis, the factorization lemma and the stopped recurrence follow IceKylin's written argument (icekylinx, PR #104) |
| Round-eight complex side (completed-core sharing) | Histogram rebuilt from the word; exact `Q(i)` replays of the word and of the sharing with dirty scratch and controls; partition checked (`independent/round8-coreshare/`, `independent/round8-complex-gate/`) |
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
make round7      # round-seven certificate and checks of both witnesses (~45 min)
make round8      # round-eight certificate, bit histogram rebuild, partition and wrapper checks (minutes)
make round8-heavy  # round-eight word replays and the all-edges bit gate (many cores, numpy, a C compiler)
make lean        # Lean 4 kernel checks of rounds five to eight
make notes       # PDF notes (tectonic)
```

The PDFs are in `artifacts/`.

## Attribution

Author: **Swapnil Jain**. Research, implementation and drafting were done with
assistance from Claude (Anthropic). This builds on OpenAI's manuscript and Douglas
Colkitt's framework. It also uses published results from open pull requests on
Colkitt's repository, cited in `NOTICE`; those are citations of prior work, not
collaborators. Licensed under Apache-2.0.
