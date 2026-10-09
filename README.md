# Integer multiplication: a conditional witness above 2^-11

**Conditional research draft by Swapnil Jain.**

This draft builds on OpenAI's *Integer multiplication below n log n* and on Douglas
Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). In
that model (a fixed finite alphabet and a fixed number of one-dimensional tapes),
it gives the conditional witness

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\boxed{\kappa=\frac{6685733}{10^{10}}\approx2^{-10.547}\approx6.68573\times10^{-4}>2^{-11}}.
$$

That is about 1.09 times our round-ten witness `6.15338e-4` (`2^-10.666`), 1.43 times our round-nine witness `4.66374e-4`,
5.30 times our round-eight witness `1.26130e-4`, 10.20 times our round-seven witness `6.55177e-5`, and about 89200 times our first
witness `7499/10^12`. These figures compare asymptotic exponents, not practical runtimes.

## Round eleven

`scripts/certificate_round11.py` runs round ten's certificate unchanged (the paired-cube /
three-stage cover assembly of icekylinx's PR #144, on measured per-side inventories). Before certifying ours it
reproduces PR #144's published `kappa = 4609169/10^10`, our round-nine `kappa = 4663738/10^10` and our round-ten
`kappa = 6153378/10^10`. Both sides are new words, and both are measured histograms from independent replays.

**Complex side: paired cubes with dead-copy recycling.** The base is icekylinx's paired-cube decoder (PR #144),
rebuilt in our code at `p = 11` (`h = 22`, `m = 66`), as in round ten. The changes:
- the triple, pair-disjoint and all-but-one modules are eumemic's PR #168 v4 modules, used as data and unchanged
  (pinned by SHA-256 in `NOTICE`). Their contracts are re-checked when the word is built;
- the cube-local circuit is configurable: the 13 local outputs are fixed by PR #144's decoder identity, and only how
  each is summed is chosen (the configurable circuit and its coordinate space are icekylinx's, in the PR #168
  lineage; chafreaky's PR #181 chose a configuration in the same class). Our choice, C1, moves the third long-diagonal
  channel to `G[0,2,1]`. It was found by our own search and is about 1.3% better than PR #168's L1 in the same stack;
- carrier links are extended under DaysSky's relaxed acceptance condition (PR #162), reimplemented, with our greedy
  scan run in reverse order;
- birth reuse (jamesyc's PR #124) with late compensation reads (eumemic's PR #143), by our exact maximum-weight
  matching. All 2,332 gauges are paired, so NDS adds nothing;
- a per-operation physical descent with the pairs fixed (ours; the bundle-descent idea is eumemic's, PR #168);
- 42 terminal-output deletions: jamesyc's lemma (PR #166), as used by chafreaky (PR #176), with our selection;
- dead-copy recycling. The word holds pairs of registers that retire holding the same node value. One gate
  `P -= Q` at a frame containing both copies' last frames leaves `P` with zero signal and only dirty contents. A fresh
  birth is then placed on `P` with no read, and the frame-0 compensation of every non-gauged role is recomputed over
  the new transcript. This is Section 2a of icekylinx's PR #184 (zero-fresh recycling through a containing frame),
  used here with no source controls. The kernel comes from the duplicated copies. Ours are the application to this
  word, the host selection (an exact maximum-weight matching between copy pairs and births), the ledger and the
  checker. 386 births are hosted. Each removes one unit of `W` and leaves the deficit unchanged.

Completed-core sharing (an664, PR #128) is unchanged from PR #144. The result has `W = 12,922` (it was
13,308 before recycling), deficit `1,320` and largest child `20`. Certified saving: `a_c = 67147467/10^11`.

**Bit side: the PR #168 v4 bit word with our module and levers.** The base is eumemic's paired-cube bit decoder
(PRs #155 and #168: transversal-triple ports, broadcast roots, partner-pair mixing), after icekylinx's PR #144, at
`p = 12` (`h = 24`, `m = 72`). With v4's own modules, our builder reproduces PR #168 v4's bit graph byte
for byte. On it:
- the pair-disjoint module is PR #168 v4's annealed bit pair module (eumemic), used as data and unchanged;
- the local-channel association is eumemic's L1 (PR #168), reimplemented as constants;
- the all-but-one module is our cyclic-interval module (plan `1,2,3,5,7,9 / 01000`), in place of v4's annealed one;
- merged face and edge reads, as in eumemic's PR #161; carrier links by our cover-weighted matching;
- every gauge candidate kept, then the omitted gauges chosen by exact minimum cut (Th0rgal's formulation, PR #146),
  with NDS on the kept gauges;
- node, group and per-operation frame descent (ours);
- birth reuse with late births (the lineage of jamesyc's PR #124), by exact maximum-weight matching: 1,762 pairs;
- 24 terminal sinks: a sink writes onto a pivot target instead of keeping a register (jamesyc's idea, PR #166;
  our lemma and selection, and every sink is re-checked by the ledger).

This gives `W = 380341/20 ≈ 19,017.05`, deficit `1,936`, the coarse saving `a* = 670100165/10^12`, and, stopped at `theta = 1/10^3`
with `a_old = 384599/10^10`, `a_bit = 669468524735/10^15`.

Together: `kappa = 6685733/10^10` (about 6.68573e-04, `2^-10.547`), bound by the bit side. The minimum
margin `eps q = 6.6857332707e-04` is rounded down to the `10^-10` grid. That is about 1.09 times our round-ten
witness.

**Proof interfaces.** Round ten's interfaces still apply: the three-stage cover lifting, the exterior rule and the
stopped recurrence (PRs #130, #144); the `2^-P 3^-K` precision grid; block-factored frame representatives; stage two
as the complement time-reversal of stage one; NDS with cube-type representatives; and a chain spliced across a birth
hand-off charged as one ascending chain (PRs #124, #143). New this round, stated here but not machine-checked:

- dead-copy recycling rests on the zero-fresh recycling argument of PR #184, Section 2a, with no source controls.
  It is checked by F2 nesting of every chain (both copies' climbs to the mix frame included), the recycling
  chronology, an exact ledger recount, and a dirty-scratch replay mod `2^61 - 1` with the compensation recomputed.
  There is no literal frame replay of the mix gate at any size;
- terminal-output deletions (complex) and terminal sinks (bit) rest on PR #166's lemma. The complex checker tests
  sink eligibility and replays; the bit ledger checks every absorbed write exactly;
- on the complex side, the literal frame replay of the stack is round ten's (nested prefix, merged reads, birth reuse
  at `p = 6` and `7`). The C1 circuit, the per-operation descent and the recycling are checked at real size by the F2
  recount and the scalar replay only;
- the paired-cube bit word's partner-pair mixing lemma and its frame theorem follow eumemic's written argument (PRs
  #155, #168), checked here through the replays below. NDS on the bit side is the per-slot construction; no literal
  NDS replay on bit frames has been done;
- the module, circuit, gauge and host searches were scored with cheaper objectives. Only the final words are
  certified under the full objective;
- the external row reserve of the finite bridge keeps PR #144's inherited bit constants (`ROW_RESERVE` = 9909 + 252: an earlier round's coarse bit family and the ordinary leaf). Neither is a term of the current bit cover, whose own row stock is borrowed and restored internally at any fixed `m` (PR #144's rows argument, used here at `m = 72`), so they carry over.

Checks (`make round11`, standard library, seconds; `make round11-heavy`, a machine with numpy, about 3 GB):

- `scripts/certificate_round11.py` runs the regressions and certifies ours, and `tests/test_round11.py` runs round
  ten's frozen-assembly tests on the round-eleven files;
- `independent/round11-bit-gate/`: `ledger.py` replays the frozen bit word from JSON alone and imports no
  construction code. It checks an exact F2 replay with dirty scratch, forward and reflected; exact frames over Q and
  their nondegeneracy; every absorbed-sink write; and the event histogram against the accounting. It runs on the
  `p = 12` word and on the same recipe at `p = 9` and `7`, and eleven mutation controls must be rejected on the
  `p = 7` word, which has sinks. `cert.py`
  is the rational certificate with the next grid point rejected, and `invcheck.py` checks that the release inventory
  equals the replay's measured histogram;
- `independent/round11-cx-recycle-gate/`: `check_word2_recycle.py` reads only the frozen complex word. It checks the
  decoder identity, sink eligibility, the recycling chronology, F2 legality of every role and target chain, the
  ledger recounted from the chains, the deficit `2v - 3 loss`, the `10^-12` certificate with the next point rejected,
  and a dirty replay mod `2^61 - 1` on two seeds. Its controls (birth, terminal, mix skipped, mix sign flipped) must
  fail, and `mutate_hosts.py` writes three mutated words (unequal copies, an early birth, a tampered claim) that must
  be rejected. It runs at `p = 11`, `9` and `7`. `README.md` there describes the mechanism and its ledger;
- `independent/round11-cx-recycle-build/`: the construction, with the pinned modules and their notice in `data/`.

The moment certificates, the stopped mix, the finite bridge and the 47-row assembly are also checked in Lean 4's
kernel (`lean/Round11.lean`, `make lean`).

## Round ten (previous witness, `6153378/10^10`)

`scripts/certificate_round10.py` runs round nine's certificate (the paired-cube / three-stage cover assembly of
icekylinx's PR #144) on cover profiles with fractional child weights and on measured per-side inventories. Each
side's profile is built in layers, and every layer is checked: PR #144's parts, or a histogram measured by a replay;
then nullity-dependent sharing (NDS), rebuilt from each gauge's type; then birth reuse, measured from the actual role
chains. Every layer must keep the telescoping deficit `D = W m - s = 2v - 3 loss`. Before certifying ours, the script
reproduces PR #144's published `kappa = 4609169/10^10` and our round-nine `kappa = 4663738/10^10`, both from the
parts and from the bare histograms.

**Complex side: paired cubes with nested-prefix modules, birth reuse and per-operation descent.** This is icekylinx's paired-cube decoder (PR #144: the signed channel graph, coordinate frames, the chronological
partial-gauge filter and copied centres), rebuilt in our code from its data and notes, at `p = 11`
(`h = 22`, `m = 66`). The changes:
- the triple-exclusion module is our own searched restriction of eumemic's `h = 24` DAG (PR #117, used as data);
- the pair-disjoint module is eumemic's annealed module G37 from PR #168, used as data and unchanged (pinned in
  `NOTICE`). Its contract is re-checked when the word is built;
- the all-but-one module is eumemic's nested prefix (PR #168), reimplemented from its description. Our version is
  byte-identical to theirs;
- face and edge outputs share merged reads, as in eumemic's PR #161, reimplemented;
- carrier links come from a maximum matching under PR #144's test. They are then extended greedily under
  DaysSky's relaxed acceptance condition (PR #162), reimplemented: an arc is added only if the graph stays acyclic and
  the span stays inside the backward intersection;
- frames descend node by node (ours; the physical-frame idea is eumemic's PR #131);
- gauges are re-selected pair-aware (ours): a pairable candidate is charged its hand-off instead of its tail;
- birth reuse (jamesyc's PR #124) with late compensation reads (eumemic's PR #143), reimplemented with our exact
  maximum-weight matching. 3,180 of 3,180 gauges are paired;
- a final per-operation physical descent with the pairs fixed, alternated with re-matching (the bundle-descent
  idea of eumemic's PR #168, reimplemented with our exact recount as its objective).

The result has `W = 14,592`, deficit `1,320` and largest child `20`. NDS adds nothing here, because every
gauge is paired. Completed-core sharing (an664, PR #128) is unchanged from PR #144. Certified saving:
`a_c = 61609676/10^11`.

**Bit side: the paired-cube bit word with our levers.** The base is eumemic's paired-cube bit decoder (PR #155: transversal-triple ports, broadcast roots, partner-pair source mixing on the X registers) on icekylinx's face and edge modules (PR #144), at `p = 12` (`h = 24`, `m = 72`). Our builder re-implements it from its description and reproduces its published graph and child histogram exactly. The pair-disjoint module is eumemic's annealed module from PR #161, used as data (pinned in `NOTICE`). On it we stack our own levers, except where credited:
- merged reads: one shared read node per pair of face and edge outputs, read at the singleton caps;
- carrier links chosen by an exact maximum-weight matching on the cover objective;
- every gauge candidate kept, then the omitted gauges chosen by exact minimum cut (Th0rgal's formulation, PR #146), with NDS on the kept gauges;
- frame descent;
- birth reuse with late births, by exact maximum-weight matching (the lineage of jamesyc's PR #124 and DaysSky's PR #150). It needs the merged reads: with unmerged reads no port label lies in any gauge's `sigma`, so no pair is compatible. Here 1,767 pairs are reused.

This gives `W = 257897/12 ≈ 21,491.42`, deficit `1,936`, the coarse saving `a* = 632988409/10^12`, and, stopped at `theta = 1/10^3` with `a_old = 384599/10^10`, `a_bit = 632393880491/10^15`.

Together: `kappa = 6153378/10^10` (about 6.15338e-04), bound by the complex side. The minimum margin
`eps q = 6.1533781105e-04` is rounded down to the `10^-10` grid.

**Proof interfaces.** Round nine's interfaces still apply: the three-stage cover lifting, the exterior rule and the
stopped recurrence (PRs #130, #144); the `2^-P 3^-K` precision grid; block-factored frame representatives; the
frozen PR #117 DAG; and stage two as the complement time-reversal of stage one. New this round, stated here but not
machine-checked:

- NDS needs cube-type (block-factored) representatives. It fails with the general frame representatives of PRs
  #130, #131 and #137. Our literal NDS replay covers cube-type `sigma`s at `h = 8`;
- on the complex side, birth reuse, merged reads, relaxed arcs and the per-operation frames are checked at real size
  by an F2 recount from the role chains and by an exact scalar replay with dirty scratch, both with controls. Nested
  role chains imply that the transitions are realisable at the stated ranks. The literal frame replay covers the
  mechanism at `p = 6` and `7` only, not at `p = 11`;
- a chain spliced across a birth hand-off (the donor's last frame inside the recipient's gauge) is charged as one
  ascending chain. This is the interface of PRs #124 and #143;
- the cube-local K block and the source histogram `[1, 2, h - 4]` are PR #144's;
- the module and gauge searches were scored with cheaper objectives. Only the final words are certified under the
  full objective;
- the paired-cube bit word's partner-pair mixing lemma and its frame theorem, from eumemic's written argument (PR #155), checked here only through the replays below;
- NDS on the bit side is the per-slot construction; no literal NDS replay on bit frames has been done;
- the external row reserve of the finite bridge keeps PR #144's inherited bit constants (`ROW_RESERVE` = 9909 + 252: an earlier round's coarse bit family and the ordinary leaf). Neither is a term of the current bit cover, whose own row stock is borrowed and restored internally at any fixed `m` (PR #144's rows argument, used here at `m = 72`), so they carry over.

Checks (`make round10`, standard library, seconds; `make round10-heavy`, a machine with numpy and a few GB):

- `scripts/certificate_round10.py` runs the regressions and certifies ours, and `tests/test_round10.py` checks it.
  The tests include a truth table for the kept-copy rule below;
- `independent/round10-cx-prefix-gate/`: `check_word.py` reads only the frozen complex word. It checks the decoder identity mod 2^61 - 1, the F2 legality of every operation frame, spliced role chain, birth pair and read order, the child-histogram and `W` recount with the deficit `2v - 3 loss`, and the `10^-12` certificate with the next point rejected. It also runs an aliased dirty-scratch replay on two seeds with three controls. `mutate_word.py` writes four mutated words (sign, frame, claim, pair), and each must be rejected. `lit_prefix.py` is the literal small-`p` replay of the nested prefix, merged reads and birth reuse, with three controls;
- `independent/round10-cx-prefix-build/`: the construction. `final2.py` rebuilds the word, and `export_word.py`, the only adapter, writes the frozen JSON the gate reads. The pinned pair module and its notice are in `data/`;
- `independent/round10-bit-gate/`: `ledger.py` replays the frozen bit word from JSON alone and imports no construction code. It checks an exact F2 replay with arbitrary dirty scratch, forward and reflected; every frame re-canonicalised over Q, with strictly nested moves; data inside the frame at every gate, read, centre copy and birth; nondegeneracy of every frame; birth hand-offs only after the donor's last gate; and the event histogram against the accounting. It runs seven mutation controls. `cert.py` is the rational certificate and the stopped mix, and `invcheck.py` checks that the release inventory equals the replay's measured histogram;
- `independent/round10-standing-gate/`: a second gate written separately from both builds. `cxguard.py` recomputes the finite semantic guard from PR #144's formulas at the complex word's `m`, and `twosided.py` rejects the next certificate point with a rigorous lower bound.

**Checker fix.** The repository's F2 checkers (`independent/joint-frame-stack/check_schedule.py` and
`independent/round9-bit-ledger/ledger.py`) used to identify a kept V copy by its (copy, source) role pair. Words
compiled with skip-suffix strips and the gm order gate that same pair a second time, later and at a higher frame.
The old rule misread that later gate as the copy, which skipped its walk and data checks. Kept copies are now
identified by op index (`jfdata.kept_copy_ops`). On witness 2 both rules pick the same 4,210 ops, so rounds seven to
nine are unchanged, and their checks still pass.

The moment certificates (fractional weights scaled to integers, with the bit fallback), the stopped mix, the finite
bridge and the 47-row assembly are also checked in Lean 4's kernel (`lean/Round10.lean`, `make lean`).

Support compute: https://buymeacoffee.com/jainswapnil138

## Round nine (previous witness, `4663738/10^10`)

`scripts/certificate_round9.py` is our own implementation of the paired-cube / three-stage cover assembly of
icekylinx's PR #144 (merged upstream through PR #149). It reads frozen inventories (`certificates/round9/`), rebuilds
each side's cover child histogram from its parts, certifies both savings with rigorous rational bounds, applies the
stopped mix, checks the finite router bridge, and checks all 47 strict constraints and 7 cost margins. Before
certifying ours, it reproduces PR #144's published `kappa = 4609169/10^10` from PR #144's published inventories
(`certificates/round9/pr144/`, read as data), matching every constraint and margin exactly. None of PR #144's code is
run.

**The cover (icekylinx, PRs #130 and #144).** Each interchange is lifted to a three-stage Cayley cover with
`m = 3h` coordinates, and one physical role per slot runs through all three stages, so `W = 2v + R` roles per
vertex. PR #144 adds paired cubes on the complex side and shared completed cores (the sharing principle of an664's
PR #128). The sequential dirty reuse across stages follows eumemic's padded triple covers (PR #137), and the
complex supplier's lineage includes jamesyc's birth-read slot reuse (PR #124).

**Bit side: our witness 2 with gauge omission.** We put our round-seven headline witness (witness 2: deferred readouts,
V leaves, lifted frames and late copies at `h = 23`, `R = 27794`) into the `m = 69` bit cover in place of the PR #97
word. A retained deferred gauge costs an exterior child of width `3f`. An omitted one reads its old value at frame 0
in a prelude instead, which costs no exterior but lengthens the target chains. We omit 1,504 of the 10,922 deferred
gauges. The selection is optimal for this accounting: at a fixed saving, the moment is a sum of per-slot terms and
nonnegative interval costs along each target chain, so the best subset is an exact s-t minimum cut. That
formulation is Th0rgal's (PR #146, for the PR #97 word), reimplemented from its description and applied to our
witness (`certificates/round9/w2_omission.json`).
This gives `W = 31336`, `s = 2160160`, deficit `2024 = 2v - 3h(h-1)`, the coarse saving `a* = 467238028/10^12`, and, stopped at
`theta = 1/1000` with `a_old = 384599/10^10`, `a_bit = 466809249872/10^15`.

**Complex side: PR #144's paired cubes, unchanged.** This is icekylinx's `h = 24` complex inventory on eumemic's frozen
PR #117 addition DAG (`m = 72`, `W = 29937`, `s = 2153528`, largest child `60`), at PR #144's certified
`a_c = 4856569/10^10`. We rebuilt its frozen `p = 12` module role by role from PR #144's published data (the PR #117 DAG
and its 6,074 matching arcs). Every transition is nested, and the recount reproduces the published histogram
exactly (`independent/round9-complex-recount/`).

Together: `kappa = 4663738/10^10` (about 4.66374e-04), bound by the bit side. The minimum margin
`eps q = 4.6637382065e-04` is rounded down to the `10^-10` grid, as PR #144 states its own kappa.

**Proof interfaces.** As in earlier rounds, stage two is assumed to be the complement time-reversal of stage one.
This witness also rests on the following, which are stated here but not machine-checked:

- the three-stage cover lifting, the arbitrary-subspace frame theorem, the exterior rule (`3f` per retained gauge,
  none for an omitted one) and the stopped recurrence, from PRs #130 and #144;
- the frozen PR #117 DAG and PR #144's `p = 12` module. Both are rebuilt and recounted here from their published
  data, but the module's matching arcs are taken as given;
- the PR #97 word's frame theorem, inherited by our witness 2 through rounds seven and eight;
- the `2^-P 3^-K` precision grid of the semantic assembly. The star scatter of the cubes divides by 3, so the grid
  is no longer dyadic. This looks benign (exact linear arithmetic, about `log2 3` more bits per level), but only
  dyadic grids were audited before;
- block-factored frame representatives `T_{k(U+O)} = P_k (T_U (x) T_O) P_k^{-1}`. The cross-stage sharing tail has
  rank `3 dim sigma` with these. With generic representatives the sharing is still correct, but the tail rank rises.
  A factored choice always exists;
- the external row reserve of the finite bridge, which keeps PR #144's conservative bit constants
  (`ROW_RESERVE` in the certificate).

Checks (`make round9`, standard library, seconds; `make round9-heavy`, a machine with numpy and a few GB):

- `scripts/certificate_round9.py` reproduces PR #144 exactly and certifies ours, and `tests/test_round9.py` checks it;
- `independent/round9-bit-ledger/`: `ledger.py` builds a source-bound ledger of the modified word. Every register
  walks an explicit path of exact frames, and the ledger replays the complete forward F2 shear (1,051,969 events) and
  the reflected word. It also checks every gate's frame and that target chains are deleted subsequences, then
  rebuilds the child histogram from the events, with negative controls (`--ctl`). `cert.py` re-certifies the saving
  with independent two-sided ln/exp bounds, and first reproduces PR #144's published bit values;
- `independent/round9-complex-recount/`: `recount144.py` follows every one of the `R = 26417` physical roles of PR
  #144's complex module through its op list. Every transition must be nested on subspaces, and the result must
  equal the published histogram;
- `independent/round9-audit144/`: `cube_alg.py` (the cube identity `I = K + H + B`, `K^2 = I`, the root-group schedule),
  `cube_lit.py` (literal `Z[i][1/2]` replays with dirty scratch, gauge omission and controls) and `share_lit.py`
  (cross-stage sharing on one bank with a `3 dim sigma` tail, and controls);
- `independent/round9-cover-e2e/`: `cover_e2e.py` replays the three-stage cover exactly over `Q(i)` on our own local
  complex word, with controls.

The moment certificates (with the fallback), the stopped mix, the finite bridge and the 47-row assembly are also
checked in Lean 4's kernel (`lean/Round9.lean`, `make lean`).

Support compute: https://buymeacoffee.com/jainswapnil138

## Round eight (previous witness, `12612978530233/10^17`)

`scripts/certificate_round8.py` certifies both savings from frozen child-width histograms
(`certificates/round8/`) and assembles them with `scripts/certificate_round3.py` (`beta = 1/1000`, crude guard).

**Bit side: opposite bank orders (icekylinx, PR #104) on our deferred witness.** The schedule, frames and
payload word of our round-seven headline witness (deferred readouts, V leaves, lifted frames and late copies at
`h = 23`, `R = 27794`, with ikeboy's producer, PR #62, compiled by eumemic's joint frame compiler,
PR #57) are unchanged; only the address geometry of each interchange is recompiled:

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

This certifies the one-run saving `a* = 25367/(2*10^8)` and the stopped `a_b = 6338633/(5*10^10)`.

**Complex side: signed reclaim and completed-core sharing on eumemic's DAG.** Our `h = 24` two-stage complex word
(retained centres, carrier matching after icekylinx, PR #104, and the normal form for alternating
residuals, icekylinx, PR #24) is built on eumemic's frozen `h = 24` addition DAG (PR #117), read as data with every
support recomputed. Two mechanisms, both reimplemented from their descriptions, go on top:

- **Signed dependence reclamation** (jamesyc, PR #112). When an auxiliary is read for the last time,
  its signal is cleared by exact signed shears from slots that are still live (a duplicate, a difference with a
  user, a sum of its arguments, or a decomposition down the DAG to depth 6), and the cleared slot is reused by a
  later fresh slot. This takes the word from 28,705 to `R = 26874` auxiliaries.
- **Completed-core sharing** (an664, PR #128). The 2024 triples are split into 87 binary-orthonormal
  groups (83 of size 24, 4 of size 8; an664's partition, after Zhang and Ge). The cores of a group run
  consecutively on one bank of auxiliaries with dirty scratch, and the per-core wrap child is replaced by one
  separate fix-up child per group and stream. This cuts `W` to `12869228` at the same deficit `D = 1862080`.

This certifies `a_c = 3945999/31250000000` (`s_c = 7410813248`, largest child `529`).

Together: `kappa = 12612978530233/10^17` (about 1.26130e-04), bound by the complex side.

**Proof interfaces.** Like rounds six and seven, this witness assumes that stage two is the complement
time-reversal of stage one. It also rests on three interfaces that are stated, not machine-checked, here:

- the normal form for alternating residuals (icekylinx, PRs #24 and #104), which the complex replay uses for every
  alternating child;
- binary address adapters for the exterior (fix-up) children of the sharing, whose cost we take as subordinate
  (eumemic's PR #129 prices them at `64 m^2`);
- the row stock of the stopped bit interchange. In the stopped-product assembly it enters only the row-product gap
  (coefficient `10377`, degree `22000`, gap `20773/25`); `scripts/certificate_round3.py` has no such row, and
  `certificate_round8.extra_rows()` is where one would go.

Checks (`make round8`, standard library, minutes; `make round8-heavy`, a many-core machine with numpy and a C
compiler, hours):

- `scripts/certificate_round8.py`: both moment certificates with rigorous rational bounds (`ln(m/w)` from its
  series; `(m/w)^a <= 1/(1 - a ln(m/w))` on the bit side, `exp x <= 1 + x + x^2/2 + x^3/(6(1 - x/4))` on the complex
  side), the grid point checked to be the largest, the stopped mix, and the assembly; `tests/test_round8.py`;
- `independent/round8-oppbank/`: `bit_hist.py` rebuilds the bit histogram from the frozen round-seven schedule and
  first reproduces round seven's certified saving; `run_gate.sh` checks, on **every** charged residual of the
  witness, that one 60-bit random integer basis gives nonzero leading minors (fast primes, every zero re-checked
  mod `2^61 - 1`), samples the reversed-run Bruhat profile, and runs negative controls; `wrap.py` checks the outer
  wrapper and remainder pairing exactly;
- `independent/round8-coreshare/`: `complex_hist.py` rebuilds the complex histogram from the word, and
  `partition_check.py` checks the partition (exact cover, Gram `I`, even overlaps);
- `independent/round8-complex-gate/`: written without importing the construction (only `build_word.py` does, to
  pickle the word). `sreplay.py` replays both stages of the shared word exactly over `Q(i)` with dirty scratch, with
  every group exterior computed literally on 576-bit vectors, and checks every pair, every restored stream and the
  histogram (two seeds); `walk2.py` walks the actual op sequence; `excheck.py` checks the exteriors; `cert2.py`
  re-certifies `a_c` and `kappa` with independent ln/exp bounds; `endpoint.py` checks the endpoint correction.

The moment certificates, the stopped mix and the assembly are also checked in Lean 4's kernel (`lean/Round8.lean`,
`make lean`). The opposite-bank compile needs only the existence of the common basis; the all-edges check exhibits
one.

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
| Round-eight bit side (opposite bank orders, stopped interchange) | Same frozen schedule as round seven; histogram rebuilt from it; all-edges common-basis check with fast-prime passes re-checked mod `2^61 - 1`, Bruhat sample and negative controls (`independent/round8-oppbank/`); Lean kernel check of the moment, stopped mix and assembly. The existence of the common basis, the factorization lemma and the stopped recurrence follow icekylinx's written argument (PR #104) |
| Round-eight complex side (signed reclaim, completed-core sharing) | Histogram rebuilt from the word; exact `Q(i)` replay of the shared word with dirty scratch on two seeds, walk, exterior and certificate checks, written without importing the construction; partition checked (`independent/round8-coreshare/`, `independent/round8-complex-gate/`). The alternating normal form, exterior address adapters and the stopped interchange's row stock are proof interfaces |
| Round-nine bit side (witness 2 in PR #144's cover, gauge omission) | Source-bound ledger of the modified word: exact frame paths, complete forward and reflected F2 shear, target chains as deleted subsequences, histogram from events, negative controls (`independent/round9-bit-ledger/`); independent re-certification; Lean kernel check of the moment with the full fallback, the stopped mix and the 47-row assembly. The cover lifting, exterior rule and stopped recurrence follow icekylinx's written argument (PRs #130, #144) |
| Round-nine complex side (PR #144's paired cubes) | PR #144's published inventory; its `p = 12` module recounted role by role from the frozen PR #117 DAG and matching arcs (`independent/round9-complex-recount/`); cube identities, literal replays and cross-stage sharing checked with controls (`independent/round9-audit144/`); finite bridge and moment in Lean. The `2^-P 3^-K` precision grid and block-factored representatives are proof interfaces |
| Round-ten bit side | The JSON-only ledger replays the p = 12 bit word, and the same pipeline at p = 9 and 7, exactly over F2 with dirty scratch, forward and reflected, on exact rational frames. The mutation controls are rejected at every size (7 of 7 at p = 12), and the measured histogram certifies a* = 632988409/10^12 on the 10^-12 grid with the next point rejected. |
| Round-ten complex side | The p = 11 paired-cube word with the nested-prefix module, birth reuse and per-operation descent is re-checked from its frozen JSON alone: decoder identity, F2 chain legality, ledger recount and the 1e-12 certificate with the next point rejected. It is replayed exactly mod 2^61 - 1 with dirty scratch on two seeds, with mutated words and replay controls rejected, and the mechanism is replayed literally at p = 6 and 7. A second checker recomputes the finite guard and rejects the next grid point with a rigorous lower bound, and also checks a p = 7 word. |
| Round-eleven bit side | The JSON-only ledger replays the p = 12 word with its 24 terminal sinks, and the same recipe at p = 9 and 7, exactly over F2 with dirty scratch, forward and reflected, on exact rational frames. Eleven mutation controls are rejected, and the measured histogram certifies a* = 670100165/10^12 on the 10^-12 grid with the next point rejected. Lean kernel check of the moment, stopped mix and assembly. Partner-pair mixing and the sink lemma are written arguments (PRs #155, #168, #166) |
| Round-eleven complex side | The p = 11 recycled word is re-checked from its frozen JSON alone: decoder identity, sink eligibility, recycling chronology, F2 chain legality, ledger recount and the 10^-12 certificate with the next point rejected; a dirty replay mod 2^61 - 1 on two seeds; replay controls and three mutated words rejected; the same at p = 9 and 7. Zero-fresh recycling follows PR #184 Section 2a (written argument); the mix gate has no literal frame replay |
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
make round9      # round-nine certificate (with the PR #144 reproduction) and tests (seconds)
make round9-heavy  # round-nine ledger, complex recount, cover replay and PR #144 mechanism checks (numpy)
make round10     # round-ten certificate (with the PR #144 and round-nine regressions) and tests (seconds)
make round10-heavy  # round-ten replays, recounts and certificate checks of both sides (numpy)
make round11     # round-eleven certificate (with the PR #144, round-nine and round-ten regressions) and tests
make round11-heavy  # round-eleven bit ledger and complex recycling checker, with controls (numpy)
make lean        # Lean 4 kernel checks of rounds five to eleven
make notes       # PDF notes (tectonic)
```

The PDFs are in `artifacts/`.

## Attribution

Author: **Swapnil Jain**. Research, implementation and drafting were done with
assistance from Claude (Anthropic). This builds on OpenAI's manuscript and Douglas
Colkitt's framework. It also uses published results from open pull requests on
Colkitt's repository, cited in `NOTICE`; those are citations of prior work, not
collaborators. Licensed under Apache-2.0.
