# Closed ideas

This file lists approaches to raising kappa in OpenAI's *Integer multiplication below n log n* (problem #109 of
openai/math), as developed in Douglas Colkitt's [integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds),
that we tried and found do not improve kappa. Each row has a stated reason. We list it so that other researchers and
automated searches do not spend time on these routes again.

**Contributing.** Please check this list before starting on a route, so effort is not duplicated. If you find a
mistake in a closure, an idea that escapes a stated hypothesis, or a closed route of your own, open an issue or a
pull request with the argument or certificate. New closures and corrections are welcome and will be credited.

**Scope.** A closure holds only inside its stated model or hypotheses: the frame calculus and role accounting of the
upstream framework, a stated class of words, producers, modules or covers, or one fixed word. Nothing here is a lower
bound on the complexity of integer multiplication, and nothing here caps kappa in general. The hypotheses of a closure
show where an escape would have to break the model. An idea that breaks a stated hypothesis is not covered by the row.

**Status labels.** [proof]: a complete argument. [cert]: exact computer certificate (exact rational or finite-field
elimination, or exhaustive enumeration). [exp]: a decisive experiment, or a floating-point computation named as such.
[arg]: a short argument, convincing but not written as a full proof. Rows marked P1 or P2 are proved, with their
hypotheses, in our papers:

- P1: *Ceilings and Impossibility Results for Finite Witnesses of Sub-n log n Integer Multiplication*,
  [doi:10.5281/zenodo.23265442](https://doi.org/10.5281/zenodo.23265442);
- P2: *No-Go Results for Finite Witnesses of Sub-n log n Integer Multiplication*,
  [doi:10.5281/zenodo.23266393](https://doi.org/10.5281/zenodo.23266393).

Here we give only the statement and model; the proofs, and the programs behind every computation, are in the papers
and their reproduction archives. Both papers were checked by adversarial review passes of our own. That is not
independent peer review.

**Notation.** As in the framework: h is the dimension of the port space, m the cover dimension, v the number of ports, W the
number of roles, s the total child rank, D = Wm - s the deficit, l the centre loss and a* the certified saving of a
side (the root of its certificate). A *frame* is the subspace at which a gate acts; each register has a *chain* of
gate frames.

**Accounting hypotheses.** Many closures below assume some of the following (P2, Section 1). A closure that uses one
says nothing about a construction that breaks it.

- H-role: each auxiliary role is one register with its own chain of frames, a fresh request costs one role, and every
  auxiliary role contributes total rank exactly m per vertex.
- H-src: original-source registers are not used as a resource to erase responses.
- H-pair: reuse is a pairwise match between a retiring register and a new use (birth reuse, absorption, carriers,
  sinks), not a collective stock.
- H-frame: the gate frames of every register form a nested chain, increasing in time, with one child per piece.
- H-id: the delivered identity splits as I = C + V, a centre part C delivered through copied centres and a visible
  part V delivered by the producer; for paired cubes this is icekylinx's identity H + K + B = I.

icekylinx's source-assisted common-frame compression (PR #184) breaks H-role and H-src, and dead-copy recycling and
collective entrance banks break H-pair. Several of our own earlier closures had assumed these without saying so,
which is why every row now names its model.

## 1. Address fields, identities and sources

| idea | why it fails | status | credit |
|---|---|---|---|
| Odd-characteristic address fields F_p in place of F_2 | Statement: for complete families of w-subsets of [h] with a polynomial star scatter of degree r realised by the copied r-stars, the loss is C(h,r)(h - r) for every p, the largest admissible port weight is 2r + 1 for p = 2 and at most floor(3r/2) + 1 for odd p, so the binary field gives the most ports per unit of star loss at every degree. Model: H-id, complete 0/1 families with w <= h/2, polynomial star scatters | [proof] P2 | |
| Address fields F_{2^k} | Every construction in the frame calculus over F_{2^k} is, after restriction of scalars, one over F_2 with k-dimensional sources, every child rank and m multiplied by k and the roles unchanged, so a* is the same. Model: frame calculus | [proof] P2 | |
| Two-dimensional sources (j = 2) by tensoring with F_2^2 | The diagonal lift doubles every child rank and m with the same W: exactly the same a*. Model: frame calculus over F_2 | [proof] P2 | |
| Two-dimensional coordinate families (the S_h-orbit of one 2-dimensional label) | No degree-one decoder for support size K >= 3 and h >= 2K - 1; degree-two decoders are infeasible for every shape at K = 4, 5, 6 (h >= 2K). Model: H-id | [proof] + [cert] P2 | |
| The complete triple design with j = 2 sources | Every linear realisation on n >= 7 atoms, with any symmetric form, needs h >= 2n (h >= 16 at n = 9), twice the dimension of one-dimensional sources. Model: linear realisations, atom-star decoder, H-id | [proof] P2 | |

## 2. Centres

| idea | why it fails | status | credit |
|---|---|---|---|
| Star centres with more ports per unit of loss | For every family of 3-subsets of [h], h >= 4, with star centres, v/l <= (h - 2)/6, so D/v <= 2 - 18/(h - 2); for h >= 5 equality holds only for the complete triple family. Model: weight-3 ports | [proof] P1 | |
| Replace the star centres of paired cubes by another centre family | Each of these families has l >= h(h - 2), which the h stars attain (l = 440 at p = 11): every family with exactly h centres; every family whose frames all have dimension at least h - 2; every family containing h - 1 of the stars (a star cannot be replaced by centres of smaller total loss). Model: paired cubes with p >= 5, frame-0 copy model; base case by computation | [proof] + [cert] P1 | paired cubes: icekylinx (PR #144) |
| Fewer distinct bit centre values | The bit minrank is mr(h) = h for h >= 4 except h = 7 (mr(7) <= 6), so every decoding family of bit centres at h = 23 has at least 23 distinct centre values. Cases h = 5, 6 by exhaustive computation | [proof] + [cert] P1 | |
| Parity bases as bit centres | Every parity basis has l >= h(h - 1) - 1; at h = 23, l >= 505 against 506 for the stars, so at most one unit is saved. Model: odd h >= 7, h != 9, descent model, per-piece accounting | [proof] P1 | |
| Weight-five ports served by producers that hold class sums | A centre part of degree at most 2 is forced to f(j) = (j - 1)(j - 3)/8, and a producer that holds every class sum and every centre in an auxiliary register needs rho = R/v >= 16 + C(h,2)/C(h,5) (16.006 at h = 24). Model: H-role, H-id, h >= 11, producers that hold every class sum | [proof] P2 | |

## 3. Producers, meets and buses

| idea | why it fails | status | credit |
|---|---|---|---|
| Producers without auxiliary registers (direct chain meets) | Meet capacity: with one source chain per source and one target chain per target, the share of star demand that can be met is at most 1/4 + O(1/h) for the triple and cube families (43.5% at h = 24). Model: one chain per port | [proof] P1 | |
| Pad the space to escape the meet bound | Intersecting padded chains with the original space keeps every met pair, so the meet-capacity bound holds with the capacities of the unpadded space. Model: orthogonal padding, nested chains | [proof] P2 | padded covers: eumemic (PR #137) |
| Two-frame buses | Every pair that a two-frame bus serves is already met; replacing all bus gates by direct gates gives a legal word with the same map and no smaller a*. Model: bus words, per-piece accounting, H-role, no fallback term in the certificate | [proof] P1 | |
| Junk riding: one call carrying the residuals of several wires | A single call (any invertible linear map on a union of blocks) that realises the targets of several wires has region at least the sum of their separate regions; scratch wires contribute 0. Model: a single call | [proof] P1 | |

## 4. Addition modules

| idea | why it fails | status | credit |
|---|---|---|---|
| All-but-one DAGs with a fork whose consumers have nested frames | Antichain theorem: in every disjoint-support addition DAG for an all-but-one instance, the frames of the consumers of any node are pairwise incomparable. Model: the annihilator condition of the all-but-one module | [proof] P1 | modules: icekylinx (PR #144), eumemic (PR #168) |
| Fewer copies per input in t-Kneser modules | Hitting-set bound: every input costs at least t copies and the instance at least t C(n,t). The cyclic Fibonacci all-but-one module attains the floor at n = 9, so in this model its copy count cannot be lowered. Model: nesting, visibility, the annihilator hypothesis and the relaxed-arc model; one property of the words (Fact F) is verified by exhaustive computation at p = 7 and 11 but not proved | [proof] + [cert] P1 | Fibonacci strips: eumemic (PR #152); interval strips and an exact integer program for t = 1: ikeboy (PR #62) |
| Labels of a cube-local circuit held by one role | Every label costs at least two copies. Model: as in the previous row, with part (F2) of Fact F | [proof] + [cert] P1 | cube-local circuit L1: eumemic (PR #168) |

## 5. Motifs and lifting

| idea | why it fails | status | credit |
|---|---|---|---|
| Motifs whose frames and code have the same characteristic | If the frames lie in Mat_m(K) and the code is affine over a field of the same characteristic, the total cost is at least Wm, so s >= Wm and there is no deficit; a defect of rank k gains at most 2mk. A positive deficit needs mixed characteristic. Model: motifs with an affine code, copy edges allowed | [proof] P1 | |
| A motif with large saving in characteristic 0 or 2 | Every compiled motif has at least 2W children, so a <= 1 - log_m 2, strictly for m >= 3. Model: shear form, any bijective gates | [proof] P1 | |
| Segments of spectral radius at most 1 (for example +-idempotent segments) | Trace lemma: then D <= 0; in general s/(Wm) >= 1/max spr(A_e). Model: frames over a subfield of C | [proof] P1 | |
| Lift the rank-one Cayley deck of GL3(F3) to characteristic 0 | No lift with L(-J) = -L(J) and rank-one edge differences over the ring of integers of any finite extension of Q3, for the reflection deck and the full deck; every realisation over Z/3^b, b >= 2, has an edge whose frame change costs nearly two fields. Model: rigid edge sets obstructed mod 9, sign symmetry included; rigidity and obstruction by exact computation | [proof] + [cert] P1 | |
| Lift affine-plane configurations from F_q to Z/q^2 | For q in {3, 5, 7}, AG(2,q), and AG(2,q) minus a point, have no lift to (Z/q^2)^2 that keeps every line collinear. Model: as stated | [cert] P1 | |

## 6. Cost model and recursion interface

| idea | why it fails | status | credit |
|---|---|---|---|
| Non-nested or non-CSS labels on a label path | A path of total distance m between transversal Lagrangians uses only CSS labels with nested subspaces; any path with a non-CSS label or a non-nested step costs at least m + 1. Model: Lagrangian label calculus | [proof] P1 | |
| Several children per edge, or batching one child across roles | The children realising an edge have total width at least the length of its image, so one child per edge is optimal for every concave increasing cost, and cross-role batching gains exactly 0. Model: the recursion interface, costs depending only on width and linear in volume | [proof] P1 | |
| Mix motifs across nodes, levels or widths (stopped and multi-supplier recursions included) | Every such procedure costs at least that of the best single motif in the set. Model: the recursion interface | [proof] P1 | a special case for a fixed family of profiles: JosephDemarest (PR #83) |
| Unequal field widths, or q-adic frames | With block-diagonal frames over width classes, the root is at most that of the best class; q-adic frames are no better than the best F_q rank motif read off their graded pieces. Model: as stated, length-of-image costs | [proof] P1 | |
| Cross-wire digit shuffles through position maps | Position-map compiles pay the full rank: the expected charge per column is at least phi(r)(1 - o(1)) as the field width grows. Model: position maps | [proof] P1 | |
| Tighter exponential majorant in the stopped bit certificate (the round-seven refinement of PR #2, degree-eight Taylor, applied to rounds ten and eleven) | The certificate already uses E3(x) = 1 + x + x^2/(2(1 - x/3)). Recomputing the bit moment with exp and ln to 80 digits, the true root exceeds the certified a* by 2.3e-13 in round eleven and 6.0e-13 in round ten, below the 1e-12 grid of a*, so no majorant for exp can move the certified saving by a grid step; the fallback term costs about 1.7e-12 of a*. Model: the frozen bit inventories of rounds ten and eleven, their 1e-12 grid; floating point named as such (`scripts/check_majorant_slack.py`) | [exp] | the Taylor refinement: SovereignSteak (PR #2 here) |
| Re-price gauge reads that split a target chain | The split charge is exact: D is unchanged and the certificate sum rises by g_a(d) + g_a(L - d) - g_a(L) > 0. Model: H-frame, a gate inserted into a chain, the other register unchanged | [proof] P2 | |
| Children that span two visits | With a gate at the exit frame of one visit and at the entrance frame of the next, the two pieces are distinct children. Model: H-frame | [proof] P2 | |
| Source-assisted common-frame compression across stage handoffs | Stock gain 0: neither the dirty rows R nor the roles W change, and under H-role neither does D; a claimed positive gain must skip a handoff restore or carry a non-zero-fresh row across a handoff. Model: every handoff a rank-0 step with restore, every crossing row zero-fresh | [proof] P1 | the compression and its stock lemma: icekylinx (PR #184) |

## 7. Covers, stages and schedules

| idea | why it fails | status | credit |
|---|---|---|---|
| Swap two banks with fewer layers, or let bank pairs share three layers | No layered call with at most two layers swaps the banks, and every three-layer call that does has layer rank at least 3n for n port pairs, so sharing three layers gains nothing in rank. Model: layered calls with at most three layers, any field | [proof] P2 | |
| Stages that deliver a non-identity coupling repaired by free maps | When free maps act only on the second bank before a call and on the first bank after it, every stage's coefficient block is +-I. Model: three block shears, free maps only at the two ends | [proof] P2 | three-stage Cayley covers: icekylinx (PR #130) |
| Cover-level changes to a fixed bit word: padding m, orders, sharing, tails, exteriors, finishing, omission | For the bit word of our round-nine witness in icekylinx's paired-cube cover, every cover-level arrangement has a* <= 4.8429e-4 < 2^-11 (4.7895e-4 for m >= 69). Model: H-role, H-src, H-pair, H-frame, the local multiset of the word kept, deficit at most 2v - 3l | [proof] + [exp] (floating-point minimum cut) P2 | a ceiling over all frame layouts of one fixed word: DaysSky (PR #192) |
| Cover, sharing, padding or exact-compile changes to the PR #144 complex word that keep its visits | The visit ceiling for that word is 9.70e-4 at the least admissible m = 70 (9.45e-4 at m = 72), below 2^-10. Model: that word only; its compiled children stay within visits, its visits are monotone, H-role, H-src, H-frame | [proof] + [exp] P1 | the word: icekylinx (PR #144); the inequality a* <= D/Lambda: DaysSky (PR #192) |
| Fuse consecutive calls: apply the next call's twiddle before the previous call's Hadamard layer finishes on a register | A twiddle conjugated by a Hadamard layer is monomial only if it is a scaled character, so two calls cannot interleave on a register at a top-level boundary. Model: the current call interface, hypotheses (I1)-(I3) of P2 | [proof] P2 | |
| Run the paired cores of a cover in lockstep on one shared auxiliary stream | The merged operator is right but the schedule is not: one core's frame steps transform the data the other has already added. Exact replays with dirty scratch: 64 of 64 target entries wrong on a small literal instance, 3,136 of 3,136 data pairs wrong on an h = 8 complex word whose consecutive schedule passes. Only merges at the boundary between consecutive cores survive, about 0.05% on our words | [exp] P2 | lockstep cores: ikeboy, who re-checked the replay, withdrew the claim and gave the sound consecutive variant; replay confirmed by GamingPuzzled |
| Run whole stage-1 and stage-2 invocations in lockstep to reuse centre copies | The source connectors alone cost at least 2(v - 1) = 4046, against at most 2l = 1104 saved. Model: whole invocations, PR #130's label assignment | [arg] P2 | label assignment: icekylinx (PR #130) |

## 8. Scratch and cleanup

| idea | why it fails | status | credit |
|---|---|---|---|
| Conditionally clean ancillae in linear words | A linear map that fixes every vector outside a proper affine subspace is the identity, except over F_2, where it is a parity-controlled transvection: linear words have no branches beyond parity-controlled gates. Model: linear words | [proof] P2 | the notion: Khattar and Gidney |
| Clean scratch inside the W roles | Conservation: at every time, the auxiliary registers free of z are at most as many as the data registers that involve z. Model: linear words; its cost reading uses H-role | [proof] P2 | |
| Cleanups that end a role below the join of its carriers' frames | Carrier lemma: the last gate frame of an auxiliary register contains every frame at which another register carries its dirty symbol. Model: transparent linear words with nested chains, any field | [proof] P2 | |
| Full low-frame source cleanup of a source's demand | If the lowest entry node holds a second source s1, the targets reached from S through source-cleaned roles number at most mu(S) = max over s != S of the shared demand of S and s. Model: H-frame (nesting only), source-cleaned roles, a second source in the lowest entry node | [proof] P2 | |

## 9. Earlier stacks

| idea | why it fails | status | credit |
|---|---|---|---|
| Tune the rows of the stack after the compressed network and fast Gaussian resampling | Within that row set (unchanged layout, CRT reversal by pairwise width-l swaps, the Gaussian precision budget of about 8d^2 fitted in O(log n)), every witness has kappa < a_b/2: the butterflies save on the number of axes d, the chunk and axis layouts save on the axis width l, and d l ~ log n. This is algebra on those rows, not a limit on new rows; see "The ceiling and how this gets past it" in the README | [proof] | the rows: eumemic (PRs #3, #5) |
| Longer digits alone | Digits of (log n)^(1+x) bits separate the precision budget from the address length, but alone they leave the Gaussian row capping eps < 1/2; they help only in combination with the other changes listed in the README (segmented inverse, one chunk per axis, fine-bit exposure, permuted selected-bit swaps) | [arg] | |
