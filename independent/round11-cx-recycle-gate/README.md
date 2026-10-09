## Complex side: dead-copy recycling (a_c = 67147467/10^11 = 6.7147467e-4, p = 11)

### Base word

The base is the paired-cube complex word on #168 v4's modules (eumemic, data, sha-pinned). It uses the configurable
cube-local circuit in our C1 configuration (localmotif), our relaxed arcs (reverse scan), per-operation descent, birth
reuse with late reads, and 42 terminal-output deletions (#166/#176).

Its own certificate is a_c = 662034051/10^12, and `check_word2.py` passes on it.

### Mechanism

The base word holds 428 pairs of registers that each retire holding the same node value: two copies of one value,
neither of which is touched again. For such a pair (pivot P, control Q):
1. One gate `P -= Q` runs at a frame V that contains the last frames of both P and Q. After it, P holds zero signal
   (no x-dependence) and only arbitrary dirty contents.
2. A fresh birth B is placed on P's register. B is a role that starts at frame 0, carries no gauge, is not a source
   and is not a sink, and its first op has frame V. The mix runs immediately before that op, and B is renamed to P.
3. No read is added. B's start content is a linear function of initial dirty values. The frame-0 compensation read
   of every non-gauged role is recomputed backward over the new transcript, and the mix enters it as an ordinary
   gate with coefficients (1, -1). Copies whose dirty part depends on a late-read (gauged) value are excluded, so
   every gauge read is unchanged.

This is Section 2a of icekylinx's #184 (zero-fresh recycling through a containing frame), used with E = 0: no source
controls. The kernel comes from duplicated copies, not from source erasure.

### Ledger

Per host, with u0, u1 the dims of the copies' last frames and V the dim of the mix frame, every term times 3 stages:

| register | removed | added |
|---|---|---|
| pivot | [h - u0] (its exit to F) and [V] (B's first step from 0) | [V - u0] |
| control | [h - u1] | [V - u1] and [h - V] (it must sit at V for the gate) |

W drops by 1 per host. The added widths sum to the removed ones minus h per stage, so the rank mass drops by m and
D = Wm - s is unchanged (deficit 2v - 3 loss = 1320 is checked).

The sign test is dN < dW at the current a. Hosts are chosen by a max-weight bipartite matching between copy pairs and
births. Each copy is used once, and pairs are disjoint.

On the p = 11 word: 386 hosts, W 13308 -> 12922. Shapes (u0, u1, V): (2,3,4) x205, (3,3,4) x138, (2,3,8) x30.

### Results (exact; certificate_round10's complex_side gives the same a_c)

| p | base a_c | recycled a_c | hosts | change |
|---|---|---|---|---|
| 11 | 662034051/10^12 | 67147467/10^11 | 386 | +1.426% |
| 9 | 133579197/2.5e11 | 270748449/5e11 | 183 | +1.344% |
| 7 | 104238523/5e11 | 105042291/5e11 | 45 | +0.771% |

### Validation

`independent/round11-cx-recycle-gate/check_word2_recycle.py` (`heavy.d/21-cx-recycle.sh`) checks:
- the decoder identity;
- sink eligibility;
- recycling chronology: both copies retired before the mix, the control never touched again, copies are not
  root or gauge roles, hosts disjoint, and the mixed copies hold equal values (asserted while values are tracked);
- F2 nesting of every slot chain (including both copies' climbs to V) and of every target chain;
- the ledger, recounted from chains, against the exporter's formula claim;
- the deficit;
- the 1e-12 certificate, with the next point rejected;
- a mod-P replay with arbitrary dirty slots and targets, on 2 seeds.

Controls that must be rejected:
- all of check_word2's birth and terminal controls;
- mix skipped and mix sign flipped (replay controls);
- three mutated words: a control holding a different value, a birth before the copies retire, and a tampered claim.

All pass at p = 7, 9 and 11. The p = 7 word rebuilt from `independent/round11-cx-recycle-build`
(`bash verify2.sh cfgC1sm.json 7 fwd; python3 hpair3.py ...`) is byte-identical in content to the released one.

### Not covered

- There is no literal frame replay of the mix gate: frames are checked by F2 nesting of the chains, as for every
  other complex op.
- Hosted roles that later die are not returned to the stock. This is conservative.
- Larger zero-signal supports (three or more retired rows) are not used. The frame-compatible ceiling is 1108
  births (bound only).
