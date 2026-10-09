# A sharper certificate for the frozen round-seven bit network

Basis: `Swapnil-jain/integer-mult-kappa` commit
`741e7aa078392553815df7926ee17ac5e25a8c38`.

The round-seven staircase bit histogram already in this repository admits
the stronger rational moment certificate

```
bit saving a = 63983013240044 / 10^18
conditional kappa = 799736495947 / 12500000000000000
                  = 0.00006397891967576.
```

This raises the pinned headline `0.00006396990758892` by exactly
`225302171/25000000000000000`, about 0.014088%. The schedule, histogram,
complex saving, assembly selection rule and all construction assumptions are
unchanged. This is a tighter certificate for the existing conditional
construction. Its stage-two, analytic, precision, and tape implementation
dependencies remain as stated in the original manuscript and notes.

## Why the bound improves

For the unchanged child histogram, put `m=529`, `W=108516254` and

```
M(a) = sum_t [n_t t/(m W)] exp(a log(m/t)).
```

The original certificate used `exp(u) <= 1/(1-u)`. We enclose `log(m/t)`
between exact rationals and use the positive Taylor series through degree
eight, with the geometric upper tail

```
exp(u) <= sum_{j=0}^8 u^j/j! + (u^9/9!)/(1-u/10)
```

for `0 <= u < 1/1000`. Directed rounding yields strict rational checks
`M(63983013240044/10^18) < 1` and
`M(63983013240045/10^18) > 1`. All child widths satisfy `0<t<m`, so the
moment is strictly increasing. The interval identifies the last feasible
point on this grid for this histogram. `ENCLOSURE_PROOF.md` gives the
independent degree-twelve derivation and its ablation of the old majorant.
The Taylor technique itself is established; CrocSwap PR65 also uses a
tighter exponential enclosure.

| Variant | Bit saving | Conditional kappa |
| --- | ---: | ---: |
| Published reciprocal, `10^-9` grid | .000063974 | .00006396990758892 |
| Same reciprocal, `10^-18` grid | .000063974328963846 | .00006397023651068 |
| Taylor, `10^-18` grid | .000063983013240044 | .00006397891967576 |

Independent exact assembly arithmetic preserves `beta=1/1000`, the complex
saving `36926111/(5*10^11)`, the original epsilon and kappa grids, and
the guard `C1=14694`, `x=14695`. Simultaneous butterflies remain binding.
The refinement changes no source outside this directory.

## Reproduce

Use Python 3.11+ with assertions enabled, and run from this repository root.
Output paths must be new. The scripts use only the standard library. Start
with the source pins:

```sh
python3 research/round7-moment-refinement/check_sources.py
ROUND7_REFINED_OUT=$(mktemp -d)
python3 research/round7-moment-refinement/refine_moment.py \
  --output "$ROUND7_REFINED_OUT/primary"
python3 research/round7-moment-refinement/independent_moment.py \
  --case accepted=exp:accept:63983013240044/1000000000000000000 \
  --case adjacent=exp:reject:63983013240045/1000000000000000000 \
  --output "$ROUND7_REFINED_OUT/independent.json"
python3 research/round7-moment-refinement/independent_assembly.py \
  --claim "$ROUND7_REFINED_OUT/primary/certificate.json" \
  --out "$ROUND7_REFINED_OUT/assembly.json"
```

The primary checker uses a base-two logarithm reduction, 40 atanh terms and
degree-eight Taylor expansion. The independent checker uses base `3/2`, 48
atanh terms, degree twelve and its own directed rounding. The assembly
checker computes the guard ceiling with rational bounds and checks every
strict constraint and cost margin, including the old reciprocal thresholds.
No floating-point comparison decides acceptance.

With Lean 4.31.0, the following compiles two finite rational statements
using the entire existing `Round7.bitHist`:

```sh
lean -o "$ROUND7_REFINED_OUT/Round7.olean" lean/Round7.lean
LEAN_PATH="$ROUND7_REFINED_OUT" lean \
  -o "$ROUND7_REFINED_OUT/Round13Moment.olean" \
  research/round7-moment-refinement/Round13Moment.lean
```

`#print axioms` reports Lean's standard `propext` axiom for both statements;
there are no `sorry` declarations or project-specific axioms. The statements
check rational predicates, not the real-series tail lemma or the complete
multiplication theorem. The original physical word,
frame and histogram checks remain the construction evidence. No native
performance claim follows from this certificate.

Credit: Swapnil Jain's round-seven construction and its attributed OpenAI and
CrocSwap dependencies remain the basis. This numerical refinement and
independent checks were prepared in the SovereignSteak research workspace
with OpenAI Codex assistance. Existing notices and licenses are unchanged.
