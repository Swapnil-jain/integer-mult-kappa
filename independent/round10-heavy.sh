#!/bin/bash
# Round-ten heavy checks (make round10-heavy). Assembled from the release gates' snippets; each one exits
# nonzero on failure. R is the repository root.
set -e
R=$(cd "$(dirname "$0")/.." && pwd)

# ---- 20-bit-gate.sh
(
# bit side, bitcube-port gate (paired-cube bit word, p = 12, with NDS, min-cut omission, descent and birth reuse):
# the independent JSON-only ledger (exact F2 replay with dirty scratch, forward and reflected, exact frames over Q,
# G-nondegeneracy, event histogram, rational certificate with the next grid point rejected) on the real word and on
# the same pipeline at p = 9; seven mutation controls at p = 9 must be rejected; the release inventory must equal the
# replay's measurement. About 2 GB and 5 minutes.
T=/tmp/round10-bit; mkdir -p $T
cd "$R/independent/round10-bit-gate"
python3 ledger.py "$R/certificates/round10/bit_word_p12.json.gz" --expect 632988409/1000000000000 --out $T/led_p12.json | grep -q '] ALL PASS' || { echo "bit ledger p=12 FAILED"; exit 1; }
python3 ledger.py "$R/certificates/round10/bit_word_p9.json.gz" --expect 1426101/3906250000 --out $T/led_p9.json | grep -q '] ALL PASS' || { echo "bit ledger p=9 FAILED"; exit 1; }
for c in drop_prelude_read wrong_coefficient desc_shrink late_readout below_level broken_mix early_birth; do
  python3 ledger.py "$R/certificates/round10/bit_word_p9.json.gz" --ctl $c | grep -q 'ALL FAIL' || { echo "bit control $c NOT rejected"; exit 1; }
done
python3 invcheck.py $T/led_p12.json "$R/certificates/round10/bit_inventory_p12.json.gz" | grep -q 'ALL PASS' || { echo "bit inventory != replay"; exit 1; }
echo "bit gate checks: PASS"
)

# ---- 20-cx-prefix.sh
(
# complex side, cube-prefix word (p = 11, a_c 15402419/25000000000): the independent checker on the frozen word
# (decoder identity mod P, F2 chain legality, ledger recount = claim, 1e-12 certificate with the next point rejected,
# aliased dirty-scratch replay on two seeds, three replay controls), four mutated words that must be rejected, and the
# literal small-p replay of the mechanism (nested-prefix module + merged reads + birth reuse) with its controls.
# About 1 GB and a few minutes.
T=/tmp/round10-cx-prefix; rm -rf $T; mkdir -p $T
cd "$R/independent/round10-cx-prefix-gate"
W="$R/certificates/round10/cx_prefix_word_p11.json.gz"
python3 check_word.py "$W" | grep -q 'RESULT check_word pf_G37_rx p=11 a_c=15402419/25000000000: ALL PASS'
python3 mutate_word.py "$W" $T > /dev/null
for c in sign frame claim pair; do
  if python3 check_word.py $T/mut_$c.json.gz > $T/mut_$c.out 2>&1; then echo "cx-prefix control $c ACCEPTED"; exit 1; fi
done
python3 lit_prefix.py 6 lit 1 gauge desc birth prefix merge | grep -q 'control=-: ALL PASS'
python3 lit_prefix.py 7 scal 1 gauge desc birth prefix merge | grep -q 'control=-: ALL PASS'
for c in dropbirth birthlate birthframe; do
  python3 lit_prefix.py 6 lit 3 gauge desc birth prefix merge control=$c | grep -q 'SOME FAIL' || { echo "cx-prefix literal control $c passed"; exit 1; }
done
echo "cx-prefix checks: PASS"
)

# ---- 40-standing-gate.sh
(
# standing gate supplement. Complex: the word checker on a smaller frozen word (p = 7) of the same pipeline, the
# finite semantic guard (#144's formulas) on the real and the small word, and the certificate's next grid point
# rejected by a LOWER bound of the moment. Standard library, about a minute.
cd "$R/independent/round10-standing-gate"
T=/tmp/round10-standing; rm -rf $T; mkdir -p $T
CW="$R/independent/round10-cx-prefix-gate/check_word.py"
python3 "$CW" "$R/certificates/round10/standing_cx_word_small_p7.json.gz" | tail -1 | grep -q 'ALL PASS' || { echo "small complex word FAILED"; exit 1; }
for w in "$R"/certificates/round10/cx_prefix_word_p*.json.gz "$R/certificates/round10/standing_cx_word_small_p7.json.gz"; do
  python3 cxguard.py "$w" | grep -q 'ALL PASS' || { echo "cxguard $w FAILED"; exit 1; }
  python3 -c "import gzip, json, sys; J = json.load(gzip.open(sys.argv[1], 'rt')); json.dump(dict(C=J['claim']['C'], W=J['claim']['W'], m=J['m'], a=J['claim']['a_c']), open(sys.argv[2], 'w'))" "$w" $T/c.json
  python3 twosided.py $T/c.json $(python3 -c "import json; d = json.load(open('$T/c.json')); print(d['m'], d['a'])") | grep -q 'ALL PASS' || { echo "two-sided certificate $w FAILED"; exit 1; }
done
echo "standing-gate checks: PASS"
)
echo "round10-heavy: ALL PASS"
