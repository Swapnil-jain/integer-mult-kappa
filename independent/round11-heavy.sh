#!/bin/bash
# Round-eleven heavy checks (make round11-heavy). Assembled from the release gates' snippets; each one exits
# nonzero on failure. R is the repository root.
set -e
R=$(cd "$(dirname "$0")/.." && pwd)

# ---- 20-bit-gate.sh
(
# bit side, round-11 gate (paired-cube bit word, p = 12): #168 v4's pair module (data), our cyclic all-but-one module,
# L1 association, merged reads, NDS, min cut, node/group/per-op descent, birth reuse, 24 terminal sinks.
# The checker is the independent JSON-only ledger:
#   - exact F2 replay with dirty scratch, forward and reflected;
#   - exact frames over Q and G-nondegeneracy;
#   - absorbed-sink events;
#   - the event histogram;
#   - the rational certificate, with the next grid point rejected.
# It runs on the real word and on the same configuration at p = 9 and p = 7. Eleven mutation controls must be
# rejected at p = 7 (the four sink controls need sinks; p = 9 has none). The release inventory must equal the
# p = 12 replay's measurement. About 3 GB and 5 minutes.
T=/tmp/round11-bit; mkdir -p $T
cd "$R/independent/round11-bit-gate"
python3 ledger.py "$R/certificates/round11/bit_word_p12.json.gz" --expect 134020033/200000000000 --out $T/led_p12.json | grep -q '] ALL PASS' || { echo "bit ledger p=12 FAILED"; exit 1; }
python3 ledger.py "$R/certificates/round11/bit_word_p9.json.gz" --expect 2289691/6250000000 --out $T/led_p9.json | grep -q '] ALL PASS' || { echo "bit ledger p=9 FAILED"; exit 1; }
python3 ledger.py "$R/certificates/round11/bit_word_p7.json.gz" --expect 91020503/500000000000 --out $T/led_p7.json | grep -q '] ALL PASS' || { echo "bit ledger p=7 FAILED"; exit 1; }
for c in drop_prelude_read wrong_coefficient desc_shrink late_readout below_level broken_mix early_birth drop_post drop_pre late_write stray_write; do
  python3 ledger.py "$R/certificates/round11/bit_word_p7.json.gz" --ctl $c | grep -q 'ALL FAIL' || { echo "bit control $c NOT rejected"; exit 1; }
done
python3 invcheck.py $T/led_p12.json "$R/certificates/round11/bit_inventory_p12.json.gz" | grep -q 'ALL PASS' || { echo "bit inventory != replay"; exit 1; }
echo "bit gate checks: PASS"
)

# ---- 21-cx-recycle.sh
(
# complex side, dead-copy recycling (p = 11, a_c 67147467/100000000000): check_word2_recycle.py is check_word2.py
# extended by recycling hosts. Two retired copies of one value are differenced (pivot -= control) at the first-op
# frame of a fresh birth, which is then born on the pivot with no read.
# It checks:
#   - the decoder identity;
#   - sink eligibility;
#   - the recycling chronology (copies retired, control never touched again, no root/gauge copies, disjoint hosts);
#   - F2 legality of every role chain incl. both copies' climbs to the mix frame and of every target chain;
#   - the ledger recounted from chains == the exporter's claim;
#   - deficit 2v - 3 loss;
#   - the 1e-12 certificate, next point rejected;
#   - a mod-P dirty-slot/dirty-target replay on two seeds, with the frame-0 compensation recomputed from the new
#     transcript.
# Controls: all of check_word2's birth/terminal controls, plus two replay controls (mix skipped, mix sign flipped)
# and three mutated words (unequal copies, birth before the copies retire, tampered claim) that must be rejected.
# The same runs on the p = 7 and p = 9 words from the same recipe.
T=/tmp/round11-cx-recycle; rm -rf $T; mkdir -p $T
cd "$R/independent/round11-cx-recycle-gate"
W="$R/certificates/round11/cx_recycle_word_p11.json.gz"
python3 check_word2_recycle.py "$W" | grep -q 'p=11 a_c=67147467/100000000000 sinks=42: ALL PASS' || { echo "recycle p=11 FAILED"; exit 1; }
python3 check_word2_recycle.py "$R/certificates/round11/cx_recycle_small_word_p9.json.gz" | grep -q 'p=9 a_c=270748449/500000000000 sinks=26: ALL PASS' || { echo "recycle p=9 FAILED"; exit 1; }
python3 check_word2_recycle.py "$R/certificates/round11/cx_recycle_small_word_p7.json.gz" | grep -q 'p=7 a_c=105042291/500000000000 sinks=27: ALL PASS' || { echo "recycle p=7 FAILED"; exit 1; }
for w in "$W" "$R/certificates/round11/cx_recycle_small_word_p7.json.gz"; do
  python3 mutate_hosts.py "$w" $T > /dev/null
  for c in unequal early claim; do
    if python3 check_word2_recycle.py $T/mut_$c.json.gz > $T/mut_$c.out 2>&1; then echo "recycle control $c ACCEPTED"; exit 1; fi
  done
done
echo "cx-recycle round11 checks: PASS"
)
echo "round11-heavy: ALL PASS"
