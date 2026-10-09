#!/bin/bash
# All-edges gate of the round-eight bit side (opposite bank orders on round-seven witness 2). Heavy: run it on a
# many-core machine (about 0.75 GB per process; J parallel processes, default 12).
#   1. export.py re-runs check_design.py's exact head on the frozen round-seven data, collects every charged residual
#      step and aux sigma, rebuilds the one-run histogram and writes the edge blob (about 100 MB) to $GATE_OUT;
#   2. C1 on every edge at ONE 60-bit random integer basis X: fast passes mod 2^29-3 (small ranks, copies) and
#      2^21-9 (aux edges, data entrances); every zero candidate is re-checked mod 2^61-1, where a nonzero minor is
#      nonzero over Q. A small basis gives genuine zero minors, so the basis must stay 60-bit;
#   3. dense cross-check, Bruhat sample (one reversed contiguous run), negative controls (neg, negI);
#   4. wrap.py (outer wrapper, remainder pairing) and acct.py (histogram == frozen file, moment, stopped mix, kappa).
# Usage: GATE_OUT=/path J=12 independent/round8-oppbank/run_gate.sh
set -e
D="$(cd "$(dirname "$0")" && pwd)"; R="$(cd "$D/../.." && pwd)"
OUT=${GATE_OUT:-/tmp/round8-bitgate}; J=${J:-12}; mkdir -p "$OUT"
[ -s "$OUT/gate_data.bin" ] || (cd "$R" && GATE_OUT="$OUT" python3 independent/round8-oppbank/export.py)
cd "$OUT"
for p in gate29 gated gate_p61; do gcc -O3 -march=native -o $p "$D/$p.c" -lm; done
chunks() {  # prog mode total nchunks tag : run [lo, hi) chunks, J at a time, skipping finished ones
  local prog=$1 mode=$2 n=$3 k=$4 tag=$5 i lo hi
  for i in $(seq 0 $((k - 1))); do
    lo=$((i * n / k)); hi=$(((i + 1) * n / k))
    grep -qx "${tag}_$i done" progress.txt 2>/dev/null && continue
    while [ "$(jobs -rp | wc -l)" -ge "$J" ]; do sleep 5; done
    (./$prog gate_data.bin $mode $lo $hi > cand_${tag}_$i.txt 2> log_${tag}_$i.txt && echo "${tag}_$i done" >> progress.txt) &
  done
  wait
}
chunks gate29 small 122504 24 small
chunks gate29 copy 1 1 copy
chunks gated aux 4231 24 aux
chunks gated entr 1771 24 entr
cat cand_small_*.txt cand_copy_*.txt cand_aux_*.txt cand_entr_*.txt | grep -v RESULT > cand_all.txt || true
./gate_p61 gate_data.bin recheck cand_all.txt > recheck_all.txt 2> log_recheck.txt
./gate_p61 gate_data.bin dense 4 > /dev/null 2> log_dense.txt
./gate_p61 gate_data.bin bruhat 2 > bruhat.txt 2> log_bruhat.txt
./gate_p61 gate_data.bin neg > neg.txt 2> log_neg.txt
./gate29 gate_data.bin small 0 3000 negI > cand_negI.txt 2> log_negI.txt
grep -h DONE log_small_*.txt log_copy_*.txt log_aux_*.txt log_entr_*.txt | sed 's/^\[[^]]*\] //'
tail -3 recheck_all.txt; grep -h DONE log_dense.txt; tail -1 log_bruhat.txt; cat neg.txt; tail -1 cand_negI.txt
cd "$R" && python3 independent/round8-oppbank/wrap.py && GATE_OUT="$OUT" python3 independent/round8-oppbank/acct.py
