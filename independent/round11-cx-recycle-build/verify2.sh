#!/bin/bash
# verify2.sh CFG.json P ORDER : full stack on a fresh compile of a recipe graph.
# final2 (NODESC=1 PD_IT=4: pdescent from caps, rematch) with dirty replay + controls -> export -> check_word ->
# t176 terminal-sink selection with its own replay/controls -> add_sinks -> check_word2 (independent, sinks included).
# ORDER = fwd|reverse is the greedy arc-extension scan; it must be the same in every step (export rebuilds the word).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NODESC=1 PD_IT=4 RX_ORDER=$3
c=$1; p=$2; t=$(basename $c .json)_$3; export TAG=_$3
mkdir -p rel
nice -n 3 python3 final2.py $p @$c replay > out_final2_${t}_p$p.txt 2>&1
nice -n 3 python3 export_word.py $p @$c plan2_@${c}_p$p$TAG.pkl rel/${t}_p$p.json.gz > out_export_${t}_p$p.txt 2>&1
nice -n 3 python3 check_word.py rel/${t}_p$p.json.gz > out_check_${t}_p$p.txt 2>&1
nice -n 3 python3 -I t176.py rel/${t}_p$p.json.gz --out t176_${t}_p$p.json > out_t176_${t}_p$p.txt 2>&1
nice -n 3 python3 add_sinks.py rel/${t}_p$p.json.gz t176_${t}_p$p.json rel/${t}_sinks_p$p.json.gz > out_sinks_${t}_p$p.txt 2>&1
nice -n 3 python3 check_word2.py rel/${t}_sinks_p$p.json.gz > out_check2_${t}_p$p.txt 2>&1
echo done $t $p
