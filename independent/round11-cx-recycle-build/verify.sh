#!/bin/bash
# verify.sh CFG.json P : final2 (descent + rematch) with dirty-scratch replay and controls, export, independent check
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
c=$1; p=$2; t=$(basename $c .json)
mkdir -p rel
nice -n 3 python3 final2.py $p @$c replay > out_final2_${t}_p$p.txt 2>&1
nice -n 3 python3 export_word.py $p @$c plan2_@${c}_p$p.pkl rel/${t}_p$p.json.gz > out_export_${t}_p$p.txt 2>&1
nice -n 3 python3 check_word.py rel/${t}_p$p.json.gz > out_check_${t}_p$p.txt 2>&1
echo done $t $p
