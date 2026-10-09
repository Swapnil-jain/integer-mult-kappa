"""Exhaustive mixed-gap pattern scan (#173's idea on our complex word): worker k of K scores patterns k, k+K, ...
with score_cp (tri11 + G37 + merged reads + mixed-gap all-but-one, caps). Usage: python3 qmix.py k K"""
import sys, json, itertools
from score_cp import score
from restrict import allbut_mixed
from make_word import TRI11, G37
k, K = int(sys.argv[1]), int(sys.argv[2])
pats = [''.join(t) for t in itertools.product('LR', repeat=7)]
for pat in pats[k::K]:
    s = score(11, dict(merge=0, tri=TRI11, pairmod=G37, allbut=allbut_mixed(9, pat)))
    print(json.dumps(dict(pat=pat, score=s)), flush=True)
