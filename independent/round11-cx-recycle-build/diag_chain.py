import os, pickle, sys
from collections import Counter
os.environ['NODESC'] = '1'
from resel_br import *
import make_word, pdescent
W = make_word.make(11, sys.argv[1]); PL = pickle.load(open(sys.argv[2], 'rb'))
W2 = with_sel(W, PL['sel']); F = dict(gfirst_frames(W2), frames=PL['frames'])
seqs = pdescent.build_seqs(W2, F, PL['pairs']); ops = F['ops']
pat = Counter(); ex = []
for q in seqs:
    d = [len(it[1]) if it[0] == 'F' else len(F['frames'][it[1]]) for it in q]
    steps = tuple(b - a for a, b in zip(d, d[1:]) if b > a)
    pat[steps] += 1
    if steps[-2:] == (4, 15) and len(ex) < 5: ex.append((d, [it if it[0] == 'O' else ('F', len(it[1])) for it in q]))
for s, c in pat.most_common(25): print(c, s)
for d, q in ex:
    print('dims', d)
    for it in q:
        if it[0] == 'O':
            i = it[1]; a, b, x = ops[i]; print('   op', i, 'roles', a, b, 'node', x, 'dim', len(F['frames'][i]), 'span', len(spans_of(W2['wit']['g'])[x]))
