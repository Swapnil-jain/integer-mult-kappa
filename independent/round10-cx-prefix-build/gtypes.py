"""Classify #144's gauge annihilators A (= sigma^perp, dim u = h - d) for nullity-dependent sharing:
nondegenerate? radical dim r, alternating (all-even, i.e. 1 in sigma)?"""
import pickle, sys
from collections import Counter
sys.path.insert(0, '.')
from cube_word import basis, perp
def rank2(vs):
    piv = {}
    for x in vs:
        while x:
            b = x.bit_length() - 1
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)
def gram_rank(B):
    return rank2([sum((bin(a & b).count('1') & 1) << j for j, b in enumerate(B)) for a in B])
def classify(A, h):
    A = basis(A); u = len(A); r = u - gram_rank(A)
    alt = all(bin(x).count('1') % 2 == 0 for x in A)
    return (u, r, alt)
if __name__ == '__main__':
    W = pickle.load(open('word144.pkl', 'rb')); h = 24
    T = Counter(classify(z['A'], h) for z in W['sel'])
    print('gauge types (u, radical r, alternating):', dict(T))
    print('first-frame dims of gauged roles:', dict(Counter(z['r'] for z in W['sel'])))
    print('targets per gauge:', dict(Counter(len(z['targets']) for z in W['sel'])))
