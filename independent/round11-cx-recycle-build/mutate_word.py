"""Negative controls for check_word.py: write mutated copies of a frozen word; each must make the checker FAIL.
  sign    one signed DAG addition flipped (decoder identity and the replay must break)
  frame   one donor's last op frame replaced by the full space (its chain then cannot nest into the recipient gauge)
  claim   the claimed a_c moved to the next 1e-12 grid point (the certificate must reject it)
  pair    one recipient re-paired to a donor whose last op is after the recipient's read clock (chronology)
Usage: python3 mutate_word.py WORD.json.gz OUTDIR"""
import sys, gzip, json, os
from fractions import Fraction as Q
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
J0 = json.load(gzip.open(src, 'rt'))
def save(name, J):
    with gzip.open(os.path.join(out, 'mut_%s.json.gz' % name), 'wt') as f: json.dump(J, f, separators=(',', ':'))
J = json.loads(json.dumps(J0)); x = next(k for k, s in enumerate(J['signs']) if s == -1); J['signs'][x] = 1; save('sign', J)
J = json.loads(json.dumps(J0))
pos = {i: k for k, i in enumerate(J['sched'])}; last = {}
for i in J['sched']:
    a, b, _ = J['ops'][i]; last[a] = i; last[b] = i
d0 = J['pairs'][0][0]; J['frames'][last[d0]] = [1 << j for j in range(J['h'] - 1, -1, -1)]; save('frame', J)
J = json.loads(json.dumps(J0)); J['claim']['a_c'] = str(Q(J['claim']['a_c']) + Q(1, 10 ** 12)); save('claim', J)
J = json.loads(json.dumps(J0)); taken = {a for a, _ in J['pairs']} | {b for _, b in J['pairs']}
roots = set(J['rootroles']); b = J['pairs'][0][1]; tb = J['tau'][str(b)]
late = next(a for a in last if a not in taken and a not in roots and pos[last[a]] >= tb)
J['pairs'][0][0] = late; save('pair', J)
print('wrote 4 mutated words to', out)
