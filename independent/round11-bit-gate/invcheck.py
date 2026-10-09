"""Compare a ledger measurement (ledger.py --out) with the release inventory (certificate_round10 schema): the
inventory's child histogram and W must be exactly the histogram the replay measured, and its claimed a* the one the
replay certified. Standard library only.  Usage: python3 invcheck.py LEDGER_OUT.json INVENTORY.json[.gz]"""
import sys, json, gzip
from fractions import Fraction as Q
led = json.load(open(sys.argv[1])); p = sys.argv[2]
inv = json.load(gzip.open(p) if p.endswith('.gz') else open(p))
k = int(inv.get('scale', 1))
H = {int(r): Q(n) / k for r, n in inv['child_histogram'].items() if Q(n)}
L = {int(r): Q(n) for r, n in led['hist'].items() if Q(n)}
ok = led['all'] and H == L and Q(inv['W_per_vertex']) / k == Q(led['W']) and Q(inv['cert']['astar']) == Q(led['a']) \
     and inv['side'] == 'bit' and 3 * int(inv['h']) == int(inv['m'])
print('RESULT invcheck: inventory == replay measurement (%d widths, W %s, a* %s): %s' % (len(H), led['W'], led['a'], 'ALL PASS' if ok else 'FAIL'))
sys.exit(0 if ok else 1)
