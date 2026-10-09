"""Round-eleven witness: round ten's certificate (scripts/certificate_round10.py), unchanged, on the round-eleven
inventories. Both sides are measured histograms (no new layer); the assembly, the stopped mix, the finite bridge and
the 47 rows are round ten's code. This file only adds the round-ten regression and writes to round eleven's folders.

Regressions, before ours: PR #144's published kappa and round nine's (through round ten's own regressions), then
round ten's frozen kappa 6153378/10^10 recomputed from certificates/round10/.

Usage:
  python3 scripts/certificate_round11.py                       # regressions + the frozen round-eleven files
  python3 scripts/certificate_round11.py BIT.json CX.json      # evaluate other inventories (nothing written)
  python3 scripts/certificate_round11.py --freeze BIT.json CX.json
        # pin the computed savings as claims; write certificates/round11/{bit,cx}_inventory.json.gz, summary.json
        # and lean/round11-histograms.json"""
import json, os, sys
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certificate_round10 as c10
from certificate_round10 import load

FROZEN = os.path.join(c10.ROOT, 'certificates', 'round11')
FROZEN10 = c10.FROZEN
LEAN_NAME = 'round11-histograms.json'
R10_KAPPA = Q(6153378, 10**10)


def regressions(verbose=True):
    out = c10.regressions(verbose)
    r10 = c10.certify(load(os.path.join(FROZEN10, 'bit_inventory.json.gz')),
                      load(os.path.join(FROZEN10, 'cx_inventory.json.gz')))['k']['kappa']
    with open(os.path.join(FROZEN10, 'summary.json')) as f:
        pinned = Q(json.load(f)['kappa'])
    out['round 10 kappa, frozen'] = r10 == R10_KAPPA == pinned
    if verbose:
        print('  %-50s %s' % ('round 10 kappa %s, frozen' % R10_KAPPA, 'PASS' if out['round 10 kappa, frozen'] else 'FAIL'))
    assert all(out.values()), [n for n, ok in out.items() if not ok]
    return out


def main(argv):
    do_freeze = '--freeze' in argv
    args = [a for a in argv if not a.startswith('--')]
    c10.c9.reproduce_pr144()
    regressions()
    if args:
        bpath, cpath = args
    else:
        bpath, cpath = os.path.join(FROZEN, 'bit_inventory.json.gz'), os.path.join(FROZEN, 'cx_inventory.json.gz')
    bit, cx = load(bpath), load(cpath)
    assert bit['side'] == 'bit' and cx['side'] == 'complex'
    # round eleven uses the ordinary stop only; the ladder stays a round-ten research switch
    assert int(bit['cert'].get('ladder', 1)) == 1, 'round eleven freezes no ladder'
    r = c10.certify(bit, cx)
    c10.report(r, bit, cx)
    assert r['k']['ok'], r['k']['bad']
    if do_freeze:
        c10.freeze(r, bit, cx, frozen=FROZEN, lean_name=LEAN_NAME)
    return r


if __name__ == '__main__':
    main(sys.argv[1:])
