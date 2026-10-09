"""Build a cube word variant (smallh-search levers: searched pair/triple restrictions, merged hyperplane reads),
compile with our max-card matching, run cube-stack's frame descent, select gauges. Returns the cbirth word dict."""
import sys, os, json, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cube_word as cw
from matching import maxcard
from descent import run as descend, apply as dapply
from restrict import pair, tri, tri_from, allbut_prefix

TRI11 = [1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23]
G37 = 'data/pmod_G37_w02_5.6098194e-4.json'
H20 = ('data/h20_g1.json.gz', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 16])
CFG = {
    'default': {},
    'merge': dict(merge=0),
    'p12best': dict(merge=0, pair={"L": [2, 4, 5, 6, 7, 8, 9, 10, 11, 15, 17], "s": 18, "r": 22},
                    tri=[1, 2, 5, 6, 8, 10, 12, 14, 17, 18, 21, 23]),
    'tri11': dict(merge=0, tri=[1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23]),     # tsearch.py p=11 seed 1 (caps 5.2301e-4)
    'tri13': dict(merge=0, tri=[1, 2, 5, 6, 8, 9, 11, 13, 15, 17, 19, 21, 22]),  # tsearch.py p=13 seed 4 (caps 5.2348e-4)
    # eumemic's annealed pair-disjoint module (fork branch paired-cube-merge-annealed, read as DATA, Apache-2.0)
    'pmod11e': dict(merge=0, tri=[1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23], pairmod='data/pmod_C35_eumemic.json'),
    'pmod11e_dt': dict(merge=0, pairmod='data/pmod_C35_eumemic.json'),
    'merge_rx': dict(merge=0, relaxed=True),
    'tri11_rx': dict(merge=0, tri=[1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23], relaxed=True),
    'tri13_rx': dict(merge=0, tri=[1, 2, 5, 6, 8, 9, 11, 13, 15, 17, 19, 21, 22], relaxed=True),
    'pmod11e_rx': dict(merge=0, tri=[1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23], pairmod='data/pmod_C35_eumemic.json', relaxed=True),
    'p12best_rx': dict(merge=0, pair={"L": [2, 4, 5, 6, 7, 8, 9, 10, 11, 15, 17], "s": 18, "r": 22},
                       tri=[1, 2, 5, 6, 8, 10, 12, 14, 17, 18, 21, 23], relaxed=True),
    'p12pair': dict(merge=0, pair={"L": [2, 4, 5, 6, 7, 8, 9, 10, 11, 15, 17], "s": 18, "r": 22}),
    # cube-prefix: #168's nested-prefix all-but-one module (own generator, equal to their qmod JSON) and, as DATA
    # (eumemic, Apache-2.0, sha-pinned in #168's SOURCE.json), their h20_g1 triple witness and G37 pair module
    'pf_pmod11e_rx': dict(merge=0, tri=TRI11, pairmod='data/pmod_C35_eumemic.json', allbut='prefix', relaxed=True),
    'pf_G37_rx': dict(merge=0, tri=TRI11, pairmod=G37, allbut='prefix', relaxed=True),
    'e168_rx': dict(merge=0, trisrc=H20, pairmod=G37, allbut='prefix', relaxed=True),
    'pf_G37_dt_rx': dict(merge=0, pairmod=G37, allbut='prefix', relaxed=True),
    'e168tree_rx': dict(merge=0, trisrc=H20, pairmod=G37, relaxed=True),
}

def build_cfg(p, cfg):
    o = (cw.mod_pairs, cw.mod_triples, cw.MERGE, cw.mod_allbut)
    if 'trisrc' in cfg: T = tri_from(*cfg['trisrc']); cw.mod_triples = lambda pp, w: T
    if cfg.get('allbut') == 'prefix': cw.mod_allbut = allbut_prefix
    if 'pair' in cfg: st = cfg['pair']; cw.mod_pairs = lambda n, w: pair(st['L'], st['s'], st['r'])
    if 'pairmod' in cfg:   # a direct pair-disjoint module (our module format), loaded as data
        M = json.load(open(os.path.join(HERE, cfg['pairmod']))) if isinstance(cfg['pairmod'], str) else cfg['pairmod']
        cw.mod_pairs = lambda n, w: M
    if 'tri' in cfg: L = cfg['tri']; cw.mod_triples = lambda pp, w: tri(L)
    cw.MERGE = cfg.get('merge')
    try: g = cw.build(p)
    finally: cw.mod_pairs, cw.mod_triples, cw.MERGE, cw.mod_allbut = o
    return g

def make(p, name, sweeps=6, cache=True):
    fn = os.path.join(HERE, 'word_%s_p%d.pkl' % (name, p))
    g = build_cfg(p, CFG[name]); h = 2 * p
    if CFG[name].get('relaxed'):
        from relaxed import relaxed_arcs
        cw.RELAXED = True
        try: prof, wit = cw.compile_frames(g, relaxed_arcs)
        finally: cw.RELAXED = False
    else: prof, wit = cw.compile_frames(g, lambda D: maxcard(D)[0])
    if cache and os.path.exists(fn):
        ann, rank, H = pickle.load(open(fn, 'rb'))
    else:
        ann, rank, H = descend(prof, wit, sweeps=sweeps, verbose=False); pickle.dump((ann, rank, H), open(fn, 'wb'))
    prof2, wit2 = dapply(prof, wit, ann, rank, H)
    pr, s = cw.select_gauges(prof2, wit2)
    return dict(p=p, h=h, m=3 * h, g=g, prof=prof2, wit=wit2, pr=pr, s=s, best=None, prof0=prof, wit0=wit)
