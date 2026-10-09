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
L1 = 'data/local_L1.json'; TB31 = 'data/tmod_TB31_L1f8_5.9717259e-4.json'
H56 = 'data/pmod_H56snap_w02_5.6251423e-4.json'; QB01 = 'data/qmod_anneal_best01.json'
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
    # own annealed all-but-one modules (pmod_anneal3.py KIND=allbut, started from the nested prefix)
    'qa_G37_rx': dict(merge=0, tri=TRI11, pairmod=G37, allbut='data/qmod_anneal_a.json', relaxed=True),
    'qa_G37_h20b_rx': dict(merge=0, trisrc=('data/h20_g1.json.gz', [0, 1, 2, 5, 7, 8, 11, 12, 14, 16, 19]), pairmod=G37, allbut='data/qmod_anneal_a.json', relaxed=True),
    'pf_small_rx': dict(merge=0, allbut='prefix', relaxed=True),   # small-p analogue (default PR117 restrictions)
    # #168 v3 (commit fd25adb7, eumemic, Apache-2.0, data): local_L1 cube-local circuit, TB31 triple, H56snap pair,
    # annealed all-but-one, 'f8:00111100' fused reads; our matching / descent / gauges / reuse on top
    'v3_rx': dict(local=L1, trimod=TB31, pairmod=H56, allbut=QB01, fuse='00111100', relaxed=True),
    # #168 v4 (commit 4a3c769, data): J0 pair, TE_TD_TB3 triple, climb3u all-but-one; local_L1 and fuse unchanged
    'v4_rx': dict(local=L1, trimod='data/tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json', pairmod='data/pmod_J0_full_6.0666810e-4.json',
                  allbut='data/qmod_climb3u_best.json', fuse='00111100', relaxed=True),
    'v4f': dict(local=L1, trimod='data/tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json', pairmod='data/pmod_J0_full_6.0666810e-4.json',
                allbut='data/qmod_climb3u_best.json', fuse='00111100', arcs='data/arcs_v4.json', nodesc=True),
    'v4fib_rx': dict(local=L1, trimod='data/tmod_TE_TD_TB3_1_1_4_full_6.0617964e-4.json', pairmod='data/pmod_J0_full_6.0666810e-4.json',
                     allbut='data/abo_fib_9.json', fuse='00111100', relaxed=True),   # defrag's cyclic Fibonacci all-but-one (ours)
    'v3small_rx': dict(local=L1, fuse='00111100', allbut='prefix', relaxed=True),     # small-p analogue
    'e168tree_rx': dict(merge=0, trisrc=H20, pairmod=G37, relaxed=True),
}

def build_cfg(p, cfg):
    o = (cw.mod_pairs, cw.mod_triples, cw.MERGE, cw.mod_allbut, cw.LOCAL, cw.FUSE, cw.LOCAL_BUG, cw.LOCALFN)
    if cfg.get('gen'):   # localmotif: generalised recipe builder
        sys.path.insert(0, os.path.dirname(HERE)); import lmlocal; cw.LOCALFN = lmlocal.gen_local
    if 'local' in cfg: cw.LOCAL = json.load(open(os.path.join(HERE, cfg['local']))) if isinstance(cfg['local'], str) else cfg['local']
    if 'fuse' in cfg: cw.FUSE = cfg['fuse']
    if 'local_bug' in cfg: cw.LOCAL_BUG = cfg['local_bug']
    if 'trimod' in cfg: TM = json.load(open(os.path.join(HERE, cfg['trimod']))); cw.mod_triples = lambda pp, w: TM
    if 'trisrc' in cfg: T = tri_from(*cfg['trisrc']); cw.mod_triples = lambda pp, w: T
    if cfg.get('allbut') == 'prefix': cw.mod_allbut = allbut_prefix
    if isinstance(cfg.get('allbut'), str) and cfg['allbut'].endswith('.json'):
        QF = json.load(open(os.path.join(HERE, cfg['allbut']))); QF = QF.get('module', QF); cw.mod_allbut = lambda n: QF
    if isinstance(cfg.get('allbut'), dict): QM = cfg['allbut']; cw.mod_allbut = lambda n: QM
    if 'pair' in cfg: st = cfg['pair']; cw.mod_pairs = lambda n, w: pair(st['L'], st['s'], st['r'])
    if 'pairmod' in cfg:   # a direct pair-disjoint module (our module format), loaded as data
        M = json.load(open(os.path.join(HERE, cfg['pairmod']))) if isinstance(cfg['pairmod'], str) else cfg['pairmod']
        M = M.get('module', M)
        cw.mod_pairs = lambda n, w: M
    if 'tri' in cfg: L = cfg['tri']; cw.mod_triples = lambda pp, w: tri(L)
    cw.MERGE = cfg.get('merge')
    try: g = cw.build(p)
    finally: cw.mod_pairs, cw.mod_triples, cw.MERGE, cw.mod_allbut, cw.LOCAL, cw.FUSE, cw.LOCAL_BUG, cw.LOCALFN = o
    return g

def make(p, name, sweeps=6, cache=True):
    if name.startswith('@'):   # localmotif: '@path.json' = {cfg} with a recipe
        CFG[name] = json.load(open(name[1:]))
    fn = os.path.join(HERE, 'word_%s_p%d.pkl' % (name, p))
    g = build_cfg(p, CFG[name]); h = 2 * p
    if CFG[name].get('arcs'):
        # frozen upstream carrier arcs (data) compiled in compile_closure's order (frame rank, plain-heap topological
        # index): with #168's graph (same sha) this reproduces its ops, roles and gauges index for index, so its
        # frozen per-op frames and pairs apply to our word unchanged (recount re-checks every nesting)
        cw.RELAXED = True; o_ = cw.ORDER_RANK; cw.ORDER_RANK = True
        try: prof, wit = cw.compile_frames(g, json.load(open(os.path.join(HERE, CFG[name]['arcs']))))
        finally: cw.RELAXED = False; cw.ORDER_RANK = o_
    elif CFG[name].get('relaxed'):
        from relaxed import relaxed_arcs
        cw.RELAXED = True
        try: prof, wit = cw.compile_frames(g, relaxed_arcs)
        finally: cw.RELAXED = False
    else: prof, wit = cw.compile_frames(g, lambda D: maxcard(D)[0])
    if os.environ.get('NODESC') or CFG[name].get('nodesc'):   # caps frames: leave node descent to pdescent
        pr, s = cw.select_gauges(prof, wit)
        return dict(p=p, h=h, m=3 * h, g=g, prof=prof, wit=wit, pr=pr, s=s, best=None, prof0=prof, wit0=wit)
    if cache and os.path.exists(fn):
        ann, rank, H = pickle.load(open(fn, 'rb'))
    else:
        ann, rank, H = descend(prof, wit, sweeps=sweeps, verbose=False); pickle.dump((ann, rank, H), open(fn, 'wb'))
    prof2, wit2 = dapply(prof, wit, ann, rank, H)
    pr, s = cw.select_gauges(prof2, wit2)
    return dict(p=p, h=h, m=3 * h, g=g, prof=prof2, wit=wit2, pr=pr, s=s, best=None, prof0=prof, wit0=wit)
