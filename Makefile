.PHONY: verify roles round4 round7 round8 round8-heavy round9 round9-heavy round10 round10-heavy round11 round11-heavy lean notes
verify:
	python3 scripts/check_identities.py
	python3 scripts/check_segmented_inverse.py
	python3 scripts/check_subblock_inverse.py
	python3 scripts/check_reuse_inverse.py
	python3 scripts/check_reuse_inverse.py 113 128 8 120 50
	! python3 scripts/check_reuse_inverse.py 241 256 17 150 60 1
	python3 scripts/check_gs_path.py
	python3 scripts/check_fine_field.py
	python3 scripts/certificate.py
	python3 -m unittest discover -s tests -v

roles:
	cd independent/pr7-role-counts && python3 export.py 26 /tmp/local26.txt && clang++ -O2 -std=c++17 global.cpp -o /tmp/global && /tmp/global 28 /tmp/local26.txt && python3 complexside.py 28

round4:
	python3 independent/complex-network/fullbatch_hist.py 28 /tmp/fb28sf.json sf
	python3 scripts/certificate_round4.py /tmp/fb28sf.json
	python3 independent/two-stage-bit/checkside.py 47

round7:
	python3 scripts/certificate_round7.py
	python3 independent/joint-frame-stack/check_frames.py
	python3 independent/joint-frame-stack/check_schedule.py
	python3 independent/joint-frame-stack/check_design.py
	python3 independent/deferred-readout/check_word.py 23
	python3 independent/deferred-readout/check_frames.py 23
	python3 independent/deferred-readout/check_lifted.py 23
	python3 independent/deferred-readout/check_stair.py 23
	python3 independent/deferred-readout/check_stair.py 23 --control | grep -q 'certified=False'

# Round eight, fast checks (standard library, minutes on a laptop): the certificate from the frozen histograms,
# the bit histogram rebuilt from the frozen round-seven schedule, the outer wrapper and remainder pairing, and the
# PR #128 partition.
round8:
	python3 scripts/certificate_round8.py
	python3 independent/round8-oppbank/bit_hist.py
	python3 independent/round8-oppbank/wrap.py | grep 'PASS$$'
	python3 -I independent/round8-coreshare/partition_check.py certificates/round8/pr128_partition_mrp24.json | grep 'PASS$$'
	python3 -m unittest tests.test_round8 -v

# Round eight, heavy checks (a many-core machine with numpy and a C compiler): the complex histogram rebuilt from the
# word (about 7 GB), the independent exact Q(i) replay of the shared word in /tmp/round8-cx (build about 9 GB,
# replay about 9 GB, two seeds), its walk, exterior and certificate checks, and the all-edges bit gate (J processes
# of about 0.75 GB each).
round8-heavy:
	python3 independent/round8-coreshare/complex_hist.py
	mkdir -p /tmp/round8-cx
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/build_word.py 24 allE,alt,links,dag117 pr128 w24.pkl
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/sreplay.py w24.pkl 1 own | grep "'PASS': True"
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/sreplay.py w24.pkl 2 own | grep "'PASS': True"
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/walk2.py w24.pkl | grep "'PASS': True"
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/excheck.py w24.pkl | grep '"PASS": true'
	cd /tmp/round8-cx && python3 $(CURDIR)/independent/round8-complex-gate/cert2.py hist_w24_1_own.pkl
	python3 independent/round8-coreshare/complex_hist.py --gate /tmp/round8-cx/hist_w24_1_own.pkl
	cd independent/round8-complex-gate && python3 endpoint.py
	independent/round8-oppbank/run_gate.sh

# Round nine, fast checks (standard library, seconds): PR #144's published kappa, constraints and margins reproduced
# from its published inventories, then our certificate from the frozen inventories, and the tests.
round9:
	python3 scripts/certificate_round9.py
	python3 -m unittest tests.test_round9 -v

# Round nine, heavier checks (numpy; about 1.5 GB and a few minutes each): the source-bound ledger of witness 2 with
# gauge omission and its negative controls, the independent re-certification, the role-by-role recount of PR #144's
# complex module, the cube algebra, literal replays and cross-stage sharing, and the three-stage cover replay.
round9-heavy:
	cd independent/round9-bit-ledger && python3 ledger.py --out /tmp/round9_ledger.json | grep 'ALL PASS'
	cd independent/round9-bit-ledger && python3 cert.py /tmp/round9_ledger.json | grep '"kappa_floor_1e10": 4663738'
	for c in drop_prelude_read wrong_coefficient chain_not_subsequence late_prelude; do \
	  ! (cd independent/round9-bit-ledger && python3 ledger.py --ctl $$c 2>&1 | grep -q 'ALL PASS') || exit 1; done
	cd independent/round9-complex-recount && python3 recount144.py | grep 'children == published: True'
	cd independent/round9-audit144 && python3 cube_alg.py 4 5 | grep 'ALL OK'
	cd independent/round9-audit144 && python3 cube_lit.py 5 lit 6 gauge omit | grep 'ALL PASS'
	cd independent/round9-audit144 && ! python3 cube_lit.py 5 lit 7 gauge omit control=noprelude | grep -q 'ALL PASS'
	cd independent/round9-audit144 && python3 share_lit.py 3 1 | grep 'ALL PASS'
	cd independent/round9-audit144 && ! python3 share_lit.py 3 1 lockstep | grep -q 'ALL PASS'
	cd independent/round9-cover-e2e && python3 cover_e2e.py 10 allE,dual+,alt,links,pasm,ivec,lift,defer,vleaf,clos 1 2 | grep "'PASS': True"

# Round ten, fast checks (standard library, seconds): PR #144's and round nine's kappa reproduced by the round-ten
# certificate (parts and bare histograms), then our certificate from the frozen inventories, and the tests.
round10:
	python3 scripts/certificate_round10.py
	python3 -m unittest tests.test_round10 -v

# Round ten, heavier checks (numpy; see the script for sizes): the gates' replays, recounts and certificate checks
# of both sides, each with its negative controls.
round10-heavy:
	bash independent/round10-heavy.sh

# Round eleven, fast checks (standard library, seconds): round ten's certificate on the round-eleven inventories, with
# PR #144's, round nine's and round ten's kappa reproduced first, and the tests.
round11:
	python3 scripts/certificate_round11.py
	python3 -m unittest tests.test_round11 -v

# Round eleven, heavier checks (numpy; about 3 GB and 10 minutes): the JSON-only bit ledger at p = 12, 9 and 7 with
# its mutation controls, and the recycling checker of the complex word at p = 11, 9 and 7 with its controls.
round11-heavy:
	bash independent/round11-heavy.sh

lean:
	python3 lean/gen.py
	lean lean/Round5.lean
	python3 lean/gen.py round6-histograms.json Round6.lean
	lean lean/Round6.lean
	python3 lean/dump7.py
	python3 lean/gen.py round7-histograms.json Round7.lean
	lean lean/Round7.lean
	python3 lean/gen.py round8-histograms.json Round8.lean
	lean lean/Round8.lean
	python3 lean/gen.py round9-histograms.json Round9.lean
	lean lean/Round9.lean
	python3 lean/gen.py round10-histograms.json Round10.lean
	lean lean/Round10.lean
	python3 lean/gen.py round11-histograms.json Round11.lean
	lean lean/Round11.lean

notes:
	tectonic -X compile --outdir artifacts notes/segmented-inverse.tex
	tectonic -X compile --outdir artifacts notes/stack-notes.tex
	tectonic -X compile --outdir artifacts notes/fine-field-lemma.tex
	tectonic -X compile --outdir artifacts notes/partial-swap-batching.tex
	tectonic -X compile --outdir artifacts notes/round3-combination.tex
	tectonic -X compile --outdir artifacts notes/data-edge-batching.tex
	tectonic -X compile --outdir artifacts notes/complex-source-frames.tex
	tectonic -X compile --outdir artifacts notes/deferred-readout.tex
