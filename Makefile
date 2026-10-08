.PHONY: verify roles notes
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

notes:
	tectonic -X compile --outdir artifacts notes/segmented-inverse.tex
	tectonic -X compile --outdir artifacts notes/stack-notes.tex
	tectonic -X compile --outdir artifacts notes/fine-field-lemma.tex
	tectonic -X compile --outdir artifacts notes/partial-swap-batching.tex
	tectonic -X compile --outdir artifacts notes/round3-combination.tex
