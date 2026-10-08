.PHONY: verify
verify:
	python3 scripts/check_identities.py
	python3 scripts/check_segmented_inverse.py
	python3 scripts/certificate.py
	python3 -m unittest discover -s tests -v
