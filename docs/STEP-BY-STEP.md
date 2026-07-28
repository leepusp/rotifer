1) Install sphinx and theme:
pip install sphinx sphinx-book-theme

2) Generate autodoc

sphinx-apidoc -f -o docs/api PATH/TO/DESIRED/PACKAGE setup.py
ex: lib/rotifer documents most of it.

3) Build the documentation:
cd docs
make html

Run it as standard html (f5 on VS CODE).

This is work on progress. Do not rely on this version of docs for research work under any circumstance.
