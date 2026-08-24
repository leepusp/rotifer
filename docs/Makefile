# Makefile for the ROTIFER documentation.
#
# Targets:
#   html      build the HTML site into _build/html
#   livehtml  serve the site with live reload (sphinx-autobuild)
#   strict    build treating every warning as an error (what CI runs)
#   doctest   run the doctest examples embedded in docstrings and pages
#   linkcheck check external links
#   clean     remove all build artifacts

SPHINXBUILD   ?= sphinx-build
SPHINXAUTO    ?= sphinx-autobuild
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: html livehtml strict doctest linkcheck clean

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

livehtml:
	$(SPHINXAUTO) "$(SOURCEDIR)" "$(BUILDDIR)/html" --open-browser

strict:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" -W --keep-going -n

doctest:
	$(SPHINXBUILD) -b doctest "$(SOURCEDIR)" "$(BUILDDIR)/doctest"

linkcheck:
	$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck"

clean:
	rm -rf "$(BUILDDIR)"
