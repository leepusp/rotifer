# Contributing to the documentation

This section is for anyone who writes ROTIFER code or improves this
site: laboratory members documenting their own modules, and outside
contributors sending patches. It covers three things.

- {doc}`documenting-code` explains how to write docstrings the way
  this project renders them: the NumPy standard, parsed by numpydoc
  and nothing else.
- {doc}`building-the-docs` explains how to build the site locally,
  what each Makefile target does, and how to fix the build failures
  this repository is known to produce.
- {doc}`design-system` explains where the visual design lives and
  how to extend the site without breaking its consistency.

The {doc}`data access section </db/index>` and the generated pages
of {py:mod}`rotifer.db` are the reference implementation: when in
doubt about how a docstring or a section page should look, copy
what `rotifer.db` does.

```{toctree}
:hidden:

documenting-code
building-the-docs
design-system
```
