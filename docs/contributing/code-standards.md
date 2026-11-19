# Code standards

For anyone who writes code for ROTIFER, lab member or outside contributor. MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

## Development environment

Follow the installation instructions in the {ghsrc}`README.md` to create the `rotifer` conda environment and put the libraries and tools on your `PYTHONPATH` and `PATH`. Working from a clone with `conda-develop` means your edits take effect immediately, without reinstalling.

## Where code goes

Library code goes in `lib/rotifer/<package>`, in the package whose topic it belongs to. Command line tools go in `bin/`. Python scripts in `bin/` are normally extensionless, and that is fine.

Nothing new goes under `lib/rotifer/devel/`. Exploratory work and one-off analyses belong in Jupyter notebooks outside the rotifer repository, importing rotifer like any other user. See {doc}`personal-codebases`.

## Reuse before writing

ROTIFER is large and spread across several usage cases, so functionality is not always where you would expect it. Before writing something new you MUST search the repository for existing functionality:

```bash
git grep -i <keyword>
```

Try a few synonyms and also search `bin/`, whose scripts have no `.py` suffix. If something similar exists but does not quite fit, extend it. Extending existing code is a welcome kind of contribution. A pull request that duplicates existing functionality will be asked to extend the existing code instead.

If you are unsure, open an issue and ask before investing time in an implementation, especially for larger changes. Discussing the design first is much less frustrating than having a finished pull request turned down.

## Documentation

Public functions, classes and command line tools MUST have NumPy-style docstrings. The conventions this project renders are described in {doc}`documentation/documenting-code`.

A change to public API MUST update the documentation in the same pull request. Documentation that lands later tends not to land.

## Tests

New behaviour MUST come with tests. Tests live under `test/` and mirror the package path: tests for `lib/rotifer/db/ncbi/entrez.py` go under `test/db/`.

Bug fixes SHOULD include a regression test that fails without the fix and passes with it.

Tests MUST NOT need network access or private data, unless they are marked and skipped by default.

Run the tests for the area you changed before opening a pull request, and check that the tools you touched still run:

```bash
python -m pytest test/<area>
```

There are no required CI checks yet, so nothing else will run them for you. A test workflow will be added later and will then become required.

## Scientific correctness

Code that looks plausible can still give wrong biology. Changes that touch coordinates, indexing, identifiers or numerical results MUST state in the pull request how correctness was checked. Examples are comparing against a known output, a worked example from a paper, or a reference tool such as HMMER or BLAST.

Say which convention you rely on: 0-based or 1-based, inclusive or exclusive ends, strand handling.

## Dependencies

A new runtime dependency MUST be justified in the pull request and approved by a maintainer. Say what the package does that the standard library or an existing dependency does not, and how well maintained it is. Prefer what the `rotifer` environment already provides.

## Moving and renaming code

Moving or renaming a public module or function breaks every script that imports it. It SHOULD only happen when there is a concrete benefit, such as removing a real duplicate or fixing a name that misleads users.

When you do it, the old name MUST keep working through a deprecated alias that points to the new one. The mechanics and the removal schedule are in {doc}`compatibility`.

## AI-assisted contributions

Using an AI assistant to help write code, tests or documentation is allowed. The author is fully responsible for the result, exactly as if they had typed it. You MUST have read the code, understood it and tested it on real data before you submit it.

Generated code in a scientific codebase carries a specific risk: it is often plausible-looking and subtly wrong. It invents function arguments that do not exist, silently changes edge-case behaviour, misuses domain conventions such as coordinate systems, sequence indexing and taxonomic identifiers, and produces results that look reasonable enough to pass a quick glance and still be scientifically incorrect. A reviewer cannot reliably catch this if the author did not check it first.

Before you open the pull request:

- Read the code you are submitting and be able to explain why each function is there.
- Verify that the APIs, arguments and file formats it uses exist and behave as assumed. Do not trust an assistant's description of ROTIFER's own internals.
- Run it on real data and check the output, not only that it runs without an exception.
- Remove code that was generated "just in case" and is not needed.

Substantial AI assistance MUST be disclosed in the pull request description. Disclosure is not held against you. It tells reviewers where to look more carefully.

Reviewers MAY close a pull request that appears to be unreviewed generated output without reviewing it in detail.
