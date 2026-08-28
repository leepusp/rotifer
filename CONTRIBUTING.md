# Contributing to ROTIFER

Thank you for your interest in contributing to ROTIFER. This project is a  
collection of libraries and command line tools for comparative
genomics and the computational analysis of biological sequences, and it grows
mostly through the work of people who use it in their own research.

There are many ways to help: reporting bugs, improving documentation, suggesting
features, and submitting code. This document explains how to do each of these in
a way that is easy for maintainers to review.

## Code of Conduct

This project has adopted the [Contributor Covenant](CODE_OF_CONDUCT.md). By
participating, you are expected to uphold it.

## Asking questions

If you have a question about how to use ROTIFER, please open an issue with
the *question* label rather than emailing the maintainers directly. Answers in
public benefit everyone who runs into the same doubt later.

## Reporting issues

### Before you open an issue

**Search first.** Please look through the
[open and closed issues](https://github.com/leepusp/rotifer/issues)
before creating a new one. If you find an existing report that matches your
problem:

- Add a comment only if you have new information: a different platform, a
  smaller reproduction case, a stack trace the original report lacks.

Duplicate issues are not a disaster, but consolidating a discussion in one place
makes it far more likely that the problem actually gets fixed.

### Writing a good bug report

Please use the issue forms. They ask for the information we need. The most
useful reports include:

- The version or commit of ROTIFER you are running.
- Your operating system and how you installed the environment (conda, and which
  package versions if relevant).
- A minimal, self-contained example that triggers the problem. A short script or
  a few lines in IPython is much easier to debug than a description of a
  pipeline.
- What you expected to happen, and what actually happened, including the full
  error message and traceback.

If the bug involves remote data (NCBI queries, for example), please mention
whether it is reproducible or intermittent.

## Suggesting features

### Check whether it already exists

ROTIFER is large and spread across several usage cases, so functionality is not
always where you would expect it. **Before writing something from scratch, please
check whether it already exists.** 

Searching the repository (`git grep -i <keyword>`) is often faster than browsing.
If something similar exists but does not quite fit your needs, extending it is
usually preferable to adding a parallel implementation, and it is a very welcome
kind of contribution.

If you are unsure, open an issue and ask before investing time in an
implementation. This is especially worthwhile for larger changes: it is much less
frustrating to discuss the design first than to have a finished pull request
turned down.

In order to make this matter, and contribution and usage  easier, our complete documentation of ROTIFER is on the way.

## Contributing code

### Setting up a development environment

Follow the installation instructions in the [README](README.md) to create the
`rotifer` conda environment and put the libraries and tools on your `PYTHONPATH`
and `PATH`. Working from a clone with `conda-develop` means your edits take
effect immediately, without reinstalling.

Before opening a pull request, please run the tests in `test/` and make sure the
tools you touched still run.

### Branches

Work on a branch in your fork rather than on `master`. Some suggestions that
tend to make review easier (treat them as guidance, not as rules):

- Keep a branch focused on one thing. A branch that fixes a parsing bug *and*
  refactors the logging setup is harder to review, and harder to revert if
  something turns out to be wrong.
- Give the branch a descriptive name, for example `fix/ipg-cursor-batch-size` or
  `doc/install-instructions`.
- Rebase or merge `master` into your branch before opening the pull request, so
  that the diff reflects only your changes.

### Commits

Likewise, atomic commits — each one a single, coherent, self-contained change —
make history much more useful. They make `git bisect` meaningful, they make
`git revert` safe, and they let a reviewer follow your reasoning step by step
instead of reading one large diff.

In practice this means:

- Avoid mixing formatting changes with behavioural changes in the same commit.
  If you need to reformat a file, do it in its own commit.
- Write a short imperative summary line (around 50 characters), a blank line,
  and then a body explaining *why* the change is being made if it is not obvious.
  What changed is visible in the diff; the reasoning usually is not.
- Reference the relevant issue in the body, for example `Fixes #42`.

Nobody is going to reject a contribution over commit granularity. But if you can
split your work into clean commits, please do. It genuinely speeds up review.

### Pull requests

- Open the pull request against `master`.
- Fill in the pull request template, including the checklist.
- Describe what the change does and how you tested it. If it changes behaviour
  that users depend on, say so explicitly.
- Be ready for review comments. Questions about a change are not a judgement of
  your work; they are how the project stays maintainable.

By submitting a pull request, you agree that your contribution is licensed under
the project's [BSD 3-Clause License](LICENSE).

## A note on AI-assisted contributions

Using an AI assistant to help write code, tests, or documentation is fine. Many
of us do. But please **do not submit code you have not read, understood, and
tested yourself.**

Generated code in a scientific codebase carries a specific risk: it is often
plausible-looking and subtly wrong. It invents function arguments that do not
exist, silently changes edge-case behaviour, misuses domain conventions
(coordinate systems, sequence indexing, taxonomic identifiers), and produces
results that look reasonable enough to pass a quick glance and still be
scientifically incorrect. A reviewer cannot reliably catch this if the author did
not check it first.

So, before you open the pull request:

- Read the code you are submitting and be able to explain why each function is there.
- Verify that the APIs, arguments, and file formats it uses actually exist and
  behave as assumed, do not trust the assistant's description of ROTIFER's own
  internals.
- Run it on real data and check the output, not just that it executes without
  an exception.
- Remove code that was generated "just in case" and is not needed.

Please mention in the pull request description if a substantial part of the
change was AI-assisted. This is not held against you; it simply tells reviewers
where to look more carefully. Pull requests that appear to be unreviewed
generated output may be closed without detailed review.

## If you have a personal codebase

- If you have a personal codebase inside ROTIFER, avoid using it whenever possible. Use a branch instead.
- If you do use it, pay extra close attention to the code you submit, especially:
  - Creation of functions that already exist on ROTIFER
  - AI-generated code
- Your codebase will eventually be fully merged into ROTIFER. Useful code will be adapted and duplicates deleted.
- Contact the maintainers if you have any questions about this.

## Documentation

A full, detailed documentation of ROTIFER is currently being built. Please check it for further information about contribution and/or the tool itself.

## Thank you

Your contributions are what
keep this project useful. We appreciate the time you put into them.