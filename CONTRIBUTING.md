# Contributing to ROTIFER

Thank you for your interest in contributing to ROTIFER. The full
contribution rules live on the documentation site:

**https://leepusp.github.io/rotifer/contributing/**

This file is a summary. Where it and the site disagree, the site wins.

## From 2026-12-01

- `master` is protected. Every change goes through a pull request
  with one approving review from someone other than the author.
- Pull requests merge by squash or rebase. History stays linear.
- The personal modules in `lib/rotifer/devel/alpha/` no longer accept
  changes. Exploratory code moves to Jupyter notebooks kept outside
  this repository, and reusable functions come in as pull requests.
  See [Personal codebases](https://leepusp.github.io/rotifer/contributing/personal-codebases.html).

## Before you start

- Open or find an issue for anything larger than a small fix, and
  search the code for existing functionality (`git grep -i <keyword>`)
  before writing something new.
- Work on a branch named with a prefix and the issue number, such as
  `fix/142-envelope-coords`. Outside contributors work from a fork.
- Keep a pull request to one concern, write atomic commits, and fill
  in the pull request template.
- Add tests and NumPy-style docstrings for new public code, and run the
  tests locally.
- AI tools are allowed, but you must have read, understood and tested
  everything you submit, and you must say in the pull request if a
  substantial part was AI-assisted.

The details are in the
[Workflow](https://leepusp.github.io/rotifer/contributing/workflow.html)
and [Code standards](https://leepusp.github.io/rotifer/contributing/code-standards.html)
pages.

## Questions, bugs and features

Use the issue forms. Search the
[open and closed issues](https://github.com/leepusp/rotifer/issues)
first. Report security problems privately as described in
[SECURITY.md](SECURITY.md), never in a public issue.

## Conduct and licence

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
By submitting a pull request, you agree that your contribution is
licensed under the project's [BSD 3-Clause License](LICENSE).
