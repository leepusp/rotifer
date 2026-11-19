# Contributing

This section is for anyone who changes ROTIFER: laboratory members
working on their own analyses, maintainers, and outside contributors
sending patches. It holds the project's contribution rules, which
apply from 2026-12-01.

```{important}
From 2026-12-01 the `master` branch is protected. Every change goes
through a pull request with one approving review, and the personal
modules in `lib/rotifer/devel/alpha/` no longer accept changes. If
you have code there, read {doc}`personal-codebases` now.
```

## Start here

- {doc}`policy` is the binding text: what is required, what is
  recommended, and what happens in an emergency. Read it once.
- {doc}`workflow` walks through a change from issue to merge, with the
  commands to run.
- {doc}`code-standards` says where code goes and what a pull request
  needs before review: tests, docstrings, a note on how correctness
  was checked.
- {doc}`personal-codebases` explains the end of personal modules and
  where that code goes instead.

## Reference

- {doc}`reviewing` is the checklist for anyone approving pull
  requests.
- {doc}`compatibility` covers deprecations and version numbers.
- {doc}`governance` lists the maintainers and how decisions are made.
- {doc}`repository-settings` shows how the rules are enforced on
  GitHub, for admins.
- {doc}`sources` links each rule to the open-source projects and
  practices it was taken from.
- {doc}`documentation/index` covers writing docstrings, building the
  site locally, the design system that keeps it consistent, and a FAQ
  on which file controls which part of the site.

```{toctree}
:hidden:

policy
workflow
code-standards
personal-codebases
reviewing
compatibility
governance
repository-settings
sources
documentation/index
```
