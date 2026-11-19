# Contribution policy

The rules for anyone who changes code, tests, documentation or
configuration in rotifer. Read it before your first pull request.

The key words MUST, MUST NOT, SHOULD, SHOULD NOT and MAY in this
documentation are to be interpreted as described in
[RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

## Scope and effective date

The policy covers everything in the `leepusp/rotifer` repository:
library code, tests, documentation, scripts in `bin/` and
repository configuration. It applies to maintainers and to every
other contributor.

It takes effect on 2026-12-01. Until then current practice
continues, and contributors SHOULD start following it now.

## Rules at a glance

1. Non-trivial changes start with an issue. See {doc}`workflow`.
2. Work happens on a branch with a conventional prefix, and each
   pull request has one concern. See {doc}`workflow`.
3. Every pull request description states what changed and why, how
   it was tested, whether user-visible behaviour changed and whether
   AI tools were used substantially. See {doc}`workflow`.
4. Commits are atomic, formatting-only changes are separate commits,
   and messages follow the 50/72 form. See {doc}`workflow`.
5. Nothing reaches `master` without a reviewed pull request. See
   {doc}`reviewing`.
6. Library code goes in `lib/rotifer/<package>`, existing
   functionality is reused, and public code has NumPy-style
   docstrings. See {doc}`code-standards`.
7. New behaviour comes with tests, and changes that touch
   coordinates, indexing or numerical results say how correctness
   was checked. See {doc}`code-standards`.
8. Moving or removing public code goes through a deprecation cycle:
   one release deprecated, removed in the next.
   See {doc}`compatibility`.
9. Personal modules in `lib/rotifer/devel/alpha` are frozen on the
   effective date. See {doc}`personal-codebases`.
10. AI tools are allowed, the author is fully responsible for the
    result, and substantial use is disclosed. See
    {doc}`code-standards`.
11. Every user-visible change adds a changelog fragment in the same
    pull request. See [Changelog](https://leepusp.github.io/rotifer/contributing/changelog.html).
12. Changes to packages, persistent formats, cross-module conventions,
    central dependencies or project infrastructure need an accepted
    architecture decision record before they merge. See [Architecture decision records](https://leepusp.github.io/rotifer/contributing/adr.html).

Roles and decision making are in {doc}`governance`, and the
precedent behind each rule is in {doc}`sources`.

## The master branch

The default branch is `master`. It is protected by a ruleset with
these rules:

- Changes arrive only through a pull request. Direct pushes are not
  allowed.
- A pull request MUST have 1 approving review from someone other
  than its author.
- The most recent push MUST be approved by someone other than the
  person who pushed it. New commits after an approval require a new
  approval.
- All review threads MUST be resolved before merge.
- Merging uses squash or rebase. Merge commits are not allowed, so
  the history stays linear. Rebase when every commit is clean and
  atomic, and squash otherwise.
- The branch SHOULD be up to date with `master` before merge. This
  is not enforced.
- Required status checks: none yet. A test workflow will be added
  later and will then become a required check. Until then, run the
  tests locally before opening a pull request.

Maintainers configure and audit these settings as described in
{doc}`repository-settings`.

## Force-push and deletion

Force-pushing to `master` and deleting `master` are forbidden for
everyone. There are no exceptions, and the emergency bypass below
does not cover them.

## Emergency bypass

Maintainers with the admin role MAY push directly to `master`, and
only for one of these reasons:

- `master` is broken.
- There is a security problem.
- The documentation deploy is broken.

Within 2 working days the maintainer MUST open an issue or a
follow-up pull request that explains what was changed and why.
Using the bypass for convenience is a policy violation.

## Security

Do not report vulnerabilities in public issues. Use GitHub private
vulnerability reporting (Security tab, then Report a
vulnerability). The details are in {ghsrc}`SECURITY.md`. Reports
are acknowledged within 7 days.

## Conduct

All participation follows the Contributor Covenant. The full text
and the reporting process are in {ghsrc}`CODE_OF_CONDUCT.md`.

## Changing this policy

A change to the policy is a pull request against the policy pages
in `docs/contributing/`. It MUST be approved by two maintainers
and announced to the lab. The objection window for policy changes
is described in {doc}`governance`.
