# Workflow

For anyone sending a change to ROTIFER, from the first issue to the merge.
MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

## Open or find an issue

Non-trivial changes start with an issue: features, behaviour changes,
refactors, and anything over about 200 changed lines. Describe the problem
and the approach you have in mind, and wait for a maintainer or reviewer to
react before you write much code. A typo or a one-line fix MAY skip the
issue.

Search the existing issues first. If one already covers your change, comment
there instead of opening a second one.

## Create a branch

Branch from an up-to-date `master`. Use one of these prefixes and, when there
is an issue, put its number right after the prefix:

- `feature/` for new functionality
- `fix/` for bug fixes
- `refactor/` for changes that keep behaviour as it is
- `docs/` for documentation only
- `chore/` for tooling, packaging and housekeeping
- `test/` for tests only

Examples: `fix/142-envelope-coords`, `docs/contribution-policy`.

Lab members MAY branch inside leepusp/rotifer. Outside contributors MUST work
from a fork.

```bash
git switch master
git pull
git switch -c fix/142-envelope-coords
```

## Make commits

Commits SHOULD be atomic: each one builds and passes the tests on its own,
so `git bisect` stays useful and `git revert` stays safe. Changes that only
reformat code MUST be separate commits from changes that alter behaviour.

Write the message as an imperative summary of at most 50 characters, a blank
line, and a body wrapped at 72 characters that explains why the change is
needed. The diff already shows what changed. End with the issue reference.

```text
Fix envelope coordinates for reversed hits

The envelope start and end were swapped when the hit was on the
reverse strand, so downstream slicing returned the wrong residues.
Swap them back before building the interval.

Fixes #142
```

Commits MUST NOT contain secrets, files over 5 MB, generated outputs, or
notebook cell outputs. Clear the outputs of a notebook before you commit it.

## Keep the branch current

Before you ask for review, bring the branch up to date with `master`, by
rebase or by merging `master` into it. The diff then shows only your change.
A branch SHOULD be up to date with `master` before it is merged, although the
repository does not enforce this.

```bash
git fetch origin
git rebase origin/master
git push --force-with-lease
```

Force-pushing your own feature branch after a rebase is fine. Force-pushing
`master` is forbidden for everyone.

## Open the pull request

Push the branch and open a pull request against `master`.

```bash
git push -u origin fix/142-envelope-coords
gh pr create --base master --fill
```

The template at {ghsrc}`.github/PULL_REQUEST_TEMPLATE.md` loads when you
open the pull request. The description MUST state:

- what the change does and why
- how you tested it
- whether user-visible behaviour changed
- whether AI tools were used substantially

Keep to one concern per pull request. If you find an unrelated problem on
the way, open an issue or a second branch for it. The remaining requirements
(tests, docstrings, correctness statement, dependencies) are in
{doc}`code-standards`.

Draft pull requests are welcome for early feedback or work in progress.
Reviewers do not review them in depth, so mark the pull request as ready when
you want a full review.

## Review

A pull request needs one approving review from someone other than its
author, and all review threads MUST be resolved before merge. After new
commits are pushed, the latest push MUST be approved by someone other than
the person who pushed it. Answer questions in the thread, push fixes as new
commits, and mark each thread resolved once it is dealt with.

Reviewers: see {doc}`reviewing`.

## Merge

The repository allows two merge methods, squash and rebase. Merge commits are
not allowed, so the history of `master` stays linear.

- Rebase when every commit in the branch is clean and atomic.
- Squash when the branch has fixup commits, work-in-progress commits, or
  messages that would not read well in the history of `master`.

When the pull request is approved and the threads are resolved, the author
merges it if they have write access. Otherwise a reviewer merges it. Delete
the branch after the merge.

```bash
gh pr merge --squash --delete-branch
```

No required CI checks exist yet. A test workflow will be added later and then
become required. Until then, run the tests locally before you open the pull
request.
