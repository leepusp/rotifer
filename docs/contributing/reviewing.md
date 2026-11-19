# Reviewing pull requests

For reviewers: any lab member with write access. MUST, SHOULD and MAY are
used as defined in the {doc}`policy`. The author's side of the process is in
{doc}`workflow`.

## Response time

A reviewer SHOULD give a first response within 5 working days of a pull
request being ready for review. If you cannot review it, say so and suggest
someone else. Draft pull requests do not need a full review.

## What to check

Scope
: The pull request has one concern. If it mixes a fix with a refactor, ask
  for a split.

Issue
: A non-trivial change links to an issue where the approach was discussed.
  A pull request over about 200 changed lines without one needs a short
  discussion before the review goes further.

Tests
: New behaviour comes with tests under `test/`, mirroring the package path.
  A bug fix SHOULD include a regression test that fails without the fix.
  Check that the tests assert something meaningful and do not need network
  access or private data unless they are marked and skipped by default.
  Ask whether the author ran them locally.

Docstrings and docs
: Public functions, classes and command-line tools have NumPy-style
  docstrings ({doc}`documentation/documenting-code`). A change to the public
  API updates the documentation in the same pull request.

Reuse
: Search for existing functionality with `git grep -i <keyword>`. If the
  change duplicates code that exists, ask the author to extend it. New
  library code goes in `lib/rotifer/<package>`, command-line tools in `bin/`,
  and nothing new under `lib/rotifer/devel/`.

Scientific correctness
: Changes that touch coordinates, indexing, identifiers or numerical results
  MUST say in the description how correctness was checked, for example
  against a known output or a reference tool. Look at that statement and
  decide whether it convinces you.

Dependencies
: A new runtime dependency is justified in the description and approved by a
  maintainer. Do not approve it on your own.

Deprecation
: Moving or renaming a public module or function leaves a deprecated alias
  behind, and a removal follows the cycle in {doc}`compatibility`. Moves
  SHOULD have a concrete benefit.

AI use
: The description says whether AI tools were used substantially. Look for
  signs of generated code the author did not review: arguments that do not
  exist, changed edge-case behaviour, wrong coordinate or indexing
  conventions, code added "just in case", comments that describe what the
  code does line by line, and results that look plausible but were never
  compared with real data. Ask the author to explain a function if in doubt.
  A reviewer MAY close a pull request that looks like unreviewed generated
  output without a detailed review.

## How to comment

Ask questions before you assert. "What happens when the hit is on the
reverse strand?" gets a better answer than "this is wrong", and sometimes
the author has a reason you did not see.

Say whether each comment blocks the merge. Prefix suggestions that the
author may take or leave with `nit:`, and keep them few. A comment without
the prefix is understood to need an answer or a change.

Comment on the code and the change, and give the reason for a requested
change. Follow the Contributor Covenant in `CODE_OF_CONDUCT.md`.

## Approve or request changes

Approve when the change meets the points above and you would be
comfortable maintaining the code.
Request changes when something blocks the merge. Say what is needed so the
author knows what a second look will check.

Approval applies to the commits you saw. When the author pushes new commits,
the latest push MUST be approved again by someone other than its pusher. If
you pushed fixes to the branch yourself, another reviewer approves them.

You MAY approve with `nit:` comments still open and trust the author to
handle them. The author still resolves every thread, `nit:` ones
included, before merging: by making the change or by replying why not.
GitHub refuses the merge while any thread is open.

## Personal modules

From 2026-12-01 reviewers MUST NOT approve pull requests that add to or
change personal modules in `lib/rotifer/devel/alpha/`. The one exception is a
maintainer's migration pull request. Point the author to
{doc}`personal-codebases` for where the code should go instead.

## Merging

Once a pull request has its approval and no open threads, the author merges
it if they have write access. Otherwise a reviewer merges it.

Choose the method when you merge. Rebase if every commit is clean and
atomic. Squash if the branch has fixups or messages that would clutter the
history, and check that the squashed message follows the 50/72 format and
carries the issue reference. Merge commits are not allowed. Delete the branch
afterwards.
