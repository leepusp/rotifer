# Governance

Who can do what in rotifer and how decisions are made.

MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

## Roles

```{list-table}
:header-rows: 1
:widths: 18 32 50

* - Role
  - Who
  - What they can do
* - Contributor
  - Anyone
  - Open issues, comment, and send pull requests from a fork. Lab
    members MAY also branch inside `leepusp/rotifer`.
* - Reviewer
  - Any lab member with write access
  - Everything a contributor can do, plus approve pull requests
    (see {doc}`reviewing`). Reviewers MUST NOT approve their own
    pull requests.
* - Maintainer
  - The people listed below
  - Everything a reviewer can do, plus approve new runtime
    dependencies, approve policy changes, manage repository
    settings, and use the emergency bypass described in the
    {doc}`policy`.
```

## Maintainers

- leepbioinfo (Robson F. de Souza), project lead:
  https://github.com/leepbioinfo
- castlelli: https://github.com/castlelli
- EdwardNSeven: https://github.com/EdwardNSeven

## How decisions are made

Most decisions use lazy consensus. A proposal is made in an issue
or a pull request, and if nobody objects in a reasonable time, it
goes ahead. Silence counts as agreement. Anyone MAY object, and an
objection MUST come with a reason and, where possible, an
alternative.

Policy changes have a fixed objection window of 5 working days,
counted from the announcement to the lab. A policy change needs the
approval of two maintainers and cannot be merged before the window
closes.

When maintainers disagree and the discussion does not settle, the
project lead (leepbioinfo) decides.

## Adding and removing maintainers

Maintainers are added or removed by agreement of the existing
maintainers. The change is recorded by a pull request to this page,
so the list above always matches the repository settings.

## Code owners

The file {ghsrc}`.github/CODEOWNERS` asks GitHub to request a
review from the maintainers automatically when a pull request is
opened. Code-owner approval is not required to merge. The merge
rules that apply are the ones in the {doc}`policy`.
