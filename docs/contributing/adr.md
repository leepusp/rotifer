# Architecture decision records

For contributors proposing a significant change and for maintainers
approving it. MUST, SHOULD and MAY are used as defined in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

An architecture decision record (ADR) is a short numbered document
that records one significant decision: the context, the decision,
the alternatives and the consequences. ADRs are published on the
site and listed in {doc}`/adr/index`. The format follows the idea in
[Documenting architecture decisions](https://cognitect.com/blog/2011/11/15/documenting-architecture-decisions);
[MADR](https://adr.github.io/madr/) is a similar template.

## When an ADR is required

A change MUST have an accepted ADR before its implementation merges
when it does any of the following.

1. Adds, removes, splits or merges a package under `lib/rotifer/`,
   for example removing `rotifer.cluster` or `rotifer.alchemy`.
2. Defines or changes a persistent format that other code or users
   depend on: on-disk layouts, database schemas and standard output
   schemas, for example session storage or the search output schema
   of #141.
3. Sets a convention used across modules: coordinate systems and
   indexing, identifier handling, logging, configuration, error
   reporting or command-line conventions.
4. Adds or replaces a central dependency or external service that
   many modules or users would rely on, for example ClickHouse or a
   storage engine. Ordinary dependencies only need maintainer
   approval, as described in [Code standards](https://leepusp.github.io/rotifer/contributing/code-standards.html).
5. Redesigns public API across several modules, or makes breaking
   changes that need migration across the codebase.
6. Changes project infrastructure: packaging and build, the release
   process, the documentation toolchain or the CI architecture.
7. Chooses between substantial alternatives where reversing the
   choice would be costly, for example translating all Perl code to
   Python.

## When an ADR is recommended

An author SHOULD write an ADR when a reviewer asks "why was it done
this way?" and the answer does not fit in a code comment.

## When an ADR is not needed

An ADR is not needed for bug fixes, refactors inside one module that
keep its interface, new functions that follow existing patterns,
performance work with unchanged interfaces, or documentation.
Changes to the contribution policy go through the policy pages, not
through an ADR.

As a rule of thumb, if someone two years from now could reasonably
undo the change without knowing why it was made, write an ADR.

## What an ADR contains

Every ADR MUST have a header with its status, date, deciders, issue
and pull request, and these sections:

- Context: the problem and the forces behind it. Facts, constraints
  and what breaks today. No solution yet.
- Decision: what will be done, in active voice, for example "ROTIFER
  stores sessions as ...". It MUST be specific enough for a reviewer
  to check the implementation against it.
- Alternatives considered: each serious option and why it was not
  chosen. Doing nothing counts as an option.
- Consequences: what becomes easier and harder, the compatibility
  impact and any deprecations, the changelog impact and follow-up
  work.

A Validation section MUST be included when the decision can change
scientific results. It says how the decision is or will be shown to
work: benchmarks, equivalence checks against the old behaviour, or
comparison with reference tools. A References section is optional.

An ADR SHOULD fit in about two pages.

## How to write one

1. Copy {ghsrc}`docs/adr/0000-template.md` to
   `docs/adr/NNNN-short-title.md`.
2. Use the next free four-digit number and a kebab-case title. Numbers
   are sequential and never reused. If two open pull requests pick
   the same number, the one merged second renumbers.
3. Fill in the header and sections, and delete the template comments.
4. Open a pull request that adds the file, with status `Proposed`,
   and link it from the issue.

An ADR MAY ship in the same pull request as its implementation. For
large changes it SHOULD come first, as its own pull request.

## Lifecycle

```{list-table}
:header-rows: 1

* - Status
  - Meaning
* - Proposed
  - The pull request adding the ADR is open.
* - Accepted
  - The pull request was approved by two maintainers and merged.
* - Rejected
  - The proposal was declined. The file is still merged, with the
    reason.
* - Deprecated
  - The decision no longer applies and nothing replaced it.
* - Superseded by ADR NNNN
  - A newer ADR replaces the decision.
```

An ADR becomes Accepted when its pull request is approved by two
maintainers. Implementation MUST NOT merge before that approval. The
approving maintainers are listed in the Deciders line.

Accepted ADRs are immutable except for the status line, links and
typos. To change a decision, write a new ADR that supersedes the old
one, and change the old status to "Superseded by ADR NNNN".

## Relation to the changelog

An ADR that changes user-visible behaviour still needs changelog
fragments, as described in {doc}`changelog`. Changelog entries MAY
link the ADR.
