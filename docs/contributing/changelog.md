# Changelog

For contributors who open pull requests. MUST, SHOULD and MAY are used as defined in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

ROTIFER follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Each pull request adds a small fragment file, and the release manager assembles the fragments into the changelog with [towncrier](https://towncrier.readthedocs.io). Maintainers cutting a release should read {doc}`releasing`. The published history is on the {doc}`/changelog` page.

## Files involved

- `changelog.d/<PR>.<type>.md` holds one fragment per change. It is added in the pull request that makes the change.
- {ghsrc}`CHANGELOG.md` at the repository root holds the full history, newest first. It is generated at release time by towncrier. Past releases are never edited by hand, except to mark a release as yanked.
- `docs/changelog.md` includes `CHANGELOG.md` so the documentation site shows it.
- {ghsrc}`towncrier.toml` defines the fragment types and the output format.

towncrier is a development tool and not a runtime dependency. Only the release manager needs it (`pip install towncrier`).

## When a fragment is required

A pull request MUST add a fragment when it makes a change that a user can observe:

- public functions, classes or modules that are added, changed, deprecated or removed;
- command line tools, their options or their defaults;
- output formats, column names or file layouts;
- results, meaning the same input gives different numbers, coordinates, hits or identifiers, including bug fixes;
- runtime dependencies, minimum versions or supported Python versions;
- security fixes.

A fragment is not required for tests only, CI, internal refactors with identical behaviour, contribution docs and typos. In that case the pull request description MUST say "No changelog: <reason>". The reviewer checks that there is either a fragment or a reason.

## Naming and types

Name the file `changelog.d/<PR>.<type>.md`, where `<PR>` is the pull request number. Always use the pull request number, never the issue number. Open the pull request first, then add the fragment.

The types are listed in the order they appear in a release.

| Type | Section heading | Use |
| --- | --- | --- |
| `breaking` | Upgrade notes | Can break existing code or scripts. MUST say what to do instead. |
| `results` | Changed results | Same input, different output. MUST name the function or tool and how results differ. |
| `added` | Added | New features. |
| `changed` | Changed | Non-breaking behaviour changes. |
| `deprecated` | Deprecated | MUST name the replacement and say it will be removed in the next release. |
| `removed` | Removed | Things removed after their deprecation release. |
| `fixed` | Fixed | Bug fixes that do not change results. |
| `security` | Security | Vulnerability fixes. |
| `dependencies` | Dependencies | New, removed or re-pinned dependencies, and Python support. |

A pull request MAY add several fragments, for example `167.added.md` and `167.deprecated.md`. A change that is both breaking and results-changing gets both fragments. For two fragments of the same type, number the second one: `167.fixed.md` and `167.fixed.1.md`.

## Writing entries

- Write for users and describe what changed from their point of view.
- Use one or two sentences in the present tense.
- Put names in code formatting, with the full import path the first time.
- Leave the pull request number out of the text. towncrier adds the link.
- For `breaking`, `results` and `deprecated` entries, say what the user should do.
- An entry MAY link an ADR (see {doc}`adr`).

For example, `changelog.d/147.results.md` could contain:

```markdown
`rotifer.seq.envelope` now returns 1-based, inclusive coordinates. Before this release the end coordinate was off by one on the reverse strand, so saved tables built with older versions should be regenerated.
```

## Preview locally

```bash
pip install towncrier
towncrier build --draft --version X.Y.Z
```

The draft prints the section that would be written and changes no files.

## Layout of a release section

Each release in `CHANGELOG.md` has this layout. Dates are ISO dates, empty sections are left out, and every entry ends with a link to its pull request.

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Highlights
### Upgrade notes
### Changed results
### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security
### Dependencies
### Contributors
```

Highlights has one to five bullets written by the release manager, and MAY be skipped for a PATCH release. Contributors lists everyone with a commit in the release.
