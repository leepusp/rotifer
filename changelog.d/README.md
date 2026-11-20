# Changelog fragments

Each pull request that changes something a user can observe adds one small file here. At release time towncrier collects the files into `CHANGELOG.md` and deletes them.

## Naming

`changelog.d/<PR>.<type>.md`, where `<PR>` is the pull request number (not the issue number). For a second fragment of the same type, use `167.fixed.1.md`.

## Types

In the order they appear in a release:

- `breaking`: can break existing code or scripts; say what to do instead
- `results`: same input, different output; name the function or tool
- `added`: new features
- `changed`: non-breaking behaviour changes
- `deprecated`: name the replacement; removal comes in the next release
- `removed`: things removed after their deprecation release
- `fixed`: bug fixes that do not change results
- `security`: vulnerability fixes
- `dependencies`: dependency and Python support changes

## Example

`changelog.d/147.results.md`:

```markdown
`rotifer.seq.envelope` now returns 1-based, inclusive coordinates. Before this release the end coordinate was off by one on the reverse strand, so saved tables built with older versions should be regenerated.
```

No fragment is needed for tests, CI, internal refactors or docs. Say "No changelog: <reason>" in the pull request instead.

Full rules: https://leepusp.github.io/rotifer/contributing/changelog.html
