# Releasing

For maintainers who cut a release. MUST, SHOULD and MAY are used as defined in [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119).

Contributors add changelog fragments as they work, as described in {doc}`changelog`. The release adds up those fragments. The result is published on the {doc}`/changelog` page.

## When releases happen

Releases happen when ready, with no fixed schedule. Any maintainer MAY propose one by opening an issue titled "Release X.Y.Z".

## Version numbers and tags

Versions follow [SemVer](https://semver.org) and are written as [PEP 440](https://peps.python.org/pep-0440/) versions. What counts as a breaking change, and how long a deprecation lasts, is set out in [Compatibility and deprecation](https://leepusp.github.io/rotifer/contributing/compatibility.html).

Tags are `vMAJOR.MINOR.PATCH`, and pre-releases are `vX.Y.ZrcN`. The old tag style `v.0.1.0.a` is not used again.

The version is set in one place: `setup.py` today, and `pyproject.toml` once packaging moves there.

## Release steps

1. Create a release branch.

   ```bash
   git switch -c release/X.Y.Z
   ```

2. Set the new version in `setup.py`.

3. Build the changelog section. This writes `CHANGELOG.md` and deletes the fragments. Preview it first with `--draft`.

   ```bash
   towncrier build --draft --version X.Y.Z
   towncrier build --version X.Y.Z --yes
   ```

4. Add the Highlights to the new section, and the Contributors list. Get the contributors from the commits since the previous tag.

   ```bash
   git shortlog -sn vPREV..HEAD
   ```

5. Check the deprecations, as described in the next section.

6. Open the release pull request. It needs one maintainer approval. The reviewer checks that every `breaking`, `results` and `deprecated` entry says what to do, and that items deprecated in the previous release are removed or have a recorded extension.

   ```bash
   gh pr create --title "Release X.Y.Z"
   ```

7. After the merge, tag the merge result on `master` and push the tag.

   ```bash
   git switch master
   git pull
   git tag -a vX.Y.Z -m "ROTIFER X.Y.Z"
   git push origin vX.Y.Z
   ```

8. Create the GitHub release, using the new section of `CHANGELOG.md` as the notes.

   ```bash
   gh release create vX.Y.Z --title "ROTIFER X.Y.Z" --notes-file <section>
   ```

## Removing deprecations due in this release

A deprecated item spends exactly one release deprecated and is removed in the next one. Before the release, find what is due:

```bash
grep -rn "DeprecationWarning" lib bin
```

Compare the result with the Deprecated section of the previous release in `CHANGELOG.md`. Every item listed there MUST be removed in this release, and each removal gets a `removed` fragment.

Only in extraordinary circumstances MAY maintainers keep an item longer. The reason MUST be recorded in the Deprecated section of the release that keeps it.

## Yanking a release

A bad release is never re-tagged or deleted. Publish a PATCH release that fixes it, and mark the bad release in its `CHANGELOG.md` heading, with the reason below it:

```markdown
## [X.Y.Z] - YYYY-MM-DD [YANKED]
```

See [PEP 592](https://peps.python.org/pep-0592/) for the same idea in Python package indexes.
