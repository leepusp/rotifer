# Repository settings

For maintainers with the admin role. The rules in the {doc}`policy`
are enforced by GitHub settings, and those settings are kept as files
in the repository so that any change to them goes through review like
everything else.

MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

## What is enforced, and where

```{list-table}
:header-rows: 1

* - Rule
  - Enforced by
* - Changes reach `master` only through a pull request, with one
    approval from someone other than the author
  - {ghsrc}`.github/rulesets/master-protection.json`
* - The most recent push to a pull request is approved by someone
    other than its pusher
  - same file, `require_last_push_approval`
* - All review threads are resolved before merging
  - same file, `required_review_thread_resolution`
* - Squash or rebase only, linear history
  - same file, `allowed_merge_methods` and `required_linear_history`;
    merge commits are also turned off in the repository settings
* - Admins can push directly in an emergency
  - same file, bypass for the `admin` repository role
* - Nobody can force-push to or delete `master`
  - {ghsrc}`.github/rulesets/master-integrity.json`, which has no
    bypass list
* - Maintainers are asked to review automatically
  - {ghsrc}`.github/CODEOWNERS`
* - Vulnerabilities are reported privately
  - private vulnerability reporting, see {ghsrc}`SECURITY.md`
```

Two rules are kept in separate rulesets on purpose. A bypass in a
GitHub ruleset applies to every rule in it, so the force-push and
deletion rules live in a ruleset with an empty bypass list. That way
the emergency bypass lets an admin push a commit but never rewrite
history.

There are no required status checks yet. A test workflow will be added
in a later pull request, and its job will then be added to
`master-protection.json` as a `required_status_checks` rule. Until
then, contributors run the tests locally (see {doc}`code-standards`).

## Applying the settings

The settings take effect on 2026-12-01. An admin runs the commands
below from a clone of the repository at that date, with the GitHub CLI
logged in as an account that has the admin role on `leepusp/rotifer`.

Create both rulesets:

```bash
gh api -X POST repos/leepusp/rotifer/rulesets \
  --input .github/rulesets/master-protection.json
gh api -X POST repos/leepusp/rotifer/rulesets \
  --input .github/rulesets/master-integrity.json
```

Turn off merge commits and delete branches once their pull request is
merged:

```bash
gh api -X PATCH repos/leepusp/rotifer \
  -F allow_merge_commit=false \
  -F allow_squash_merge=true \
  -F allow_rebase_merge=true \
  -F delete_branch_on_merge=true
```

Turn on private vulnerability reporting:

```bash
gh api -X PUT repos/leepusp/rotifer/private-vulnerability-reporting
```

Check the result:

```bash
gh api repos/leepusp/rotifer/rulesets --jq '.[] | "\(.id) \(.name) \(.enforcement)"'
gh api repos/leepusp/rotifer/rules/branches/master --jq '.[].type'
gh api repos/leepusp/rotifer/private-vulnerability-reporting
```

## Changing the settings

The JSON files are the source of truth. To change a rule, open a pull
request that edits the file, get it approved as described in
{doc}`governance`, and after merging apply it to the existing ruleset
by its ID:

```bash
gh api repos/leepusp/rotifer/rulesets --jq '.[] | "\(.id) \(.name)"'
gh api -X PUT repos/leepusp/rotifer/rulesets/<id> \
  --input .github/rulesets/master-protection.json
```

A setting changed only in the GitHub web interface drifts from the
file, and the next `PUT` silently reverts it. Maintainers SHOULD NOT
edit rulesets in the web interface.

## Who can bypass

The bypass list names the `admin` repository role. At the time of
writing the admins are leepbioinfo and castlelli. EdwardNSeven and
andreaniml are maintainers with the `write` role, so they cannot use
the emergency bypass until an admin changes their role in the
repository's collaborator settings. Any use of the bypass is subject to the
follow-up rule in the {doc}`policy`.
