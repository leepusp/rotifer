# Where these rules come from

For anyone who wants to know which established practice each rule follows.
The rules are adapted from other projects and from ROTIFER's own history.

| Rule | Precedent |
|---|---|
| MUST, SHOULD and MAY vocabulary, see {doc}`policy` | [RFC 2119](https://www.rfc-editor.org/rfc/rfc2119) (IETF), also used in Python PEPs |
| Roles, lazy consensus, project lead breaks ties, see {doc}`governance` | Apache Software Foundation lazy consensus; scikit-learn and NumPy governance documents |
| Changes reach `master` only through pull requests, see {doc}`policy` | CPython devguide (core developers also use pull requests); GitHub Flow |
| Every change reviewed by someone other than the author, one approval, see {doc}`reviewing` | [Google Engineering Practices code review guide](https://google.github.io/eng-practices/review/); one approval is the usual size for small teams such as SciPy and NumPy |
| Code owners, see {doc}`governance` | Kubernetes OWNERS files; GitHub CODEOWNERS |
| `master` always passes tests (planned required checks), see {doc}`policy` | Graydon Hoare's "not rocket science rule" (Rust) |
| Squash or rebase, linear history, see {doc}`workflow` | NumPy and Rust practice (bisect, clean reverts) |
| Discuss in an issue before large changes, see {doc}`workflow` | pandas and scikit-learn contributing guides |
| Atomic commits, 50/72 messages, see {doc}`workflow` | Linux kernel [Submitting patches](https://docs.kernel.org/process/submitting-patches.html) (separate each logical change); Tim Pope's commit message note |
| Reuse, docstrings, tests with features, regression tests, see {doc}`code-standards` | pandas and scikit-learn pull request checklists |
| Scientific correctness note in the PR, see {doc}`code-standards` | ROTIFER's own history (#142 coordinate fixes) |
| Deprecation cycle, see {doc}`compatibility` | [NumPy NEP 23](https://numpy.org/neps/nep-0023-backwards-compatibility.html) (backwards compatibility and deprecation policy); pandas deprecation policy |
| Versioning, see {doc}`compatibility` | [semver.org](https://semver.org) |
| AI disclosure and author responsibility, see {doc}`code-standards` | LLVM AI tool use policy; Fedora AI-assisted contributions policy; curl's stance on unverified AI reports |
| Emergency bypass, see {doc}`policy` | Google SRE "break-glass" practice; GitHub ruleset bypass lists |
| Experimental code outside the core, see {doc}`personal-codebases` | scikit-learn-contrib; astropy affiliated packages |
| Security reporting, see {doc}`policy` | GitHub SECURITY.md and private vulnerability reporting; OpenSSF Best Practices |
| Policy changes need more approval than code changes, see {doc}`policy` | NumPy governance |
| Changelog categories, newest release first, ISO dates, see {doc}`changelog` | [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/) |
| One fragment per pull request, collected at release, see {doc}`changelog` | [towncrier](https://towncrier.readthedocs.io), used by pytest, pip, attrs and Twisted; astropy and NumPy collect per-change release note files the same way |
| "No changelog" reason in the pull request, see {doc}`changelog` | pytest and astropy "skip changelog" labels |
| Changed results section, see {doc}`changelog` | scikit-learn release notes, "Changed models" section |
| Upgrade notes listed first, with what to do, see {doc}`changelog` | Django release notes, "Backwards incompatible changes" |
| Highlights and contributors in each release, see {doc}`releasing` | NumPy, SciPy and astropy release notes |
| Tags and pre-release names, yanking instead of re-tagging, see {doc}`releasing` | [PEP 440](https://peps.python.org/pep-0440/); [PEP 592](https://peps.python.org/pep-0592/) (yanked releases) |
| One release deprecated, removed in the next, see {doc}`compatibility` | pandas deprecation policy (deprecate in one release, remove in a later one), shortened for a small team |
| ADR format, numbering, immutability once accepted, see {doc}`adr` | [Michael Nygard, "Documenting Architecture Decisions"](https://cognitect.com/blog/2011/11/15/documenting-architecture-decisions) (2011) |
| Alternatives considered as a required ADR section, see {doc}`adr` | [MADR](https://adr.github.io/madr/) (Markdown Any Decision Records) |
| ADR status lifecycle, superseding, approval by maintainers, see {doc}`adr` | Python PEP 1; Rust RFCs; Kubernetes KEPs |
| Decision records in scientific software, see {doc}`adr` | NumPy NEPs; astropy APEs; scikit-learn SLEPs |
| Validation section for decisions that can change results, see {doc}`adr` | ROTIFER's equivalence and benchmark checks in refactors (#144, #151) |
| Conduct, see {doc}`policy` | [Contributor Covenant](https://www.contributor-covenant.org) |

## Where ROTIFER deviates on purpose

Three rules are looser than the practice they come from. The team is
small, so administrators have an emergency bypass for a broken
`master`, a security problem or a broken docs deploy, and every use has
to be explained within 2 working days. Required checks are absent
until a test workflow is added to the repository, and then they become
required. A change needs 1 approving review where larger projects
often ask for 2.
