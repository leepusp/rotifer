# Personal codebases

For anyone who keeps code in a personal module under `lib/rotifer/devel/alpha/`.
MUST, SHOULD and MAY are used as defined in the {doc}`policy`.

```{important}
On 2026-12-01 the personal modules are frozen. From that date the
repository no longer accepts changes to them. Until then current
practice continues, and owners SHOULD start moving their code now.
```

## Why the practice ends

Personal modules such as `epsoares.py`, `gian_func.py`, `malu.py`,
`trsantos.py`, `mvroliveira.py`, `rodolfo.py`, `rfsouza.py` and `igem.py`
grew inside `lib/rotifer/devel/alpha/` because that folder accepted code
regardless of quality and testing. Over time this produced three problems.

- Several functions exist in more than one module, and some repeat what
  the library packages already provide.
- Other people import code that nobody reviewed, so a change in one
  personal file can alter someone else's results.
- The code has no tests and no documentation, which makes it hard to
  tell what it is meant to do or whether it is correct.

## Timeline

| Date | What happens |
|---|---|
| Now to 2026-11-30 | Owners move exploratory code out of rotifer and tell the maintainers, through an issue, which functions they rely on and think belong in the library. |
| 2026-12-01 | `master` is protected. Changes to personal modules are no longer accepted and the files are frozen as they stand. Reviewers MUST NOT approve pull requests that add to or change personal modules, except migration pull requests from maintainers. |
| After 2026-12-01 | Maintainers merge useful code into the proper packages (adapted, tested and documented), delete duplicates and remove the personal modules. |

## What owners should do before 2026-12-01

1. Copy exploratory code and analyses into Jupyter notebooks that live
   outside the rotifer repository, in your own repository or a project
   folder. See the next section for the setup.
2. Open an issue with the label `refactor`. List the functions you use
   and think belong in the library, and say where each one is today.
3. Stop adding code to your personal module. New reusable functions go
   through a normal pull request, as described in {doc}`workflow` and
   {doc}`code-standards`.

A helper that lives in your notebook is fine. When other people need
it, propose it as a pull request.

## Setting up a notebook outside the repository

Create a project folder anywhere except inside your rotifer checkout,
activate the environment and start Jupyter there:

```bash
mkdir ~/projects/my-analysis
cd ~/projects/my-analysis
conda activate rotifer
jupyter lab
```

In the first cell, import rotifer as any user would:

```python
import rotifer
```

From there you work with the installed library.

## What happens to the files after the cutoff

The personal modules stay in the repository, unchanged, until the
maintainers have gone through them. Functions that other people use
are moved into the proper packages under `lib/rotifer/<package>`. Every
removed function that had users gets a deprecated alias for one
release, following {doc}`compatibility`. Everything else is deleted
together with the module.

## The retired alpha and beta procedure

The procedure in {ghsrc}`lib/rotifer/devel/README.md` (add code to
alpha, discuss it, move it to beta, then to its final place) is also
retired for new code. Nothing new goes under `lib/rotifer/devel/`,
neither in `alpha` nor in `beta`. New code goes into a package under
`lib/rotifer/` or, for command-line tools, into `bin/`, with an issue
first when the change is not trivial.

Questions about any of this can go to the maintainers listed in
{doc}`governance`.
