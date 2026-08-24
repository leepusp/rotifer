# Documenting code

ROTIFER docstrings follow the [NumPy docstring standard](https://numpydoc.readthedocs.io/en/latest/format.html).
The site parses them with **numpydoc only**. There is no napoleon
fallback: Google-style docstrings (`Args:`, `Returns:` with colons)
will render as unformatted text. If you find a Google-style
docstring, convert it.

numpydoc is not a replacement for autodoc. It preprocesses
NumPy-style sections into reStructuredText and autodoc still
renders the result, which is why both are enabled in `conf.py` and
why napoleon must never be added next to them: running two
docstring parsers produces duplicate Parameters sections.

## A fully annotated example

Every section that this project uses, in the order numpydoc
expects. Include a section only when it carries real information;
never emit an empty section.

```python
def fetchall(self, accessions, ipgs=None):
    """
    Fetch all gene neighborhoods as a single dataframe.

    The one-line summary above is imperative ("Fetch", not
    "Fetches") and fits on one line. This extended summary adds
    context: what the method is for and what makes it different
    from its siblings.

    Parameters
    ----------
    accessions : list of str
        NCBI protein accessions. The type after the colon matches
        the real signature. Optional parameters state their
        default, like the next one.
    ipgs : pandas.DataFrame, optional
        Precomputed identical protein group reports, used to avoid
        downloading IPGs again.

    Returns
    -------
    rotifer.genome.data.NeighborhoodDF
        The concatenated neighborhoods. Empty when nothing could
        be retrieved. Use ``Yields`` instead of ``Returns`` for
        generators.

    Raises
    ------
    ValueError
        Only list exceptions the caller can reasonably expect to
        handle.

    See Also
    --------
    fetchone : the lazy, one-result-at-a-time variant

    Notes
    -----
    Anything that does not fit the other sections: algorithmic
    details, rate limits, caching behavior.

    References
    ----------
    .. [1] Sayers E. "E-utilities Quick Start." NCBI Help Manual.
       https://www.ncbi.nlm.nih.gov/books/NBK25500/

    Examples
    --------
    Doctest format, runnable offline, or marked with a skip
    directive when it reaches a live database:

    >>> gnc = GeneNeighborhoodCursor()  # doctest: +SKIP
    >>> df = gnc.fetchall(["WP_012291365.1"])  # doctest: +SKIP
    """
```

Module docstrings state what biological resource the module talks
to and any rate limit, authentication or caching behavior the
reader needs.

## Good and bad, side by side

```{list-table}
:header-rows: 1

* - Bad
  - Good
* - `query: list of string` (type invented, colon glued to name)
  - `query : list of str` (space before the colon, real type)
* - `progress: boolean, deafult False`
  - `progress : bool, default False`
* - Summary line "This function fetches sequences."
  - Summary line "Fetch sequences from local or remote databases."
* - `Returns` section saying only "A dataframe"
  - `Returns` naming the type and what each column means
* - A `Usage` or `Rational` section (not part of the standard)
  - `Examples` and `Notes`, which numpydoc actually renders
* - Example code that cannot run (`from rotifer.db as ncbi`)
  - Copy-pasteable examples, skipped when they need the network
```

## Cross-referencing

- Link Python objects with `` {py:func}`rotifer.db.proteins` ``,
  `` {py:class}`rotifer.db.core.BaseCursor` `` and
  `` {py:meth}`~rotifer.db.core.BaseCursor.fetchall` `` (the tilde
  shows only the last name segment). Inside docstrings use the
  reST forms `` :func:`...` ``, `` :class:`...` ``,
  `` :meth:`...` ``.
- External libraries resolve through intersphinx: writing
  ``pandas.DataFrame`` in a type field, or `` :class:`python:dict` ``
  in text, links to the pandas and Python manuals. Configured
  targets: python, numpy, pandas, biopython.
- The build runs in nitpicky mode: a reference that does not
  resolve fails the strict build. Refer only to documented
  objects, or use plain literals (double backticks) for attribute
  names that have no documented target.

## When to use See Also

Use `See Also` to point at real alternatives: the faster variant
of the same operation, the delegator that wraps a backend, the
class that consumes this function's output. Do not list every
neighbor in the module; two or three entries the reader might
actually want are worth more than ten they must filter.

## Writing doctest examples

`make doctest` executes examples written inside explicit
`` ```{doctest} `` blocks (and, once the legacy docstrings are
converted, bare `>>>` examples can be enabled again through the
`DOCTEST_ONLY_EXPLICIT_BLOCKS` tunable in `conf.py`). Two rules
keep the target green:

1. Examples that run entirely offline stay executable. The build
   environment provides numpy and pandas for real; the heavy
   scientific dependencies (Biopython, ete3, tqdm and the rest of
   the mocked list in `conf.py`) are replaced by stubs.
2. Every line that reaches a live database endpoint carries
   `# doctest: +SKIP`.

The block below is a real `doctest` directive, so `make doctest`
executes it on every run. It doubles as proof that the offline
harness works: `rotifer.db.core` imports cleanly with the mocked
dependencies in place.

```{doctest}
>>> from rotifer.db.core import BaseCursor
>>> sorted(BaseCursor().parse_ids("a,b"))   # runs offline
['a', 'b']
>>> BaseCursor().parse_ids(["x", "x", "y"]) == {"x", "y"}
True
```

A line that reaches a live endpoint is marked instead of executed:

```python
>>> cursor.fetchall(["WP_063732599.1"])     # doctest: +SKIP
```

## Referencing a paper

Use a numbered reference in the `References` section and cite it
from the text as `[1]_`. Give authors, title, journal or resource
name, and a DOI or stable URL. Never invent a citation: if you do
not have the reference at hand, leave the section out and note the
gap in `docs/OPEN_QUESTIONS.md`.
