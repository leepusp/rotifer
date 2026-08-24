# Open questions and known gaps

Working notes for the documentation. Each entry names the file and
line of the code it refers to. None of these can be resolved from
the docs side without editing executable code, which is out of
scope for documentation changes.

## Modules excluded from the API reference (import failures)

- `lib/rotifer/db/neighbors.py:131`: unterminated f-string
  expression (`{cursor.fetchone()[0]` without the closing brace)
  makes the module unimportable (SyntaxError), so autodoc cannot
  document it. Excluded from `docs/api/index.rst`. Docstrings have
  since been written for all of its definitions, so it is ready to
  be added to the reference once the syntax error is fixed.
- `lib/rotifer/db/local/core.py:14`: calls `loadConfig` without
  importing it, so importing the module raises NameError. Excluded
  from `docs/api/index.rst`.
- `lib/rotifer/io/base.py:6`: unclosed parenthesis (SyntaxError).
  Excluded.
- `lib/rotifer/alchemy/db_classes.py`: uses `MetaData` without
  importing it (NameError at import). Excluded.
- `lib/rotifer/seq/` (entire subtree) and
  `lib/rotifer/tools/search.py`: import the module
  `rotifer.table`, which does not exist in the repository
  (ModuleNotFoundError). Excluded.
- `lib/rotifer/devel/`: personal sandbox modules with tab
  indentation, Google-style docstrings and broken imports.
  Deliberately not documented.

## Legacy docstrings with malformed reST

The strict build suppresses the `docutils` warning category
(`SUPPRESS_LEGACY_DOCSTRING_WARNINGS` in `docs/conf.py`) because
these docstrings emit warnings that only a code edit can fix:

- `lib/rotifer/core/cli.py` (`version`, `parser.add`,
  `action.autoload`): stray indentation and short title
  underlines.
- `lib/rotifer/core/functions.py` (`vmsg`): definition list ends
  without a blank line.
- `lib/rotifer/genome/data.py` (`NeighborhoodDF.vicinity`,
  `NeighborhoodDF.series_to_compact_frequency`): block quote and
  indentation problems.
- `lib/rotifer/genome/database.py` and
  `lib/rotifer/genome/db/clickhouse.py` (`submit`): indentation
  problems and short title underlines.
- `lib/rotifer/genome/utils.py` (`seqrecords_to_dataframe`):
  unexpected indentation.
- `lib/rotifer/io/hhsuite.py` (`read_hhr`, `parse_hhr`):
  unbalanced `*` characters parsed as emphasis markers.
- `lib/rotifer/pandas/functions.py` (`print_everything`):
  unexpected indentation.

Once a file is fixed, rebuild with
`SUPPRESS_LEGACY_DOCSTRING_WARNINGS = False` to confirm it is
clean.

## Doctests disabled for bare example blocks

Legacy docstrings contain bare `>>>` blocks that are not runnable;
some are not even valid Python (for example the unclosed
`print(SeqIO.write(...)` in `lib/rotifer/db/ncbi/__init__.py:83`).
`DOCTEST_ONLY_EXPLICIT_BLOCKS` in `docs/conf.py` therefore limits
`make doctest` to explicit `.. doctest::` directives. When the
docstrings are rewritten to the NumPy standard (runnable examples,
`# doctest: +SKIP` on lines that reach live endpoints), flip the
toggle so docstring examples are tested again.

## Behavior unclear from the code

- `lib/rotifer/db/__init__.py:68` (`proteins`): after the first
  retrieval method, `targets = cursor.missing` assigns a pandas
  DataFrame (the `missing` property) where a set of accessions is
  expected by the next iteration. The fallback path therefore
  looks broken; the docs describe the evident intent (pass
  unresolved identifiers to the next method) without promising the
  current code does it.
- `lib/rotifer/db/ncbi/__init__.py:530` (GeneNeighborhoodCursor):
  the `save` backend is registered as a reader but the writer
  registration is commented out, so `save` reads existing SQLite3
  stores and does not write new ones. Documented as read only.
- `lib/rotifer/db/sql/sqlite3.py` (`BaseSQLite3Cursor.stored`):
  the `column` parameter defaults to the string `'block_id'` but is
  iterated with `for col in column`, expecting a list. The only
  caller that relies on the default,
  `GeneNeighborhoodCursor.insert`'s `self.stored(data)`, therefore
  iterates over the characters of `'block_id'` instead of the
  column name and would raise `KeyError` on the first character.
  Flagged with a `Note` in `stored`'s docstring rather than fixed.
- `lib/rotifer/db/uniprot/__init__.py` and
  `lib/rotifer/db/uniprot/webapi/__init__.py`: both define a
  `local_database_path` configuration default (`fadb/nr/nr` under
  `ROTIFER_DATA`), copied from `rotifer.db.local.__init__`. Nothing
  in either package reads this key; UniProt access does not use a
  local FASTA database. Flagged with a `Note` in each module's
  docstring rather than removed.
- `lib/rotifer/db/neighbors.py` carries three further errors found
  while documenting it, all of which would raise at runtime:
  `post_ids` (line ~121) and `runRneighbors` (line ~142) select
  their branch with two consecutive `if` statements where
  `if`/`elif` is meant, so the `asm` filter is discarded whenever
  `gacc` is empty and the `missing` output format also runs the
  neighborhood branch; and `query` (line ~188) calls
  `str(x[0], lost)` and `join(lost)` where `str(x[0])` and
  `'\n'.join(lost)` are meant. Documented as `Notes` on the
  affected methods.
- `lib/rotifer/db/local/core.py` (`FileCollection.__init__`):
  `checksum` defaults to True and appends a `checksum` column name,
  but no checksum value is ever computed (`hashlib` is imported and
  unused), so the row tuples are one element shorter than the
  column list and building the dataframe raises. Documented as a
  `Notes` entry on the class.
- `lib/rotifer/db/uniprot/webapi/idmapping.py`: the functions
  through `get_id_mapping_results_stream` are adapted from
  UniProt's own published ID mapping REST client script. The file
  ends with a large triple-quoted block of commented-out
  exploratory/debug code (duplicate
  `get_data_frame_from_tsv_results`, ad hoc test calls). It is
  inert and was left in place; removing dead code is a cleanup
  task, not a documentation one.

## Site-level notes

- The module-level `@classmethod from_genbank` in
  `lib/rotifer/genome/__init__.py:4` cannot be documented by
  autodoc (`not a callable object`); it is skipped through
  `AUTOSUMMARY_SKIP_MEMBERS` in `docs/conf.py`.
- The Read the Docs PDF build (`formats: [pdf]`) has not been
  exercised locally; if LaTeX chokes on the generated pages, drop
  `pdf` from `.readthedocs.yaml` or fix the offending markup.
- sphinx-iconify's web component is self-hosted and the icons used
  by the site are preloaded in
  `docs/_static/vendor/iconify-preload.js`, so pages render
  offline. An icon used without extending the preload file will
  only render when the viewer's browser can reach
  `api.iconify.design`.
