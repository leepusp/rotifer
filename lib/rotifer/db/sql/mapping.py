"""
The SQL of an identifier mapping, without a dialect.

Translating identifiers through UniProt is one query whichever engine
answers it, and the parts that differ between engines are small and
local: how a value is bound into a statement, and what to do when
there are too many values to write into one. Everything else -- which
branches a query needs, how the two ends are joined, which shapes need
no join at all -- follows from what was asked and is the same
everywhere.

That shared part lives here, so a backend supplies a table name and a
way of binding values and gets the statement back.

The shape of the query
----------------------
A mapping matches identifiers in one set of databases, finds the
UniProtKB accessions they belong to, and returns the identifiers those
accessions carry in another set. Both ends accept
:data:`UNIPROTKB`, which is not a row of the table but its join key,
and that is what decides the shape:

``source`` is only the accession
    The accessions are the key, so every row of a matching accession
    already carries the answer: no join.
``target`` is only the accession
    The rows matching the queried identifiers already carry the
    accession: no join either.
anything else
    Two ends to satisfy, so the source side is built first and the
    target side is restricted to the accessions it found. Restricted,
    not merely joined: a join whose right hand side is the whole table
    builds a hash of every row in it, which on a few billion rows does
    not return.

See Also
--------
rotifer.db.uniprot.clickhouse : this SQL, bound for ClickHouse
rotifer.db.uniprot.sqlite3 : the same, bound for SQLite3
"""

import types

import rotifer
logger = rotifer.logging.getLogger(__name__)

#: Name standing for a UniProtKB accession where a database is
#: expected. Repeated from :class:`rotifer.db.methods.MappingCursor` so
#: this module can be read on its own.
UNIPROTKB = 'UniProtKB-AC'

#: Columns of the mapping table these statements read.
TABLE_COLUMNS = ('accession', 'id_type', 'id')


class Binder:
    """
    Collect the values a statement binds, and render placeholders.

    Dialects differ in how a value is referred to from SQL and in what
    they accept for a list of them. A binder hides both, so the SQL
    above can be written once: it asks for a placeholder and gets back
    whatever this dialect writes there, while the value is kept for
    the driver.

    This class is abstract; :class:`NamedBinder` and
    :class:`QmarkBinder` are the two spellings in use.
    """

    def scalar(self, value):
        """
        Bind one value and return the placeholder standing for it.

        Parameters
        ----------
        value : object

        Returns
        -------
        str
        """
        raise NotImplementedError

    def collection(self, values):
        """
        Bind several values and return the ``IN`` list standing for them.

        Parameters
        ----------
        values : iterable

        Returns
        -------
        str
            Including its parentheses, so it follows ``IN`` directly.
        """
        raise NotImplementedError

    @property
    def parameters(self):
        """
        Whatever the driver should be handed alongside the statement.

        Returns
        -------
        dict or list
        """
        raise NotImplementedError


class NamedBinder(Binder):
    """
    Bind values by name, as ``%(name)s``, with a dictionary.

    This is what ClickHouse's driver expects, and it takes a list as
    one bound tuple rather than as a placeholder each.
    """

    def __init__(self, prefix='p'):
        self._prefix = prefix
        self._values = {}
        self._next = 0

    def _name(self):
        name = f'{self._prefix}{self._next}'
        self._next += 1
        return name

    def scalar(self, value):
        name = self._name()
        self._values[name] = value
        return f'%({name})s'

    def collection(self, values):
        name = self._name()
        self._values[name] = tuple(values)
        return f'%({name})s'

    @property
    def parameters(self):
        return dict(self._values)


class QmarkBinder(Binder):
    """
    Bind values by position, as ``?``, with a list.

    This is what the standard library's SQLite3 driver expects. A list
    becomes one placeholder per value, so a caller must keep batches
    within the engine's limit on how many a statement may carry.
    """

    def __init__(self):
        self._values = []

    def scalar(self, value):
        self._values.append(value)
        return '?'

    def collection(self, values):
        values = list(values)
        self._values.extend(values)
        return '(' + ", ".join(['?'] * len(values)) + ')'

    @property
    def parameters(self):
        return list(self._values)


def _databases(names, open_includes_pivot=True):
    """
    Split a list of databases into the pivot and the rest.

    The two ends read an omitted list differently, and deliberately.
    An omitted ``source`` means an identifier may be anything, the
    accession included, since the caller has not said what they are
    holding. An omitted ``target`` means every database the table
    cross-references, which does not include the accession: that is
    the table's key rather than one of its rows, and returning it
    would answer a question nobody asked.

    Parameters
    ----------
    names : list of str or None
        Databases asked for, or None for every one of them.
    open_includes_pivot : bool, default True
        Whether an omitted list covers the accession itself. True for
        the source end, False for the target.

    Returns
    -------
    tuple of (bool, list or None)
        Whether the accession itself was asked for, and the real
        databases named, which is None when every one was wanted.
    """
    if isinstance(names, types.NoneType):
        return open_includes_pivot, None
    return UNIPROTKB in names, [ x for x in names if x != UNIPROTKB ]


def source_subquery(table, restrict, source, binder, release=None):
    """
    Build the subquery naming the identifiers a query asked about.

    Two kinds of identifier can be asked for and they live in
    different columns, so each becomes its own branch: UniProtKB
    accessions are the table's join key, everything else is a row of
    it. Both branches are restricted to the identifiers given, which
    is what keeps the query from reading the table.

    Parameters
    ----------
    table : str
        Qualified name of the mapping table.
    restrict : callable
        Called with a column name, returns the condition selecting the
        queried identifiers in that column. It is called once per
        place the condition appears, and binds its values each time:
        a dialect binding by position needs its values in the order
        the statement mentions them, so a condition cannot be a string
        written once and repeated.
    source : list of str or None
        Databases the identifiers belong to, or None for any.
    binder : Binder
        Collects the values bound along the way.
    release : str, optional
        Restrict to one release of the data.

    Returns
    -------
    str
        A subquery yielding ``value``, ``db`` and ``accession``.
    """
    wants_accession, others = _databases(source)
    branches = []

    if wants_accession:
        where = [restrict('accession')]
        if release:
            where.append(f'f.release = {binder.scalar(release)}')
        branches.append(
            f"SELECT accession AS value, '{UNIPROTKB}' AS db, accession"
            f' FROM {table} AS f WHERE {" AND ".join(where)}'
        )
    if isinstance(others, types.NoneType) or others:
        where = [restrict('id')]
        if others:
            where.append(f'f.id_type IN {binder.collection(others)}')
        if release:
            where.append(f'f.release = {binder.scalar(release)}')
        branches.append(
            f'SELECT id AS value, id_type AS db, accession'
            f' FROM {table} AS f WHERE {" AND ".join(where)}'
        )
    return " UNION ALL ".join(branches)


def mapping_query(table, restrict, source, target, binder, release=None):
    """
    Build the whole statement, both ends included.

    Every fragment is built in the order the finished statement
    mentions it, because a dialect binding by position offers no other
    way of matching values to placeholders.

    Parameters
    ----------
    table : str
        Qualified name of the mapping table.
    restrict : callable
        Called with a column name, returns the condition selecting the
        queried identifiers. See :func:`source_subquery`.
    source, target : list of str or None
        The two ends of the mapping. None means every database.
    binder : Binder
        Collects the values bound along the way.
    release : str, optional
        Restrict to one release of the data.

    Returns
    -------
    str
        A SELECT returning ``source``, ``source_type``, ``accession``,
        ``target`` and ``target_type``.
    """
    source_pivot, source_others = _databases(source)
    target_pivot, target_others = _databases(target, open_includes_pivot=False)
    only_source_pivot = source == [UNIPROTKB]
    only_target_pivot = target == [UNIPROTKB]

    # Both ends pinned to a single kind of identifier need no join: one
    # end is then the table's own key, already carried by every row the
    # other end matches. These are the two commonest questions, and
    # joining for them costs an order of magnitude more than reading
    # the rows once.
    if only_source_pivot and not only_target_pivot:
        selects = []
        where = [restrict('accession')]
        if target_others:
            where.append(f'f.id_type IN {binder.collection(target_others)}')
        if release:
            where.append(f'f.release = {binder.scalar(release)}')
        selects.append(
            f"SELECT DISTINCT accession AS source, '{UNIPROTKB}' AS source_type,"
            f' accession AS accession, id AS target, id_type AS target_type'
            f' FROM {table} AS f WHERE {" AND ".join(where)}'
        )
        if target_pivot:
            where = [restrict('accession')]
            if release:
                where.append(f'f.release = {binder.scalar(release)}')
            selects.append(
                f"SELECT DISTINCT accession AS source, '{UNIPROTKB}' AS source_type,"
                f" accession AS accession, accession AS target, '{UNIPROTKB}' AS target_type"
                f' FROM {table} AS f WHERE {" AND ".join(where)}'
            )
        return " UNION ALL ".join(selects)

    if only_target_pivot:
        selects = []
        if isinstance(source_others, types.NoneType) or source_others:
            where = [restrict('id')]
            if source_others:
                where.append(f'f.id_type IN {binder.collection(source_others)}')
            if release:
                where.append(f'f.release = {binder.scalar(release)}')
            selects.append(
                f'SELECT DISTINCT id AS source, id_type AS source_type, accession AS accession,'
                f" accession AS target, '{UNIPROTKB}' AS target_type"
                f' FROM {table} AS f WHERE {" AND ".join(where)}'
            )
        if source_pivot:
            where = [restrict('accession')]
            if release:
                where.append(f'f.release = {binder.scalar(release)}')
            selects.append(
                f"SELECT DISTINCT accession AS source, '{UNIPROTKB}' AS source_type,"
                f" accession AS accession, accession AS target, '{UNIPROTKB}' AS target_type"
                f' FROM {table} AS f WHERE {" AND ".join(where)}'
            )
        return " UNION ALL ".join(selects)

    # Both ends are real. The source side is built first because it
    # comes first in the statement, and the target side is restricted
    # to the accessions it found -- restricted, not merely joined: a
    # join whose right hand side is the whole table hashes every row
    # in it, which on a few billion rows does not return.
    selects = []
    if isinstance(target_others, types.NoneType) or target_others:
        outer = source_subquery(table, restrict, source, binder, release=release)
        inner_source = source_subquery(table, restrict, source, binder, release=release)
        where = [f'accession IN (SELECT accession FROM ({inner_source}))']
        if target_others:
            where.append(f'id_type IN {binder.collection(target_others)}')
        if release:
            where.append(f'release = {binder.scalar(release)}')
        inner = f'SELECT accession, id, id_type FROM {table} WHERE {" AND ".join(where)}'
        selects.append(
            f'SELECT DISTINCT s.value AS source, s.db AS source_type, s.accession AS accession,'
            f' t.id AS target, t.id_type AS target_type'
            f' FROM ({outer}) AS s'
            f' INNER JOIN ({inner}) AS t ON t.accession = s.accession'
        )
    if target_pivot:
        outer = source_subquery(table, restrict, source, binder, release=release)
        selects.append(
            f'SELECT DISTINCT s.value AS source, s.db AS source_type, s.accession AS accession,'
            f" s.accession AS target, '{UNIPROTKB}' AS target_type"
            f' FROM ({outer}) AS s'
        )
    return " UNION ALL ".join(selects)
