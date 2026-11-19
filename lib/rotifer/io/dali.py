"""Read a DALI server TSV export as a star alignment.

Public function
---------------
dali_tsv
    Build the alignment, with optional DSSP rows, from one TSV file.

How the private functions fit together
--------------------------------------
``dali_tsv`` calls ``_read_dali_tsv`` to get the hit table, the query
sequence and the query DSSP string. ``_row_blocks`` turns the ``qstarts``,
``sstarts`` and ``lengths`` columns of a row into aligned blocks
(``_parse_int_list`` reads each column). ``_collect_insert_sizes`` goes over
all blocks of all hits and finds how many insertion columns each query
position needs. ``_make_star_query`` writes the query (or its DSSP string)
with those columns left as gaps. ``_make_star_row`` writes one hit, or its
DSSP string (decompressed by ``_decompress_dssp``), in the same columns, and
formats inserted residues with ``_format_unaligned``.
``_get_subject_description`` picks the description of a hit.
"""

import re
from pathlib import Path

import pandas as pd


def _parse_int_list(x):
    """
    Return every run of digits in ``x`` as an integer.

    Parameters
    ----------
    x : object
        A cell such as ``"1,5,9"`` or ``"[1 5 9]"``; it is converted with
        ``str`` first, so a single number also works.

    Returns
    -------
    list of int
        The integers in the order they appear (empty if there are none).
        Signs are ignored.
    """
    return [int(n) for n in re.findall(r"\d+", str(x))]


def _decompress_dssp(x):
    """
    Expand a run-length encoded DSSP string.

    Parameters
    ----------
    x : str or None
        The ``compressed-sbjct-dssp`` cell, pairs of a count and a letter
        such as ``"3H2-4E"``. A missing value (NaN or None) is accepted.

    Returns
    -------
    str
        The expanded string (``"HHH--EEEE"`` for the example), or ``""`` when
        ``x`` is missing.
    """
    if pd.isna(x):
        return ""

    parts = re.findall(r"(\d+)([A-Za-z-])", str(x))
    return "".join(char * int(n) for n, char in parts)


def _format_unaligned(text, mode="lower"):
    """
    Write residues that fall outside the aligned blocks as ``mode`` asks.

    Parameters
    ----------
    text : str
        The residues of one insertion.
    mode : {'lower', 'keep', 'gap', 'remove'}, default 'lower'
        ``keep`` returns ``text`` as it is, ``lower`` lowercases it, ``gap``
        replaces each residue with ``-`` and ``remove`` returns an empty
        string.

    Returns
    -------
    str
        The formatted text.

    Raises
    ------
    ValueError
        If ``mode`` is not one of the four modes.
    """
    if mode == "keep":
        return text
    if mode == "lower":
        return text.lower()
    if mode == "gap":
        return "-" * len(text)
    if mode == "remove":
        return ""

    raise ValueError("unaligned must be one of: 'keep', 'lower', 'gap', 'remove'")


def _read_dali_tsv(path, max_hits=None):
    """
    Read the hit table and the query's sequence and DSSP from a DALI TSV.

    Parameters
    ----------
    path : str or os.PathLike
        The TSV file. Lines starting with ``#`` are comments; the query's
        data is in the ``# SEQUENCE:`` and ``# DSSP:`` ones.
    max_hits : int, optional
        Keep only the first ``max_hits`` rows of the table (all if None).

    Returns
    -------
    dali : pandas.DataFrame
        One row per hit.
    query_seq : str
        The query sequence.
    query_dssp : str or None
        The query's DSSP string, None if the file has no ``# DSSP:`` line.

    Raises
    ------
    ValueError
        If the file has no ``# SEQUENCE:`` line.
    """
    path = Path(path)

    query_seq = None
    query_dssp = None

    with path.open() as f:
        for line in f:
            if line.startswith("# SEQUENCE:"):
                query_seq = line.split(":", 1)[1].strip()
            elif line.startswith("# DSSP:"):
                query_dssp = line.split(":", 1)[1].strip()

    if query_seq is None:
        raise ValueError("Could not find '# SEQUENCE:' line in DALI file.")

    dali = pd.read_csv(path, sep="\t", comment="#")

    if max_hits is not None:
        dali = dali.head(max_hits).copy()

    return dali, query_seq, query_dssp


def _row_blocks(row):
    """
    Return the aligned blocks of one hit as 1-based coordinates.

    The code treats the ``qstarts`` and ``sstarts`` of the DALI TSV as
    0-based and adds 1 to each. (This has not been checked against the DALI
    server's own documentation.)

    Parameters
    ----------
    row : pandas.Series
        One row of the DALI table, with the columns ``qstarts``,
        ``sstarts`` and ``lengths`` (lists of integers, see
        :func:`_parse_int_list`).

    Returns
    -------
    list of tuple of int
        One ``(q_start, s_start, length)`` per block: the 1-based start in
        the query, the 1-based start in the hit and the block length.
    """
    qstarts = [x + 1 for x in _parse_int_list(row["qstarts"])]
    sstarts = [x + 1 for x in _parse_int_list(row["sstarts"])]
    lengths = _parse_int_list(row["lengths"])

    return list(zip(qstarts, sstarts, lengths))


def _collect_insert_sizes(dali, query_len):
    """
    Find how many insertion columns each query position needs.

    For every hit, the residues between one block and the previous one (or
    the start of the hit) go before the first query position of the block;
    the residues after the last block go after the last query position. The
    widest insertion over all hits sets the width of each position.

    Parameters
    ----------
    dali : pandas.DataFrame
        The hit table, with the columns ``qstarts``, ``sstarts``,
        ``lengths`` and ``sbjct-sequence``.
    query_len : int
        Length of the query sequence.

    Returns
    -------
    list of int
        ``query_len + 1`` widths: the columns before each query position,
        then the columns after the last one.
    """
    insert_sizes = [0] * (query_len + 1)

    for _, row in dali.iterrows():
        subject_seq = row["sbjct-sequence"]
        prev_s_end = 1

        for q_start, s_start, length in _row_blocks(row):
            q_insert_pos = q_start - 1
            subject_insertion_len = s_start - prev_s_end

            insert_sizes[q_insert_pos] = max(
                insert_sizes[q_insert_pos],
                subject_insertion_len,
            )

            prev_s_end = s_start + length

        tail_len = len(subject_seq) - prev_s_end + 1
        insert_sizes[query_len] = max(
            insert_sizes[query_len],
            max(0, tail_len),
        )

    return insert_sizes


def _make_star_query(query_seq, insert_sizes):
    """
    Write the query (or a string as long as it) on the star alignment.

    Parameters
    ----------
    query_seq : str
        The query sequence or its DSSP string.
    insert_sizes : list of int
        ``len(query_seq) + 1`` widths, see :func:`_collect_insert_sizes`.

    Returns
    -------
    str
        The residues uppercased, with ``insert_sizes[i]`` gaps before
        position ``i`` and ``insert_sizes[-1]`` gaps at the end.
    """
    aln = []

    for i, char in enumerate(query_seq):
        aln.append("-" * insert_sizes[i])
        aln.append(char.upper())

    aln.append("-" * insert_sizes[len(query_seq)])

    return "".join(aln)


def _make_star_row(row, residues, query_len, insert_sizes, unaligned="lower"):
    """
    Lay one hit's residues (or its DSSP string) on the star alignment.

    Residues the DALI blocks align to the query go in the query's columns,
    uppercase; residues between blocks, before the first block and after
    the last one go in the insertion columns that precede each query
    position (and in the tail after the last one), formatted by
    ``unaligned``, and are padded with ``-`` to the width that
    :func:`_collect_insert_sizes` reserved for that position.

    Parameters
    ----------
    row : pandas.Series
        One row of the DALI table; only ``qstarts``, ``sstarts`` and
        ``lengths`` are read (0-based starts, see :func:`_row_blocks`).
    residues : str
        The hit's sequence (``sbjct-sequence``) or its decompressed DSSP
        string; both are indexed by the hit's residue numbers.
    query_len : int
        Length of the query sequence.
    insert_sizes : list of int
        ``query_len + 1`` widths: the insertion columns before each query
        position, then after the last one. All zeros with
        ``unaligned='remove'``.
    unaligned : {'lower', 'keep', 'gap', 'remove'}, default 'lower'
        How residues outside the aligned blocks are written (see
        :func:`_format_unaligned`).

    Returns
    -------
    str
        The aligned row, ``query_len + sum(insert_sizes)`` characters long.
        Query positions no block covers are ``-``.

    Raises
    ------
    ValueError
        If ``unaligned`` is not one of the four modes.
    """
    on_query = ["-"] * query_len
    inserts = [""] * (query_len + 1)

    prev_s_end = 1

    for q_start, s_start, length in _row_blocks(row):
        q_insert_pos = q_start - 1

        insertion = residues[prev_s_end - 1:s_start - 1]
        inserts[q_insert_pos] += _format_unaligned(
            insertion,
            mode=unaligned,
        )

        for offset in range(length):
            q_idx = q_start - 1 + offset
            s_idx = s_start - 1 + offset

            if 0 <= q_idx < query_len and 0 <= s_idx < len(residues):
                on_query[q_idx] = residues[s_idx].upper()

        prev_s_end = s_start + length

    tail = residues[prev_s_end - 1:]
    inserts[query_len] += _format_unaligned(tail, mode=unaligned)

    aln = []

    for i, char in enumerate(on_query):
        insertion = inserts[i]
        aln.append(insertion + "-" * (insert_sizes[i] - len(insertion)))
        aln.append(char)

    tail = inserts[query_len]
    aln.append(tail + "-" * (insert_sizes[query_len] - len(tail)))

    return "".join(aln)


def _get_subject_description(row):
    """
    Pick the description of a hit.

    Parameters
    ----------
    row : pandas.Series
        One row of the DALI table.

    Returns
    -------
    str
        ``description`` if the row has it and it is not missing, else
        ``sbjct-description``, else the hit's id (``sbjct``).
    """
    if "description" in row.index and pd.notna(row["description"]):
        return row["description"]

    if "sbjct-description" in row.index and pd.notna(row["sbjct-description"]):
        return row["sbjct-description"]

    return row["sbjct"]


def dali_tsv(
    path,
    max_hits=None,
    n_structures=None,
    add_dssp=True,
    unaligned="lower",
):
    """
    Build a star alignment, with DSSP rows, from a DALI server TSV file.

    The DALI server's TSV export has one row per structure aligned to the
    query (columns ``query``, ``sbjct``, ``qstarts``, ``sstarts``,
    ``lengths``, ``sbjct-sequence``, ``compressed-sbjct-dssp`` and,
    optionally, ``description``/``sbjct-description`` and
    ``query-description``), and the query's sequence and secondary
    structure in ``# SEQUENCE:`` and ``# DSSP:`` comment lines. Every hit is
    laid on the query's columns; residues a hit has between its aligned
    blocks open insertion columns in all rows, so every row has the same
    length.

    Parameters
    ----------
    path : str or os.PathLike
        The TSV file.
    max_hits : int, optional
        Keep only the first ``max_hits`` rows of the table (all if None).
    n_structures : int, optional
        Write DSSP rows only for the first ``n_structures`` hits (all if
        None). Sequence rows are written for every kept hit.
    add_dssp : bool, default True
        Write DSSP rows: the query's (when the file has a ``# DSSP:`` line)
        and the hits'.
    unaligned : {'lower', 'keep', 'gap', 'remove'}, default 'lower'
        How residues outside the aligned blocks are written: lowercase, as
        they are, as gaps, or not at all (with ``remove`` the alignment has
        exactly the query's columns).

    Returns
    -------
    rotifer.devel.beta.sequence.sequence
        An alignment whose ``df`` has the columns ``id``, ``sequence``,
        ``type`` and ``description``. DSSP rows come first, with ids
        ``ss_from:<id>`` and type ``residue_annotation``: the query's, then
        the hits'. Sequence rows (type ``sequence``) follow: the query,
        then the hits in file order. A hit's description is ``description``,
        else ``sbjct-description``, else its id.

    Raises
    ------
    ValueError
        If the file has no ``# SEQUENCE:`` line, if no hit is left after
        ``max_hits``, or if ``unaligned`` is not one of the four modes.
    ImportError
        If ``rotifer.devel.beta.sequence`` cannot be imported; it is loaded
        when this function runs, not when the module is imported.

    Examples
    --------
    >>> from rotifer.io.dali import dali_tsv
    >>> aln = dali_tsv("test/io/data/toy_dali.tsv", add_dssp=False)
    >>> print(aln.df[["id", "sequence"]].to_string(index=False))
    id        sequence
    q1 ACDEF-----GHIK-
    s1 ACD--xxxxxGHIKz
    s2 --DEF-----GH---
    """
    from rotifer.devel.beta.sequence import sequence

    dali, query_seq, query_dssp = _read_dali_tsv(path, max_hits=max_hits)

    if dali.empty:
        raise ValueError("No DALI hits available after filtering.")

    # sequences included in the alignment
    dali_seq = dali.copy()

    # structures with DSSP annotation rows
    if n_structures is not None:
        dali_dssp = dali.head(n_structures).copy()
    else:
        dali_dssp = dali.copy()

    if unaligned == "remove":
        insert_sizes = [0] * (len(query_seq) + 1)
    else:
        insert_sizes = _collect_insert_sizes(dali_seq, len(query_seq))

    query_id = dali_seq["query"].iloc[0]

    query_description = (
        dali_seq["query-description"].iloc[0]
        if "query-description" in dali_seq.columns
        else query_id
    )

    query_aln = _make_star_query(query_seq, insert_sizes)

    annotation_rows = []
    sequence_rows = []

    if add_dssp and query_dssp is not None:
        annotation_rows.append({
            "id": f"ss_from:{query_id}",
            "sequence": _make_star_query(query_dssp, insert_sizes),
            "type": "residue_annotation",
            "description": query_description,
        })

    sequence_rows.append({
        "id": query_id,
        "sequence": query_aln,
        "type": "sequence",
        "description": query_description,
    })

    # output all sequence rows from max_hits
    for _, row in dali_seq.iterrows():
        subject_id = row["sbjct"]
        subject_description = _get_subject_description(row)

        subject_aln = _make_star_row(
            row,
            row["sbjct-sequence"],
            len(query_seq),
            insert_sizes,
            unaligned=unaligned,
        )

        sequence_rows.append({
            "id": subject_id,
            "sequence": subject_aln,
            "type": "sequence",
            "description": subject_description,
        })

    # output DSSP only for n_structures
    if add_dssp:
        for _, row in dali_dssp.iterrows():
            subject_id = row["sbjct"]
            subject_description = _get_subject_description(row)

            subject_dssp = _decompress_dssp(row["compressed-sbjct-dssp"])

            annotation_rows.append({
                "id": f"ss_from:{subject_id}",
                "sequence": _make_star_row(
                    row,
                    subject_dssp,
                    len(query_seq),
                    insert_sizes,
                    unaligned=unaligned,
                ),
                "type": "residue_annotation",
                "description": subject_description,
            })

    df = pd.DataFrame(annotation_rows + sequence_rows)

    aln = sequence()
    aln.df = df
    aln._reset()

    return aln