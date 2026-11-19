__doc__ = """
Parser for DaliLite structure files (``.dat``).

DaliLite stores each protein chain of its database as a ``.dat`` file.
``parse_dat`` extracts the chain identifier, the compound description and
the amino acid sequence, which the ``dat2fa`` command writes as FASTA
(used by ``sbin/update_pdb`` to build the DaliLite sequence database).

This module is separate from ``rotifer.io.dali`` (DALI server results) to
stay free of heavy imports: ``dat2fa`` runs once per batch of files.
"""

import re

_ID = re.compile(r'^>>>> *(\S+)')
_SEQUENCE = re.compile(r'^-sequence *"(\S+)$')
_COMPOUND = re.compile(r'^-compnd *"(.+)"$')

def parse_dat(handle):
    """
    Read the chains of DaliLite ``.dat`` files.

    Parameters
    ----------
    handle : iterable of str
        Lines of one or more concatenated ``.dat`` files.

    Yields
    ------
    tuple of (str, str, str)
        ``(identifier, description, sequence)`` for each chain: the id on
        the ``>>>>`` line, the ``-compnd`` text and the ``-sequence``
        residues.

    Notes
    -----
    A record is emitted at each ``-compnd`` line, with the last identifier
    and sequence seen before it. A chain whose file has no ``-sequence``
    line before ``-compnd`` gets the previous chain's sequence (None for the
    first chain), as the former Perl ``dat2fa`` did.

    Examples
    --------
    >>> from rotifer.io.dalidat import parse_dat
    >>> with open("test/data/golden/inputs/dali/a0rrA.dat") as handle:
    ...     identifier, description, sequence = next(parse_dat(handle))
    >>> identifier, sequence[:10]
    ('a0rrA', 'MSSSLYLVFE')
    """
    identifier = sequence = None
    for line in handle:
        line = line.rstrip('\n')
        match = _ID.match(line)
        if match:
            identifier = match.group(1)
            continue
        match = _SEQUENCE.match(line)
        if match:
            sequence = match.group(1)
            continue
        match = _COMPOUND.match(line)
        if match:
            yield identifier, match.group(1), sequence

def dat2fasta(handle):
    """
    Write the chains of DaliLite ``.dat`` lines as FASTA records.

    Parameters
    ----------
    handle : iterable of str
        Lines of one or more concatenated ``.dat`` files.

    Yields
    ------
    str
        One record per chain, ``">id description\nSEQUENCE\n"`` (see
        :func:`parse_dat`).

    Examples
    --------
    >>> from rotifer.io.dalidat import dat2fasta
    >>> with open("test/data/golden/inputs/dali/a0rrA.dat") as handle:
    ...     next(dat2fasta(handle)).split()[0]
    '>a0rrA'
    """
    for identifier, description, sequence in parse_dat(handle):
        yield f">{identifier} {description}\n{sequence}\n"
