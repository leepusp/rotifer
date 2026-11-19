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
    Yield ``(identifier, description, sequence)`` for each chain in the
    lines of one or more concatenated DaliLite ``.dat`` files.

    A record is emitted at each ``-compnd`` line, using the last
    identifier and sequence seen before it.
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
    Yield FASTA records (strings) for the chains in DaliLite ``.dat`` lines.
    """
    for identifier, description, sequence in parse_dat(handle):
        yield f">{identifier} {description}\n{sequence}\n"
