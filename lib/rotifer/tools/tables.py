__doc__ = """
Helpers to inspect text tables.

``tdesc`` replaces the shell function of the same name that the ROTIFER
Lmod module used to define; it is also a command.
"""

def tdesc(path, sep='\t'):
    """
    Number the columns of a delimited text table.

    Parameters
    ----------
    path : str
        Table file name, or ``-`` for standard input.
    sep : str, default TAB
        Column delimiter (a regular expression).

    Returns
    -------
    list of (int, str)
        Zero-based index and name of each column in the first line.

    Examples
    --------
    >>> from rotifer.tools.tables import tdesc
    >>> tdesc("table.tsv")
    [(0, 'ID'), (1, 'architecture'), ...]
    """
    from rotifer.io import open_input, split_fields
    with open_input(path) as handle:
        header = handle.readline().rstrip('\n')
    return list(enumerate(split_fields(header, sep)))
