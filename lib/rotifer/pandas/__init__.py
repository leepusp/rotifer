"""
Text rendering, block reports and display helpers for pandas objects.

The package works on :class:`pandas.DataFrame` and :class:`pandas.Series`
objects and needs pyarrow, which it uses to pad and join columns of text.

Submodules:

:mod:`rotifer.pandas.text`
    Text rendering and padding: :func:`to_string` writes a frame as aligned
    plain text, :func:`pad` pads strings to a common width and
    :func:`to_color` adds terminal colour codes to each character.
:mod:`rotifer.pandas.blocks`
    Block reports: :func:`to_blocks` splits a frame by group and writes each
    group as a text block, using :func:`default_header` and
    :func:`default_block_rows`.
:mod:`rotifer.pandas.yaml_blocks`
    Tables in Markdown and YAML: :func:`md_table`, :func:`read_md_table`,
    :func:`to_yaml` and :func:`read_yaml`.
:mod:`rotifer.pandas.display`
    Memory and display helpers: :func:`optimize_memory_usage`,
    :func:`print_everything` and :func:`show_display_options`.
:mod:`rotifer.pandas.functions`
    Re-exports some of the above for code that imports the old module name.

The most common calls::

    >>> import pandas as pd
    >>> import rotifer.pandas as rp
    >>> df = pd.DataFrame({'gene': ['a', 'bb'], 'length': [10, 200]})
    >>> print(rp.to_string(df, index=False))
    gene : length
    a    : 10    
    bb   : 200   
    >>> rp.pad(df['gene']).tolist()
    [' a', 'bb']

The attribute ``rotifer.pandas.to_network`` is deprecated and raises a
:class:`FutureWarning` when accessed (a module-level ``__getattr__``, PEP
562). Use :meth:`rotifer.genome.data.NeighborhoodDF.to_network` instead.
"""

from .text import to_string, pad, pad_columns, to_color
from .blocks import default_header, default_block_rows, to_blocks
from .display import optimize_memory_usage, print_everything, show_display_options
from .yaml_blocks import md_table, read_md_table, to_yaml, read_yaml

__all__ = [
    'to_string', 'default_header', 'default_block_rows', 'to_blocks',
    'pad', 'to_color', 'optimize_memory_usage', 'print_everything',
    'show_display_options', 'md_table', 'read_md_table', 'to_yaml', 'read_yaml',
]


def _to_network_shim(df, target=['pfam'], ftype=['CDS'], interaction=True, ignore=[], strand=True,
                     separator='+', separator_regex=False, stats=None, replace=None,
                     replace_regex=False):
    """
    Build a network from a neighbourhood table.

    .. deprecated::
       Will be removed on 2026-11-30. Use
       :meth:`rotifer.genome.data.NeighborhoodDF.to_network` instead.

    Returned by the deprecated ``rotifer.pandas.to_network`` attribute. A
    plain frame is wrapped in a :class:`~rotifer.genome.data.NeighborhoodDF`
    first, and the call is passed on to its ``to_network`` method with the
    same arguments (``target`` is the first positional argument there too).
    An empty ``ignore`` is sent as None.

    Parameters
    ----------
    df : pandas.DataFrame
        Neighbourhood table.
    target : list of str, default ['pfam']
        Columns whose values are the nodes of the network.
    ftype : list of str, default ['CDS']
        Feature types to keep.
    interaction : bool, default True
        Passed on to ``NeighborhoodDF.to_network``.
    ignore : list, default []
        Values to leave out of the network.
    strand : bool, default True
        Passed on to ``NeighborhoodDF.to_network``.
    separator : str, default '+'
        Separator of values joined in one cell.
    separator_regex : bool, default False
        Whether ``separator`` is a regular expression.
    stats : optional
        Statistics to compute, passed on unchanged.
    replace : optional
        Replacement of node names, passed on unchanged.
    replace_regex : bool, default False
        Whether ``replace`` holds regular expressions.

    Returns
    -------
    Whatever ``NeighborhoodDF.to_network`` returns.
    """
    from rotifer.genome.data import NeighborhoodDF
    if not isinstance(df, NeighborhoodDF):
        df = NeighborhoodDF(df)
    return df.to_network(target, strand=strand, ftype=ftype, interaction=interaction,
                         ignore=ignore or None, separator=separator,
                         separator_regex=separator_regex, replace=replace,
                         replace_regex=replace_regex, stats=stats)


def __getattr__(name):
    """Serve the deprecated ``to_network`` attribute, with a FutureWarning (PEP 562)."""
    if name == 'to_network':
        import warnings
        warnings.warn(
            "rotifer.pandas.to_network is deprecated and will be removed on 2026-11-30 (issue #160); "
            "use NeighborhoodDF.to_network or rotifer.graph.network.adjacency_edges instead",
            FutureWarning,
            stacklevel=2,
        )
        return _to_network_shim
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
