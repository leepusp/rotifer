"""Pandas helpers: text rendering, block output, YAML blocks and display options."""

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
    from rotifer.genome.data import NeighborhoodDF
    if not isinstance(df, NeighborhoodDF):
        df = NeighborhoodDF(df)
    return df.to_network(target, strand=strand, ftype=ftype, interaction=interaction,
                         ignore=ignore or None, separator=separator,
                         separator_regex=separator_regex, replace=replace,
                         replace_regex=replace_regex, stats=stats)


def __getattr__(name):
    if name == 'to_network':
        import warnings
        warnings.warn(
            "rotifer.pandas.to_network is deprecated and will be removed on 2026-11-30 (issue #TBD); "
            "use NeighborhoodDF.to_network or rotifer.graph.network.adjacency_edges instead",
            FutureWarning,
            stacklevel=2,
        )
        return _to_network_shim
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
