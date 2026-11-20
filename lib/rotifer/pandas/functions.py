"""
Names kept for code that imports ``rotifer.pandas.functions``.

This module re-exports :func:`pad`, :func:`to_color` and ``config`` from
:mod:`rotifer.pandas.text`, and :func:`optimize_memory_usage`,
:func:`print_everything` and :func:`show_display_options` from
:mod:`rotifer.pandas.display`. New code can import them from
:mod:`rotifer.pandas` or from those modules. ``config`` is the configuration
loaded under the key ``rotifer.pandas.functions``.
"""

from .text import pad, to_color, config
from .display import optimize_memory_usage, print_everything, show_display_options
