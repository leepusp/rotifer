"""
Memory and display helpers for pandas objects.

Contains :func:`optimize_memory_usage`, :func:`print_everything` and
:func:`show_display_options`. They are also importable from
:mod:`rotifer.pandas` and :mod:`rotifer.pandas.functions`.
"""

import sys
import pandas as pd

def optimize_memory_usage(df):
    """
    Reduce the memory a DataFrame uses.

    Integer columns are downcast to the smallest integer type that holds
    their values. Float columns are downcast to float32 when every value fits,
    which can change values in the last digits. Text columns in which at most
    half of the values are distinct (missing values count as one value)
    become categorical. Other columns are left as they are.

    Parameters
    ----------
    df : pandas.DataFrame
        Frame to shrink. It is not modified.

    Returns
    -------
    pandas.DataFrame
        A new frame with the same index and column labels and the converted
        columns. An empty frame is returned as a copy.

    Notes
    -----
    The values of integer and text columns are unchanged, but the dtypes are
    not, so code that depends on int64, float64 or object columns may behave
    differently on the result.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas import optimize_memory_usage
    >>> df = pd.DataFrame({'n': [1, 2, 3, 4], 'tag': ['a', 'a', 'a', 'b']})
    >>> small = optimize_memory_usage(df)
    >>> str(small['n'].dtype)
    'int8'
    >>> str(small['tag'].dtype)
    'category'
    >>> str(df['n'].dtype)
    'int64'
    """

    df = df.copy()
    if df.empty:
        return df

    def convert(s):
        if pd.api.types.is_integer_dtype(s.dtype):   # numpy and nullable Int*
            return pd.to_numeric(s, downcast='integer')
        if pd.api.types.is_float_dtype(s.dtype):
            return pd.to_numeric(s, downcast='float')
        # object or pandas' 'str' dtype
        if pd.api.types.is_string_dtype(s.dtype) and s.nunique(dropna=False) / len(s) <= 0.5:
            return s.astype('category')
        return s

    # Positional, so duplicate or non-string column labels work too
    for i in range(df.shape[1]):
        df.isetitem(i, convert(df.iloc[:, i]))
    return df

def print_everything(max_rows=10000000000, max_columns=1000000000, max_colwidth=1000000000, width=100000, verbose=False):
    """
    Set the pandas display options so that whole frames are printed.

    The options ``display.max_rows``, ``display.max_columns``,
    ``display.max_colwidth`` and ``display.width`` are set for the rest of
    the session, so frames are no longer truncated to a few rows and columns.
    Use :func:`pandas.option_context` for a temporary change.

    Parameters
    ----------
    max_rows : int, default 10000000000
        Maximum number of DataFrame rows to print.
    max_columns : int, default 1000000000
        Maximum number of DataFrame columns to show.
    max_colwidth : int, default 1000000000
        Maximum width of a column. Wider columns are
        truncated.
    width : int, default 100000
        Maximum width of a row.
    verbose : bool, default False
        Warn on STDERR when changing options.

    See Also
    --------
    show_display_options : List the current display options.

    Notes
    -----
    This changes global pandas options, which affects all code in the
    process. Printing a very large frame can take a long time.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas import print_everything
    >>> print_everything()
    >>> pd.options.display.max_rows
    10000000000
    >>> pd.reset_option('display.max_rows')
    """

    pd.options.display.max_rows = max_rows
    pd.options.display.max_columns = max_columns
    pd.options.display.max_colwidth = max_colwidth
    pd.options.display.width = width

    if verbose:
        print('Pandas options have been reset!', file=sys.stderr)

def show_display_options():
    """
    List the pandas display options and their current values.

    Returns
    -------
    pandas.DataFrame
        One row per option in ``pandas.options.display``, with columns
        ``option`` and ``value``.

    See Also
    --------
    print_everything : Set the options so that whole frames are printed.
    """
    return pd.DataFrame([ (x,getattr(pd.options.display,x)) for x in dir(pd.options.display) ], columns=('option','value'))
