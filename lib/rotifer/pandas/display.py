"""Pandas memory and display-option helpers."""

import sys
import pandas as pd

def optimize_memory_usage(df):
    '''
    Optimize memory usage of a Pandas DataFrame.
    
    Usage:
      from rotifer.pandas.functions import optimize_df
      df = optimize_memory_usage(df)

    Returns:
      A memory optimized version of the input Pandas
      DataFrame with the exact same content as the input.

    Parameters:
      df : a Pandas DataFrame
    '''

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
    A function to set pandas display options to print entire
    DataFrames instead of truncating output to a few columns
    and rows.

    Parameters
    ----------
    max_rows : int, default 10000000000
        Maximum number of DataFrame rows to print.
    max_columns : int, default 1000000000
        Maximum number of DataFrame columns to show.
    max_colwidth : int, default 1000000000
        Maximum length of columns. All columns wider than
        this will be truncated.
    width : int, default 100000
        Maximum length for rows.
    verbose : bool, default False
        Warn on STDERR when changing options.

    Warnings
    --------
    Very large Pandas DataFrames could take too long to print!

    Examples
    --------
    >>> import rotifer.pandas.functions as rpf
    >>> rpf.print_everything()
    """

    pd.options.display.max_rows = max_rows
    pd.options.display.max_columns = max_columns
    pd.options.display.max_colwidth = max_colwidth
    pd.options.display.width = width

    if verbose:
        print('Pandas options have been reset!', file=sys.stderr)

def show_display_options():
    '''
    List all Pandas display options that control printed output
    formatting (see pandas.options.display).
    '''
    return pd.DataFrame([ (x,getattr(pd.options.display,x)) for x in dir(pd.options.display) ], columns=('option','value'))
