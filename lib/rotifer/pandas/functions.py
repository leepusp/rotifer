#!/usr/bin/env python3

import os
import sys
import pandas as pd
from collections import OrderedDict
from rotifer.core import functions as rcf

# Defaults
config = rcf.loadConfig(__name__, defaults = {
    'colors': rcf.findDataFiles(":colors/aminoacids.yml"),
})

def pad(data, side='left', fillchar=' '):
    data = data.astype(str)
    if isinstance(data,pd.DataFrame):
        side = { x: side[x] if x in side else 'left' for x in data.columns }
        data = data.apply(lambda x: x.str.pad(x.str.len().max(), side=side[x.name], fillchar=fillchar))
    elif isinstance(data,pd.Series):
        data = data.str.pad(data.str.len().max(), side=side, fillchar=fillchar)
    return data

def to_color(data, colors=config['colors'], style='fg', padding=None, fillchar=' '):
    """
    Convert a DataFrame or Series to DType str and add
    ASCII color codes to each character.

    Usage
    -----
    >>> from rotifer.core import functions as rcf
    >>> colored = rcf.to_color(pd.Series("ACGT"))

    Parameters
    ----------
    colors: str, os.PathLike or dict, default None
      A YAML file or a dictionary, mapping characters to
      color codes.
    style: str, default 'fg'
      Whether apply colors to the background ('bg') or toi
      characters ('fg').
    padding: {'left', 'right', 'both'}, default None
      If set, data will be padded with whitespaces to make
      sure all columns have the same length.
    fillchar: str, default ' '
      Additional character for filling, default is whitespace.
    """
    import yaml
    from rotifer.core import functions as rcf

    # Load colors
    default_color = "000"
    if colors == None or colors == False:
        colors = config["colors"]
    if not isinstance(colors,dict):
        if isinstance(colors,str):
            if os.path.exists(colors):
                colors = yaml.load(open(colors), Loader=yaml.Loader)
            else:
                default_color = colors
                colors = dict()

    def color_res(s, cs):
        if s in colors:
            return cs(s, colors[s])
        else:
            return cs(s, default_color)

    def color_bg(s, color):
        if color == default_color:
            color = f'49;5;{default_color}'
        else:
            color = f'48;5;{color}'
        return f'\033[{color}m{s}\033[0m'

    def color_fg(s, color):
        if color == default_color:
            color = f'39;5;{default_color}'
        else:
            color = f'38;5;{color}'
        return f'\033[{color}m{s}\033[m'

    # Forcing all data to be of the same length
    if padding:
        data = pad(data, side=padding, fillchar=fillchar)

    # Choose the coloring function
    color_switch = {'background':color_bg, 'bg':color_bg, 'foreground':color_fg, 'fg':color_fg}

    # Color the data
    if isinstance(data,pd.DataFrame):
        data = data.applymap(lambda x: ''.join([color_res(y, color_switch[style]) for y in x]))
    elif isinstance(data,pd.Series):
        data = data.map(lambda x: ''.join([color_res(y, color_switch[style]) for y in x]))

    return data

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

    # Make a copy of the input DataFrame
    df = df.copy()

    # int optimization
    converted_int = df.select_dtypes(include=['int'])
    converted_int = converted_int.apply(pd.to_numeric, downcast = "integer")

    # Float
    converted_float = df.select_dtypes(include=['float'])
    converted_float = converted_float.apply(pd.to_numeric,downcast='float')

    # Object
    df_obj = df.select_dtypes(include=['object'])
    converted_obj = pd.DataFrame()

    for col in df_obj.columns:
        num_unique_values = len(df_obj[col].unique())
        num_total_values = len(df_obj[col])
        if num_unique_values / num_total_values <= 0.5:
            converted_obj.loc[:,col] = df_obj[col].astype('category')
        else:
            converted_obj.loc[:,col] = df_obj[col]

    # Adding converted objects to a better df
    df[converted_int.columns] = converted_int
    df[converted_float.columns] = converted_float
    df[converted_obj.columns] = converted_obj
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
