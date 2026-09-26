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
