"""Plain-text rendering of DataFrames: to_string, pad and to_color."""

import os
import pandas as pd
from rotifer.core import functions as rcf

# Defaults
config = rcf.loadConfig('rotifer.pandas.functions', defaults = {
    'colors': rcf.findDataFiles(":colors/aminoacids.yml"),
})

# The name above is the config key users know (a user file named after
# rotifer.pandas.functions), so it must not follow this module's name.

def pad_columns(arrays, names=None, align='left', fillchar=' '):
    """Pad pyarrow string arrays to one common width per column.

    align is where the text sits: 'left' pads on the right, 'right' on the left,
    'center' on both sides (extra fill goes where str.center puts it). The width
    of a column is its longest string, at least len(name) when names are given;
    an empty or all-null column has width 0. Nulls are padded as empty strings.
    Returns (padded arrays, widths)."""
    import pyarrow as pa
    import pyarrow.compute as pc

    if align not in ('left', 'right', 'center'):
        raise ValueError(f"align must be 'left', 'right' or 'center', not {align!r}")
    padded, widths = [], []
    for k, arr in enumerate(arrays):
        arr = pc.fill_null(arr, '')
        width = pc.max(pc.utf8_length(arr)).as_py() or 0
        if names is not None:
            width = max(width, len(names[k]))
        if align == 'center': # arrow's utf8_center splits odd padding the other way round
            arr = pa.array([s.center(width, fillchar) for s in arr.to_pylist()], type=pa.string())
        else:
            arr = (pc.utf8_lpad if align == 'right' else pc.utf8_rpad)(arr, width, padding=fillchar)
        padded.append(arr)
        widths.append(width)
    return padded, widths

def to_string(df, sep=" : ", index=True, na_rep=" ", insert=None,
                      float_format=None, align="left", **kwargs):
    import pyarrow as pa
    import pyarrow.compute as pc

    just = str.rjust if align == "right" else str.ljust

    if insert:
        df = df.copy()
        for data in insert:
            if len(data) != 3:
                raise ValueError(f"insert entries must be (loc, column, value) tuples, got {data!r}")
            df.insert(loc=data[0], column=data[1], value=data[2])

    columns = []
    if index:
        columns.append((str(df.index.name or ""), df.index.to_series().astype("string")))
    for c in df.columns:
        s = df[c]
        if float_format is not None and pd.api.types.is_float_dtype(s.dtype):
            body = s.map(lambda v: na_rep if pd.isna(v) else float_format % v).astype("string")
        else:
            body = s.astype("string")
        columns.append((str(c), body.fillna(na_rep)))

    names = [name for name, _ in columns]
    arrays = [pa.Array.from_pandas(body, type=pa.string()) for _, body in columns]
    padded, widths = pad_columns(arrays, names, align=align)
    header = sep.join(just(name, w) for name, w in zip(names, widths))
    if not padded:
        return header
    lines = pc.binary_join_element_wise(*padded, pa.scalar(sep))
    return header + "\n" + "\n".join(lines.to_pylist())

# pad's side is where the fill goes; pad_columns' align is where the text sits
_PAD_ALIGN = {'left': 'right', 'right': 'left', 'both': 'center'}

def _pad_series(s, side, fillchar):
    import pyarrow as pa
    if side not in _PAD_ALIGN:
        raise ValueError("Invalid side: side must be 'left', 'right', or 'both'")
    arr = pa.Array.from_pandas(s, type=pa.string())
    (padded,), _ = pad_columns([arr], align=_PAD_ALIGN[side], fillchar=fillchar)
    return pd.Series(padded.to_pylist(), index=s.index, name=s.name, dtype=s.dtype)

def pad(data, side='left', fillchar=' '):
    data = data.astype(str)
    if isinstance(data, pd.DataFrame):
        if isinstance(side, str):
            side = { x: side for x in data.columns }
        else:
            side = { x: side[x] if x in side else 'left' for x in data.columns }
        out = data.copy()
        for i, name in enumerate(data.columns):
            out.isetitem(i, _pad_series(data.iloc[:, i], side[name], fillchar))
        return out
    if isinstance(data, pd.Series):
        return _pad_series(data, side, fillchar)
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
                with open(colors) as fh:
                    colors = yaml.safe_load(fh)
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
        data = data.map(lambda x: ''.join([color_res(y, color_switch[style]) for y in x]))
    elif isinstance(data,pd.Series):
        data = data.map(lambda x: ''.join([color_res(y, color_switch[style]) for y in x]))

    return data
