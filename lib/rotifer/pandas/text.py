"""Plain-text rendering of DataFrames: to_string, pad and to_color."""

import os
from itertools import chain
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

_ALIGNS = ('left', 'right', 'center')

def _check_align(align):
    if align not in _ALIGNS:
        raise ValueError(f"align must be 'left', 'right' or 'center', not {align!r}")

def _justify(name, width, align):
    """A header set the way pad_columns sets its values (str.center for 'center')."""
    return {'left': str.ljust, 'right': str.rjust, 'center': str.center}[align](name, width)

def _text_columns(df, na_rep, float_format, index):
    """Text of every cell of `df` as (names, pyarrow string arrays), one pair per
    column, the index first when asked. Columns are walked by position, so
    duplicate labels are fine. The one conversion shared by to_string and to_blocks."""
    import pyarrow as pa
    columns = []
    if index:
        columns.append((str(df.index.name or ""), df.index.to_series().astype("string")))
    for i, c in enumerate(df.columns):
        s = df.iloc[:, i]
        if float_format is not None and pd.api.types.is_float_dtype(s.dtype):
            body = s.map(lambda v: na_rep if pd.isna(v) else float_format % v).astype("string")
        else:
            body = s.astype("string")
        columns.append((str(c), body.fillna(na_rep)))
    names = [name for name, _ in columns]
    return names, [pa.Array.from_pandas(body, type=pa.string()) for _, body in columns]

def _join_rows(names, arrays, sep, align, taken=None):
    """Header and rows of converted columns, padded and joined. `taken` selects rows
    of every array first (an arrow take), so a block can be cut from a frame that
    was converted once."""
    import pyarrow as pa
    import pyarrow.compute as pc
    _check_align(align)
    if taken is not None:
        arrays = [a.take(taken) for a in arrays]
    padded, widths = pad_columns(arrays, names, align=align)
    header = sep.join(_justify(name, w, align) for name, w in zip(names, widths))
    if not padded:
        return header
    lines = pc.binary_join_element_wise(*padded, pa.scalar(sep))
    return header + "\n" + "\n".join(lines.to_pylist())

def to_string(df, sep=" : ", index=True, na_rep=" ", insert=None,
                      float_format=None, align="left"):
    _check_align(align)

    if insert:
        df = df.copy()
        for data in insert:
            if len(data) != 3:
                raise ValueError(f"insert entries must be (loc, column, value) tuples, got {data!r}")
            df.insert(loc=data[0], column=data[1], value=data[2])

    names, arrays = _text_columns(df, na_rep, float_format, index)
    return _join_rows(names, arrays, sep, align)

# pad's side is where the fill goes; pad_columns' align is where the text sits
_PAD_ALIGN = {'left': 'right', 'right': 'left', 'both': 'center'}

def _pad_series(s, side, fillchar, na_rep='', name=None):
    """Pad one Series; returns (padded Series, width). With name given, the width
    is at least len(name), so the column can be padded to fit its header."""
    import pyarrow as pa
    import pyarrow.compute as pc
    if side not in _PAD_ALIGN:
        raise ValueError("Invalid side: side must be 'left', 'right', or 'both'")
    arr = pc.fill_null(pa.Array.from_pandas(s, type=pa.string()), na_rep)
    (padded,), (width,) = pad_columns([arr], names=None if name is None else [name],
                                      align=_PAD_ALIGN[side], fillchar=fillchar)
    return pd.Series(padded.to_pylist(), index=s.index, name=s.name, dtype=s.dtype), width

def _justify_name(name, width, side, fillchar):
    """Set a column name on the same side as its values: 'left' pads on the left
    (name right-justified), 'right' on the right, 'both' centred."""
    if side == 'left':
        return name.rjust(width, fillchar)
    if side == 'right':
        return name.ljust(width, fillchar)
    return name.center(width, fillchar)

def pad(data, side='left', fillchar=' ', na_rep='', header=False):
    """Pad strings to a common width. side is where the fill goes ('left', 'right'
    or 'both'); on a DataFrame it may also be a dict of column -> side (missing
    columns default to 'left'). Missing values are first replaced by na_rep.
    With header=True (DataFrame only) the column names are padded too, each column
    to max(value width, name width), with the name on the same side as its values."""
    if not isinstance(data, (pd.DataFrame, pd.Series)):
        raise TypeError(f'pad expects a DataFrame or Series, not {type(data).__name__}')
    data = data.astype(str)
    if isinstance(data, pd.DataFrame):
        if isinstance(side, str):
            side = { x: side for x in data.columns }
        else:
            side = { x: side[x] if x in side else 'left' for x in data.columns }
        out = data
        names = []
        for i, name in enumerate(data.columns):
            s, width = _pad_series(data.iloc[:, i], side[name], fillchar, na_rep,
                                   name=str(name) if header else None)
            out.isetitem(i, s)
            names.append(_justify_name(str(name), width, side[name], fillchar) if header else name)
        if header:
            out.columns = names
        return out
    return _pad_series(data, side, fillchar, na_rep)[0]

def to_color(data, colors=None, style='fg', padding=None, fillchar=' '):
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

    if style not in ('fg', 'foreground', 'bg', 'background'):
        raise ValueError(f"style must be 'fg', 'foreground', 'bg' or 'background', not {style!r}")

    # Load colors
    default_color = "000"
    if colors is None or colors is False:
        colors = config["colors"]
    if not isinstance(colors,dict):
        if isinstance(colors,str):
            if os.path.exists(colors):
                with open(colors) as fh:
                    colors = yaml.safe_load(fh)
            else:
                default_color = colors
                colors = dict()

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

    if not isinstance(data, (pd.DataFrame, pd.Series)):
        return data

    # One escape-coded replacement per distinct character, built once per call;
    # characters not in colors get the default color.
    paint = color_bg if style in ('bg', 'background') else color_fg
    values = data.to_numpy().ravel()
    table = {ord(ch): paint(ch, colors[ch] if ch in colors else default_color)
             for ch in set(chain.from_iterable(values))}
    return data.map(lambda x: x.translate(table))
