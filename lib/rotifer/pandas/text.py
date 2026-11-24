"""
Plain-text rendering of DataFrames and Series.

:func:`to_string` writes a frame as aligned text with a separator between
columns, :func:`pad` pads strings to a common width and :func:`to_color`
adds terminal colour codes to each character. :func:`pad_columns` is the
padding routine they share. :func:`rotifer.pandas.blocks.to_blocks` builds on
the same column conversion. The module reads its configuration under the key
``rotifer.pandas.functions``, which holds the default colour map of
:func:`to_color`.
"""

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
    """
    Pad pyarrow string arrays to one common width per column.

    Parameters
    ----------
    arrays : list of pyarrow.Array
        One string array per column. Nulls are padded as empty strings.
    names : list of str, optional
        Column names. When given, the width of a column is at least the length
        of its name, so a header fits above the values.
    align : {'left', 'right', 'center'}, default 'left'
        Where the text sits in the field: 'left' pads on the right, 'right'
        pads on the left and 'center' pads on both sides (an odd amount of
        fill goes where :meth:`str.center` puts it).
    fillchar : str, default ' '
        Single character used for padding.

    Returns
    -------
    padded : list of pyarrow.Array
        The padded arrays, in the order given.
    widths : list of int
        Width of each column: its longest string, or the length of its name
        when that is longer. An empty or all-null column without a name has
        width 0.

    Raises
    ------
    ValueError
        If ``align`` is not one of the three values above.

    Examples
    --------
    >>> import pyarrow as pa
    >>> from rotifer.pandas.text import pad_columns
    >>> padded, widths = pad_columns([pa.array(['a', 'bbb', None])],
    ...                              names=['name'], align='right')
    >>> padded[0].to_pylist()
    ['   a', ' bbb', '    ']
    >>> widths
    [4]
    """
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
    """Raise ValueError unless align is 'left', 'right' or 'center'."""
    if align not in _ALIGNS:
        raise ValueError(f"align must be 'left', 'right' or 'center', not {align!r}")

def _justify(name, width, align):
    """A header set the way pad_columns sets its values (str.center for 'center')."""
    return {'left': str.ljust, 'right': str.rjust, 'center': str.center}[align](name, width)

def _text_columns(df, na_rep, float_format, index):
    """Text of every cell of `df` as (names, pyarrow string arrays), one pair per
    column, the index first when asked. Columns are walked by position, so
    duplicate labels work. This conversion is shared by to_string and to_blocks."""
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
    """
    Render a DataFrame as aligned plain text.

    The first line holds the column names and each following line one row.
    Cells are padded to the width of their column and joined with ``sep``.

    Parameters
    ----------
    df : pandas.DataFrame
        Frame to render. It is not modified.
    sep : str, default ' : '
        Text placed between columns.
    index : bool, default True
        Whether to write the index as the first column. Its header is the
        name of the index, or empty when it has none.
    na_rep : str, default ' '
        Text written for missing values.
    insert : list of tuple, optional
        Columns to add before rendering, as ``(loc, column, value)`` tuples
        with the meaning of :meth:`pandas.DataFrame.insert`. They are added in
        order to a copy of ``df``, so later entries see the earlier ones.
    float_format : str, optional
        A printf-style format such as ``'%.2f'``, applied to the values of
        float columns. By default floats are written with their usual text.
    align : {'left', 'right', 'center'}, default 'left'
        Alignment of headers and values within each column.

    Returns
    -------
    str
        The text, without a final newline. A frame with no columns and no
        index yields the empty header line only.

    Raises
    ------
    ValueError
        If ``align`` is not one of the three values above, or an ``insert``
        entry does not have exactly three items.

    See Also
    --------
    pad : Pad the strings of a frame or series to a common width.
    rotifer.pandas.md_table : Write a frame as a Markdown table.
    rotifer.pandas.to_blocks : Write one text block per group of rows.

    Notes
    -----
    Headers are aligned like the values of their column, and a column is as
    wide as its longest value or its name. Missing values become ``na_rep``
    in every column. With ``float_format``, missing values in float columns are
    not formatted but replaced by ``na_rep``, so ``na_rep`` need not be a
    number. Columns are handled by position, so repeated column names work.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas import to_string
    >>> df = pd.DataFrame({'gene': ['a', 'bb'], 'score': [1.5, None]})
    >>> print(to_string(df, float_format='%.2f', na_rep='NA', align='right'))
      : gene : score
    0 :    a :  1.50
    1 :   bb :    NA
    """
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
    """
    Pad strings to a common width.

    Every value is converted to text first. Padding the left side of a value
    right-justifies it, and padding the right side left-justifies it.

    Parameters
    ----------
    data : pandas.Series or pandas.DataFrame
        Values to pad. Each column of a frame is padded to its own width.
    side : {'left', 'right', 'both'} or dict, default 'left'
        Where the fill goes. ``'left'`` puts the fill on the left, so the text
        is right-justified (``'  a'``), ``'right'`` puts it on the right, so
        the text is left-justified (``'a  '``), and ``'both'`` centres the
        text. For a DataFrame, a dict maps column names to one of these
        values, and columns missing from it use ``'left'``.
    fillchar : str, default ' '
        Single character used for padding.
    na_rep : str, default ''
        Text that replaces missing values before padding.
    header : bool, default False
        For a DataFrame, also pad the column names. Each column is then as wide
        as its longest value or its name, and the name sits on the same side
        as the values. Ignored for a Series.

    Returns
    -------
    pandas.Series or pandas.DataFrame
        Same type and index as ``data``, with string values. A DataFrame
        keeps its columns, relabelled with the padded names when
        ``header=True``.

    Raises
    ------
    TypeError
        If ``data`` is neither a DataFrame nor a Series.
    ValueError
        If a side is not ``'left'``, ``'right'`` or ``'both'``.

    See Also
    --------
    to_string : Render a whole frame as aligned text.
    rotifer.pandas.md_table : Write a frame as a Markdown table.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas import pad
    >>> pad(pd.Series(['a', 'bbb'])).tolist()
    ['  a', 'bbb']
    >>> pad(pd.Series(['a', 'bbb']), side='right', fillchar='.').tolist()
    ['a..', 'bbb']
    >>> pad(pd.DataFrame({'id': [1, 22]}), header=True)
       id
    0   1
    1  22
    """
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
    Add terminal colour codes to each character of a frame or series.

    Every value becomes a string in which each character is wrapped in an
    ANSI escape sequence. Foreground colours use ``\\x1b[38;5;Nm`` and
    background colours use ``\\x1b[48;5;Nm``, where N is a colour of the
    256-colour palette, and the character is followed by a reset sequence.

    Parameters
    ----------
    data : pandas.DataFrame or pandas.Series
        Strings to colour. Anything else is returned unchanged (after
        padding, which also needs a frame or series).
    colors : dict, str or path-like, optional
        Colour of each character. It can be a dict from character to colour
        code, the path of a YAML file holding such a mapping, or a string that
        is not an existing path, taken as a single colour code used for every
        character. By default (None or False) the map in the ``colors`` entry
        of the configuration is used, which is the data file
        ``colors/aminoacids.yml`` and colours amino acid letters.
    style : {'fg', 'foreground', 'bg', 'background'}, default 'fg'
        Colour the characters (``'fg'``) or their background (``'bg'``).
    padding : {'left', 'right', 'both'}, optional
        If given, pad the values with :func:`pad` first so that all values of
        a column have the same width.
    fillchar : str, default ' '
        Character used for padding.

    Returns
    -------
    pandas.DataFrame or pandas.Series
        Same shape and index as ``data``, holding the coloured strings. A
        missing value (NaN or None) is an empty string, padded to the column
        width when ``padding`` is given.

    Raises
    ------
    ValueError
        If ``style`` is not one of the four values above.

    See Also
    --------
    pad : Pad strings to a common width.

    Notes
    -----
    A single colour code given as ``colors`` colours every character with that
    code (``38;5;N`` or ``48;5;N``). A character that is missing from a colour
    map gets the terminal default colour instead: it is written with the
    sequence ``39;5;000`` (or ``49;5;000`` for a background), not ``38;5;`` or
    ``48;5;``, so a terminal shows its own default colour or ignores the code.
    A map entry whose code is ``000`` is an explicit colour and is written
    with ``38;5;000``.

    Examples
    --------
    >>> import pandas as pd
    >>> from rotifer.pandas import to_color
    >>> print(repr(to_color(pd.Series(['A']), colors={'A': '196'}).iloc[0]))
    '\\x1b[38;5;196mA\\x1b[m'
    >>> print(repr(to_color(pd.Series(['A']), colors={'A': '196'}, style='bg').iloc[0]))
    '\\x1b[48;5;196mA\\x1b[0m'
    """
    import yaml

    if style not in ('fg', 'foreground', 'bg', 'background'):
        raise ValueError(f"style must be 'fg', 'foreground', 'bg' or 'background', not {style!r}")

    # Load colors
    default_color = "000"
    single_code = None
    if colors is None or colors is False:
        colors = config["colors"]
    if not isinstance(colors,dict):
        if isinstance(colors,str):
            if os.path.exists(colors):
                with open(colors) as fh:
                    colors = yaml.safe_load(fh)
            else:
                single_code = colors
                colors = dict()

    def color_of(ch):
        # (code, explicit): explicit characters get 38/48;5;N, the others the default
        if single_code is not None:
            return single_code, True
        if ch in colors:
            return colors[ch], True
        return default_color, False

    def color_bg(s, color, explicit):
        if explicit:
            return f'\033[48;5;{color}m{s}\033[0m'
        return f'\033[49;5;{default_color}m{s}\033[0m'

    def color_fg(s, color, explicit):
        if explicit:
            return f'\033[38;5;{color}m{s}\033[m'
        return f'\033[39;5;{default_color}m{s}\033[m'

    # Missing values (NaN or None) become '' first, so that padding fills them to the
    # column width like any other empty value.
    if isinstance(data, (pd.DataFrame, pd.Series)):
        data = data.mask(data.isna(), '')

    # Forcing all data to be of the same length
    if padding:
        data = pad(data, side=padding, fillchar=fillchar)

    if not isinstance(data, (pd.DataFrame, pd.Series)):
        return data

    # One escape-coded replacement per distinct character, built once per call.
    paint = color_bg if style in ('bg', 'background') else color_fg
    values = data.to_numpy().ravel()
    table = {ord(ch): paint(ch, *color_of(ch))
             for ch in set(chain.from_iterable(values))}
    return data.map(lambda x: x.translate(table))
