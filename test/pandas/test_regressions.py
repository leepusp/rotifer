"""Regression tests for bugs fixed in rotifer.pandas; each test names the bug it pins."""
import pandas as pd
import pytest
import yaml
import rotifer.pandas as rp
from rotifer.pandas import functions as rpf
from rotifer.pandas import yaml_blocks as yb


def test_to_blocks_writes_to_path_buf(tmp_path):
    # Bug: to_blocks used os.PathLike without importing os -> NameError for a path buf.
    df = pd.DataFrame({'g': ['A', 'A', 'B'], 'pid': [1, 2, 3]})
    kwargs = dict(groupby='g', columns=['pid'], stats={'pid': 'nunique'}, sortrows={})
    target = tmp_path / 'out.txt'
    assert rp.to_blocks(df, buf=target, **kwargs) is None
    assert target.read_text() == rp.to_blocks(df, **kwargs)


def test_to_color_dataframe_uses_map_not_applymap():
    # Bug: to_color called DataFrame.applymap, removed in pandas 3 -> AttributeError.
    df = pd.DataFrame({'s': ['AC', 'GT']})
    assert rpf.to_color(df)['s'].tolist() == rpf.to_color(df['s']).tolist()


def test_pad_dataframe_string_side_applies_to_every_column():
    # Bug: a string side was tested with `x in side` (substring check), so a column
    # named 'e' or 'l' indexed the string and raised TypeError.
    df = pd.DataFrame({'e': ['a', 'bb'], 'z': ['ccc', 'd']})
    expected = pd.DataFrame({c: rpf.pad(df[c], side='right') for c in df.columns})
    pd.testing.assert_frame_equal(rpf.pad(df, side='right'), expected)


def test_to_color_loads_colors_yaml_safely(tmp_path):
    # Bug: colors YAML was read with yaml.Loader (unsafe) from an unclosed file handle.
    path = tmp_path / 'colors.yml'
    path.write_text('A: 196\n')
    assert '38;5;196' in rpf.to_color(pd.Series(['A']), colors=str(path)).iloc[0]

    unsafe = tmp_path / 'unsafe.yml'
    unsafe.write_text('A: !!python/object/apply:os.getcwd []\n')
    with pytest.raises(yaml.YAMLError):
        rpf.to_color(pd.Series(['A']), colors=str(unsafe))


def test_to_string_insert_rejects_malformed_entries():
    # Bug: to_string silently skipped insert entries whose length was not 3.
    df = pd.DataFrame({'a': [1]})
    with pytest.raises(ValueError, match="insert entries must be"):
        rp.to_string(df, insert=[(0, 'x')])


def test_to_yaml_rejects_merged_as_table_name():
    # Bug: the to_yaml docstring left `merged` out of the reserved names; it is reserved.
    df = pd.DataFrame({'g': ['A'], 'pid': [1]})
    with pytest.raises(ValueError):
        yb.to_yaml(df, groupby='g', merged={'pid': 'nunique'})


def test_pad_renders_missing_values_as_empty():
    # Bug: pad(pd.Series(['a', None])) raised TypeError (the width was NaN).
    assert rpf.pad(pd.Series(['a', None, 'ccc'])).tolist() == ['  a', '   ', 'ccc']
    assert rpf.pad(pd.Series(['a', None]), side='right').tolist() == ['a', ' ']
    out = rpf.pad(pd.DataFrame({'x': ['a', None], 'y': [1.5, float('nan')]}))
    assert out['x'].tolist() == ['a', ' ']
    assert out['y'].tolist() == ['1.5', '   ']
    assert rpf.pad(pd.Series([], dtype=str)).tolist() == []


def test_to_color_unknown_style_names_the_accepted_values():
    # Bug: an unknown style raised KeyError from an internal dict lookup.
    with pytest.raises(ValueError, match="style must be 'fg', 'foreground', 'bg' or 'background'"):
        rpf.to_color(pd.Series(['A']), style='underline')


def test_to_color_single_code_colours_every_character():
    # Bug: a single colour code (not a path) was never used; every character got 39;5;N
    # (default foreground), so no colour showed.
    out = rpf.to_color(pd.Series(['AC']), colors='196')
    assert out.tolist() == ['\x1b[38;5;196mA\x1b[m\x1b[38;5;196mC\x1b[m']
    out = rpf.to_color(pd.Series(['AC']), colors='196', style='bg')
    assert out.tolist() == ['\x1b[48;5;196mA\x1b[0m\x1b[48;5;196mC\x1b[0m']


def test_to_color_map_entry_equal_to_default_code_is_explicit():
    # Bug: a map entry whose code equalled the default code (000) was treated as
    # "not in the map" and written as 39;5;000 instead of 38;5;000.
    out = rpf.to_color(pd.Series(['AB']), colors={'A': '000'})
    assert out.tolist() == ['\x1b[38;5;000mA\x1b[m\x1b[39;5;000mB\x1b[m']


def test_to_color_renders_missing_values_as_empty_string():
    # Bug: a NaN or None value raised TypeError ('float' object is not iterable) without
    # padding. Without padding a missing value is ''; with padding it is filled to the
    # column width like any other empty value, coloured the same way as its padding.
    series = pd.Series(['A', None, float('nan')])
    out = rpf.to_color(series, colors={'A': '196'})
    assert out.tolist() == ['\x1b[38;5;196mA\x1b[m', '', '']
    frame = rpf.to_color(pd.DataFrame({'x': ['A', None]}), colors={'A': '196'})
    assert frame['x'].tolist() == ['\x1b[38;5;196mA\x1b[m', '']
    padded = rpf.to_color(pd.Series(['AB', None]), colors={'A': '196'}, padding='right')
    assert padded.tolist() == rpf.to_color(pd.Series(['AB', '  ']), colors={'A': '196'}).tolist()
    assert padded.iloc[1] == '\x1b[39;5;000m \x1b[m\x1b[39;5;000m \x1b[m'


def test_default_block_rows_warns_once_at_the_caller():
    # Bug: the internal default_header(df, ...) call emitted a second FutureWarning
    # pointing inside the library, so the caller saw two warnings, one of them misplaced.
    import warnings
    df = pd.DataFrame({'g': ['A', 'A', 'B'], 'v': [1, 2, 3]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        rp.default_block_rows(df, groupby='g', columns=['v'], stats={'v': 'sum'})
    warns = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert len(warns) == 1
    assert warns[0].filename == __file__
    assert '2026-11-30 (issue #TBD)' in str(warns[0].message)


def test_to_blocks_buf_path_round_trips_non_ascii_as_utf8(tmp_path):
    # Bug: _io opened files with the locale encoding, so non-ASCII text was not UTF-8.
    from rotifer.pandas._io import _read_source
    df = pd.DataFrame({'g': ['Ação', 'Ação', 'Pão'], 'nome': ['José', 'Zoë', 'Ångström']})
    kwargs = dict(groupby='g', columns=['nome'], stats={'nome': 'nunique'}, sortrows={})
    target = tmp_path / 'blocks.txt'
    assert rp.to_blocks(df, buf=target, **kwargs) is None
    expected = rp.to_blocks(df, **kwargs)
    assert target.read_bytes() == expected.encode('utf-8')
    assert _read_source(str(target)) == expected
