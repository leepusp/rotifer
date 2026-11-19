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
