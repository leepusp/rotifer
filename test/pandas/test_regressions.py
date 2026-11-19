"""Regression tests for bugs fixed in rotifer.pandas; each test names the bug it pins."""
import pandas as pd
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
