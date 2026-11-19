# Block reports

A block report shows the rows of a table in groups. Each group is a block: a
header line that names the group, optionally followed by a few statistics about
it, and then the rows that belong to the group. Take a table of proteins grouped
by sequence cluster `c80e3`. Every cluster gets its own block, so the proteins of
one cluster are read together.

ROTIFER writes block reports in two formats. Both are planned by the same code, so
the blocks, their order, the order of the rows and any sampling are the same in
either format.

- Text, written by {func}`rotifer.pandas.to_blocks`. Columns are aligned, which
  makes it a quick read in a terminal or a pager. It is not meant to be parsed
  back.
- YAML, written by {func}`rotifer.pandas.to_yaml`. The file is hand-editable, so
  you can correct values, delete rows or add notes. {func}`rotifer.pandas.yaml_blocks.read_yaml`
  turns it back into DataFrames.

Use the text format to look at a result. Use the YAML format to change the report
by hand or to load the edited report into Python.

The examples on this page use this table of four proteins:

```python
import pandas as pd
import rotifer.pandas as rp

df = pd.DataFrame({
    'c80e3': ['A', 'A', 'A', 'B'],
    'pid': ['p1', 'p2', 'p3', 'p4'],
    'lineage': ['Bacteria;Firmicutes', 'Bacteria;Firmicutes',
                'Bacteria;Proteobacteria', 'Archaea;Euryarchaeota'],
    'plen': [310, 95, 1200, 45],
    'product': ['ABC transporter|ATPase', 'hypothetical protein', 'kinase',
                'small protein'],
})
```

## Text format

```python
print(rp.to_blocks(df, 'c80e3', stats={'lineage': 'nunique'}))
```

```text
#
# c80e3=A; lineage=2
  : pid : lineage                 : plen : product
0 : p1  : Bacteria;Firmicutes     : 310  : ABC transporter|ATPase
1 : p2  : Bacteria;Firmicutes     : 95   : hypothetical protein
2 : p3  : Bacteria;Proteobacteria : 1200 : kinase
#
# c80e3=B; lineage=1
  : pid : lineage               : plen : product
3 : p4  : Archaea;Euryarchaeota : 45   : small protein
```

Each block starts with a line holding `#`, then a line with the group and the
statistics. In `# c80e3=A; lineage=2`, `c80e3=A` is the group and `lineage=2` says
that the group holds two distinct lineages. The rows follow. The first column
is the index of the row in the original table. The line under the block header
gives the column names, and the cells are separated by ` : `.

The main keyword arguments are:

- `columns`: the columns shown in each row. All columns except the grouping
  columns when omitted.
- `stats`: per-block statistics for the header, as `{column: aggregation}` or
  `{name: (column, aggregation)}`.
- `sortrows` and `sortblocks`: the order of the rows inside a block and the order
  of the blocks.
- `sample`: the number of rows kept per block, chosen at random.

The full list is in the reference of {func}`rotifer.pandas.to_blocks`.

## YAML format

```python
print(rp.to_yaml(df, 'c80e3', stats={'lineage': 'nunique'}))
```

```yaml
- parameters:
    groupby: [c80e3]

# -----
- key: [A]
  stats: |
    | lineage |
    | ------- |
    | 2       |
  rows: |
    | pid | lineage                 | plen | product                 |
    | --- | ----------------------- | ---- | ----------------------- |
    | p1  | Bacteria;Firmicutes     | 310  | ABC transporter\|ATPase |
    | p2  | Bacteria;Firmicutes     | 95   | hypothetical protein    |
    | p3  | Bacteria;Proteobacteria | 1200 | kinase                  |

# -----
- key: [B]
  stats: |
    | lineage |
    | ------- |
    | 1       |
  rows: |
    | pid | lineage               | plen | product       |
    | --- | --------------------- | ---- | ------------- |
    | p4  | Archaea;Euryarchaeota | 45   | small protein |
```

The document is a YAML sequence with these elements:

- The first element, `parameters`, holds the settings of the whole file. Its
  `groupby` lists the grouping columns. They are stored once here and not repeated
  in every block.
- Each following element is one block. Its `key` gives the values of the group,
  one for each grouping column in the order of `groupby`.
- Every table of the block is a literal block scalar (the `|` after the name)
  holding a markdown table. The `stats` table holds the statistics of the block.
  The `rows` table holds the rows of the block. The names come from the keywords
  given to `to_yaml`, and the main table is called `rows` unless `name` says
  otherwise. The names `key`, `parameters` and `merged` are reserved and cannot be
  used for a table.
- A literal pipe inside a cell is written `\|`. Newlines inside a cell are
  replaced by a space.
- The padding of the tables and the spaces around the cells are for reading only.
  When the document is read back, the spaces are removed, so a table that you
  edit by hand with uneven pipes still reads correctly.
- Blocks are separated by a blank line and a `# -----` comment line.

Tables that are computed for each block, such as a table of the most frequent
domain, are given as functions. `to_yaml` calls each function with the block as a
DataFrame and writes the table it returns. This is slower than the statistics
written from a mapping, because a DataFrame is built for each block.

## Reading a YAML document

{func}`rotifer.pandas.yaml_blocks.read_yaml` takes a path, an open file or the
text of a document. The grouping columns are read from `parameters`, and their
values are taken from each `key`, so they do not need to be stored in the tables.

```python
text = rp.to_yaml(df, 'c80e3', stats={'lineage': 'nunique'},
                  columns=['pid', 'product'])
print(rp.read_yaml(text))
```

```text
  c80e3 pid                 product  lineage
0     A  p1  ABC transporter|ATPase        2
1     A  p2    hypothetical protein        2
2     A  p3                  kinase        2
3     B  p4           small protein        1
```

By default, `read_yaml` works in two steps. It concatenates the tables of the
same name across all blocks, so the `stats` of every block become one `stats`
table and the `rows` of every block become one `rows` table. It then merges the
tables on the grouping columns, which is why each row above has its `lineage`
count. When a column other than the join columns occurs in more than one table,
the copy from a later table is renamed with the table's name as a suffix. In the
example above, a `lineage` column in both tables would become `lineage_stats`.

The `merge`, `drop`, `keep` and `how` arguments change this. With `merge=False`,
the function returns a dictionary with one frame per table:

```python
tables = rp.read_yaml(text, merge=False)
print(tables['stats'])
```

```text
  c80e3  lineage
0     A        2
1     B        1
```

Merging only some tables returns a dictionary as well, with the merged frame
under the key `merged` and the tables that were not merged under their own names.
The full rules for `merge`, `drop` and `keep` are in the reference of
{func}`rotifer.pandas.yaml_blocks.read_yaml`.

## Gene neighborhoods

{class}`rotifer.genome.data.NeighborhoodDF` has two shortcuts, `to_blocks` and
`to_yaml`. They call the functions above with presets from
{mod}`rotifer.genome.blocks`, so a neighborhood table can be written with one
call. The presets are plain dictionaries:

- {data}`rotifer.genome.blocks.NEIGHBORHOOD_BLOCKS` groups by `c80e3`. It shows the
  protein id, the 70% and 100% clusters, the length, the Pfam and Aravind domains,
  the product and the taxonomy. Each block header holds the number of distinct
  lineages, neighborhoods, 70% clusters and proteins. Rows are sorted by lineage,
  classification and organism.
- {data}`rotifer.genome.blocks.NEIGHBORHOOD_YAML` groups by `c80e3` and shows the
  same columns. Rows are sorted by taxonomy, then by the 70% cluster. Blocks are
  sorted by the number of phyla and of 70% clusters, most first. Each block has an
  `annotation` table with the most frequent Pfam domain, a `stats` table and a
  `dist` table with the number of 70% clusters for each domain.

A keyword given to the method replaces the preset entry with the same name, and
the other entries stay as they are. The presets need the columns they name. A
missing column raises `KeyError` when it is a grouping or statistics column, and
is ignored when it is a column to show or to sort by.

The example below uses a `NeighborhoodDF` with the columns of the presets. Only
the first block of the YAML output is shown.

```python
from rotifer.genome.data import NeighborhoodDF

nei = NeighborhoodDF({
    'c80e3': ['A', 'A', 'B'],
    'c80i70': ['x1', 'x2', 'x3'],
    'c100i100': ['y1', 'y2', 'y3'],
    'pid': ['p1', 'p2', 'p3'],
    'plen': [300, 120, 80],
    'pfam': ['PF00001', 'PF00001', 'PF00002'],
    'aravind': ['', '', ''],
    'product': ['ABC transporter', 'protein', 'kinase'],
    'organism': ['Org one', 'Org two', 'Org three'],
    'lineage': ['Bacteria;Firmicutes', 'Bacteria;Proteobacteria',
                'Archaea;Euryarchaeota'],
    'classification': ['c1', 'c2', 'c3'],
    'block_id': ['b1', 'b1', 'b2'],
    'superkingdom': ['Bacteria', 'Bacteria', 'Archaea'],
    'phylum': ['Firmicutes', 'Proteobacteria', 'Euryarchaeota'],
    'taxid': [1, 2, 3],
})
print(nei.to_blocks())
print(nei.to_yaml())
```

```text
#
# c80e3=A; lineage=2; block_id=1; c80i70=2; pid=2
  : pid : c80i70 : c100i100 : plen : pfam    : aravind : product         : organism : lineage                 : classification
0 : p1  : x1     : y1       : 300  : PF00001 :         : ABC transporter : Org one  : Bacteria;Firmicutes     : c1
1 : p2  : x2     : y2       : 120  : PF00001 :         : protein         : Org two  : Bacteria;Proteobacteria : c2
#
# c80e3=B; lineage=1; block_id=1; c80i70=1; pid=1
  : pid : c80i70 : c100i100 : plen : pfam    : aravind : product : organism  : lineage               : classification
2 : p3  : x3     : y3       : 80   : PF00002 :         : kinase  : Org three : Archaea;Euryarchaeota : c3
```

The YAML output starts like this, with the first block:

```yaml
- parameters:
    groupby: [c80e3]

# -----
- key: [A]
  annotation: |
    | pfam    | curated |
    | ------- | ------- |
    | PF00001 | PF00001 |
  stats: |
    | superkingdom | lineage | phylum | taxid | c80i70 | pid |
    | ------------ | ------- | ------ | ----- | ------ | --- |
    | 1            | 2       | 2      | 2     | 2      | 2   |
  dist: |
    | pfam    | c80i70 |
    | ------- | ------ |
    | PF00001 | 2      |
  rows: |
    | pid | c80i70 | c100i100 | plen | pfam    | aravind | product         | organism | lineage                 | classification |
    | --- | ------ | -------- | ---- | ------- | ------- | --------------- | -------- | ----------------------- | -------------- |
    | p1  | x1     | y1       | 300  | PF00001 |         | ABC transporter | Org one  | Bacteria;Firmicutes     | c1             |
    | p2  | x2     | y2       | 120  | PF00001 |         | protein         | Org two  | Bacteria;Proteobacteria | c2             |
```

To change one setting for a single call, pass it as a keyword. This writes the
same blocks with two columns and no statistics:

```python
print(nei.to_blocks(columns=['pid', 'product'], stats=None))
```

To reuse a preset with other choices, copy it and edit the copy:

```python
from rotifer.genome.blocks import NEIGHBORHOOD_BLOCKS
import rotifer.pandas as rp

options = dict(NEIGHBORHOOD_BLOCKS, groupby='c100i100')
print(rp.to_blocks(nei, **options))
```

The reference of the methods is in {class}`rotifer.genome.data.NeighborhoodDF`,
and the reference of the presets is in {mod}`rotifer.genome.blocks`.
