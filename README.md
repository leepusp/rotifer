# ROTIFER

## Rapid Open-source Tools and Infrastructure For data Exploration and Research

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.

Documentation: https://leepusp.github.io/rotifer/

## Installation

ROTIFER needs Python 3.10 or newer.

### Python library only

```Bash
pip install "rotifer[all] @ git+https://github.com/leepusp/rotifer"
```

`rotifer[all]` installs every optional dependency. To keep the install small,
pick only the groups you need instead (e.g. `rotifer[seq,db]`):

| extra | for |
|-------|-----|
| `seq` | `rotifer.seq`, `rotifer.genome` |
| `db` | `rotifer.db`, `rotifer.alchemy`, `rotifer.tools` |
| `viz` | `rotifer.cluster` and plots |
| `network` | neighborhood graphs in `rotifer.genome` |
| `cli` | the Python tools in `bin/` |

### Full toolset (library, command line tools and their programs)

The tools in `bin/` call external programs (BLAST, HMMER, MMseqs2, ...) and
Perl modules that pip cannot install. `environment.yml` installs them all with
conda, plus the library in editable mode:

```Bash
git clone https://github.com/leepusp/rotifer.git
cd rotifer
conda env create -f environment.yml
conda activate rotifer
cpanm Module::Find
```

Then put the tools on your `PATH`, e.g. in
`$CONDA_PREFIX/etc/conda/activate.d/rotifer.sh`:

```Bash
export PATH="/path/to/rotifer/bin:$PATH"
```

## Examples of use

### colecting igps of huge number os pids:

```python
from rotifer.db.ncbi import entrez
ic = entrez.IPGCursor(progress=True, batch_size=200)
for df in ic.fetchone([query_list]):
    if os.path.exists("ipg.tsv"):
        df.to_csv("ipg.tsv", sep="\t", index=False, mode="a", header=False)
    else:
        df.to_csv("ipg.tsv", sep="\t", index=False)
ipg = pd.read_csv('ipg.tsv',sep="\t")
```

Remember to check ic.missing to see if there are any missing ipg, you could further collect it with:
SequenceCursor + ic.\_seqrecords_to_ipg

### Collect the assembly reports files:

```python
from rotifer.db import ncbi
assembly_reports = ncbi.assemblies(taxonomy=True)
```

### Merging ipgs with assembly_reports to collect taxonomy of eaach pid:

assembly_reports.query('assembly in @ipg.assembly.to_list()').merge(ipg, on='assembly').pid.nunique()

### Activating debbug function

```python
import rotifer
rotifer.logger.setLevel(rotifer.logging.DEBUG)
```
