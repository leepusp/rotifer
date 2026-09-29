# ROTIFER

## Rapid Open-source Tools and Infrastructure For data Exploration and Research

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.

Documentation: https://leepusp.github.io/rotifer/

Installing ROTIFER package with conda enviroments:

1. Clone this git reposytory.
2. Create a rotifer's conda enviroment:

```Bash
conda create --name rotifer biopython pandas numpy termcolor PyYAML tqdm \
  ipython ascii_graph matplotlib ete3 prettytable seaborn pygraphviz \
  bcbio-gff dash-bio gspread leidenalg igraph networkx oauth2client sqlalchemy \
  argcomplete perl-bioperl perl-moose perl-ipc-system-simple \
  perl-config-any perl-appconfig perl-carp-clan \
  perl-file-which perl-module-list perl-statistics-descriptive perl-yaml \
  hhsuite hmmer mmseqs2 foldseek mafft famsa blast parallel datamash
conda activate rotifer
pip install dash-cytoscape
cpan i Module::Find
```

3. Activating the rotifer enviroment:

```Bash
conda activate rotifer
```

4. Setting up the rotifer tools at the conda rotifer enviroment:

```Bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
conda-develop /path_to_git/rotifer/lib

```

5. Write the following path to eh activate.d file:

```Bash
  export PATH="path_to_git/rotifer/bin:$PATH"
  export PYTHONPATH="$PYTHONPATH:/path_to_git/rotifer/lib"
```

`

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
