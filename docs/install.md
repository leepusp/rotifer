# Install

ROTIFER isn't packaged for pip or conda yet. You run it from a copy of
the repository, inside a conda environment that holds its dependencies
and the bioinformatics programs it calls. Setting it up takes five
steps, and each command below can be copied with the button in the
corner of its block.

You need [conda](https://docs.conda.io/) (or mamba) and git. The
commands are written for a Unix shell such as bash or zsh.

## 1. Clone the repository

Download the source code and move into it:

```bash
git clone https://github.com/leepusp/rotifer.git
cd rotifer
```

The folder you just created is the path that the later steps call
`/path_to_git/rotifer`. Run `pwd` inside it to print the full path, and
keep it at hand.

## 2. Create the conda environment

This one command creates an environment called `rotifer` with
everything ROTIFER uses: its libraries and command-line programs such
as HMMER, HH-suite, MMseqs2, Foldseek, MAFFT, FAMSA and BLAST. It can
take a while.

```bash
conda create --name rotifer biopython pandas numpy termcolor PyYAML tqdm \
  ipython ascii_graph matplotlib ete3 prettytable seaborn pygraphviz \
  bcbio-gff dash-bio gspread leidenalg igraph networkx oauth2client sqlalchemy \
  argcomplete perl-bioperl perl-moose perl-ipc-system-simple \
  perl-config-any perl-appconfig perl-carp-clan \
  perl-file-which perl-module-list perl-statistics-descriptive perl-yaml \
  hhsuite hmmer mmseqs2 foldseek mafft famsa blast parallel datamash
conda activate rotifer
pip install dash-cytoscape
cpan -i Module::Find
```

The last three lines finish the job from inside the new environment:
`dash-cytoscape` and `Module::Find` aren't in the conda channels, so
they are installed separately.

## 3. Activate the environment

Every new terminal starts outside the environment. Activate it whenever
you want to use ROTIFER:

```bash
conda activate rotifer
```

## 4. Prepare the environment's activation scripts

Conda runs the scripts in `activate.d` each time the environment is
activated, and those in `deactivate.d` when it is deactivated. These
commands create one empty script in each folder, then register ROTIFER's
Python library with the environment:

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
conda-develop /path_to_git/rotifer/lib
```

Replace `/path_to_git/rotifer` with the path from step 1.
`$CONDA_PREFIX` needs no change: conda sets it to the active
environment's folder, which is why step 3 comes first.

`conda-develop` is part of the `conda-build` package. If your shell says
the command doesn't exist, install it with `conda install conda-build`
and run the line again.

## 5. Put ROTIFER on your paths

The activation script from step 4 is still empty. It needs two lines:
one adds ROTIFER's command-line tools (`rotifer/bin`) to `PATH`, the
other adds its Python library (`rotifer/lib`) to `PYTHONPATH`:

```bash
export PATH="/path_to_git/rotifer/bin:$PATH"
export PYTHONPATH="$PYTHONPATH:/path_to_git/rotifer/lib"
```

To write them into the script in one go, copy the block below, replace
`/path_to_git/rotifer` with your path in both lines, and run it:

```bash
cat >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh <<'EOF'
export PATH="/path_to_git/rotifer/bin:$PATH"
export PYTHONPATH="$PYTHONPATH:/path_to_git/rotifer/lib"
EOF
```

The quotes around `'EOF'` keep `$PATH` and `$PYTHONPATH` as written, so
they are expanded each time the environment is activated, not once
now.

## Check the installation

Deactivate and activate the environment so the new script runs, then
ask Python where it finds ROTIFER:

```bash
conda deactivate
conda activate rotifer
python -c "import rotifer; print(rotifer.__file__)"
```

The path printed should point inside your `rotifer/lib` folder. If it
does, you're ready: the {doc}`documentation overview <overview>` is a
good place to start.
