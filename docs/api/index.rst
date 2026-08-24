API reference
=============

This reference is generated from the docstrings of the installed
``rotifer`` package. The :doc:`data access section </db/index>` is the
curated entry point for the ``rotifer.db`` subpackage documented below.

Modules that cannot currently be imported are not listed here; see the
``docs/OPEN_QUESTIONS.md`` file in the repository for the details.

Data access (rotifer.db)
------------------------

.. autosummary::
   :toctree: generated

   rotifer.db
   rotifer.db.core
   rotifer.db.methods
   rotifer.db.delegator
   rotifer.db.parallel
   rotifer.db.cli
   rotifer.db.ncbi
   rotifer.db.ncbi.entrez
   rotifer.db.ncbi.ftp
   rotifer.db.ncbi.mirror
   rotifer.db.ncbi.utils
   rotifer.db.local
   rotifer.db.local.easel
   rotifer.db.local.ete3
   rotifer.db.sql
   rotifer.db.sql.sqlite3
   rotifer.db.uniprot
   rotifer.db.uniprot.webapi
   rotifer.db.uniprot.webapi.idmapping

Core infrastructure
-------------------

.. autosummary::
   :toctree: generated

   rotifer
   rotifer.core
   rotifer.core.cli
   rotifer.core.functions
   rotifer.core.io.ftp
   rotifer.core.loadpath
   rotifer.core.log
   rotifer.core.logger
   rotifer.core.methods
   rotifer.core.time

Genomes and annotation
----------------------

.. autosummary::
   :toctree: generated

   rotifer.genome
   rotifer.genome.data
   rotifer.genome.database
   rotifer.genome.db.clickhouse
   rotifer.genome.db.postgres
   rotifer.genome.io
   rotifer.genome.utils

Sequence and interval utilities
-------------------------------

.. autosummary::
   :toctree: generated

   rotifer.interval.utils
   rotifer.io.dali
   rotifer.io.fileinput
   rotifer.io.hhsuite

Data manipulation and taxonomy
------------------------------

.. autosummary::
   :toctree: generated

   rotifer.pandas
   rotifer.pandas.functions
   rotifer.taxonomy
   rotifer.taxonomy.utils

Other subpackages
-----------------

.. autosummary::
   :toctree: generated

   rotifer.alchemy
   rotifer.alchemy.connect
   rotifer.cluster
   rotifer.cluster.cluster
   rotifer.pipeline
   rotifer.view
   rotifer.view.functions
