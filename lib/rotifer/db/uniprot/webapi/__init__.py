"""
Query UniProt's REST web services.

:mod:`rotifer.db.uniprot.webapi.idmapping` implements the identifier
mapping workflow: submit a batch of accessions, poll until the job
finishes and download the results.

Note
----
The ``local_database_path`` configuration default below is not used
anywhere in this package; see ``docs/OPEN_QUESTIONS.md``.
"""

import os
import rotifer
from rotifer import GlobalConfig
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

# Configuration
config = loadConfig(__name__.replace('rotifer.',':'), defaults = {
    'local_database_path': os.path.join(GlobalConfig['data'],"fadb","nr","nr"),
})

# FUNCTIONS

