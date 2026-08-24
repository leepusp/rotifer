"""
Access data hosted by UniProt.

This package groups the tools that query UniProt's web services.
Currently :mod:`rotifer.db.uniprot.webapi` wraps UniProt's REST API
for mapping identifiers between databases.

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

