# Import external modules
import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import datetime, timedelta

# Ete3
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa

# Rotifer
import rotifer
import rotifer.db.core
from rotifer import GlobalConfig
from rotifer.taxonomy.utils import lineage
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)
config = loadConfig(__name__, defaults = {
    "taxdump_file": None,
})

# Classes
class TaxonomyCursor(rotifer.db.core.BaseCursor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ete3 = NCBITaxa()
        self.taxcols = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']

    def update_database(self):
        """
        Update Ete3 taxonomy database.
        """
        from tempfile import TemporaryDirectory
        cwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            self.ete3.update_taxonomy_database()
            os.chdir(cwd)

    def getids(self,obj):
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            if "taxid" in o.columns:
                ids.update(o.taxid.astype(str))
            if "alternative_taxids" in o.columns:
                ids.update(o.alternative_taxids.astype(str))
        return ids

    def __getitem__(self, accessions):
        targets = self.parse_ids(accessions)
        targetsAsIntegers = pd.Series(list(targets)).astype(int)
        #logger.info(f'Loading {len(targets)} taxids from Ete3 local database')

        # Find replacements for phased-out taxids
        oldToNewDF = ",".join(targets)
        oldToNewDF = f'SELECT * from merged where taxid_old in ({oldToNewDF})'
        oldToNewDF = pd.read_sql(oldToNewDF, self.ete3.db)
        idmap = oldToNewDF.set_index('taxid_old').taxid_new.to_dict()
        targetsAsIntegers.replace(idmap, inplace=True)

        # Fetch data from database
        tl = self.ete3.get_lineage_translator(targetsAsIntegers.tolist())
        if len(tl) == 0:
            return pd.DataFrame(columns=self.taxcols)
        ti = pd.Series(tl.values()).explode().unique()
        ti = self.ete3.get_taxid_translator(ti)
        li = [[x,ti[x],np.nan,np.nan,"; ".join([ ti[y] for y in tl[x] if ti[y] not in ["root","cellular organisms"] ]),x] for x in tl ]
        li = pd.DataFrame(li, columns=self.taxcols)
        li['superkingdom'] = li.classification.str.split("; ", expand=True)[0]
        li['lineage'] = lineage(li.classification)

        # Shuffle alternative taxids and query accessions
        idmap = oldToNewDF.set_index('taxid_new').taxid_old.to_dict()
        li.alternative_taxids = li.alternative_taxids.replace(idmap)
        li = li.astype(str)

        # Register missing entries
        missing = targets - self.getids(li)
        if len(missing) > 0:
            # A local database does not change between attempts, but
            # another backend may well know the taxon
            self.update_missing(missing,"Accession not found in database",retry=False)

        #logger.info(f'Loaded {len(targets.intersection(self.getids(li)))} taxids from Ete3 database')
        return li

    def fetchone(self,accessions):
        # Process accessions
        li = self.__getitem__(accessions)
        for taxid, rows in li.groupby('taxid'):
            yield rows.copy()

    def fetchall(self,accessions):
        return self.__getitem__(accessions)
